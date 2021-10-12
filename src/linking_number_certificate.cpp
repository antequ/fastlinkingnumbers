/**
 * Copyright (C) 2021 Ante Qu <antequ@cs.stanford.edu>.
 * Copyright (C) 2018 Jonathan Leaf.  All rights reserved.
 * Originally "sim_model.cpp" by Jonathan Leaf. Modified by Ante Qu.
 */

#include "linking_number_certificate.h"

#include "compute_link.h"
#include "discretize.h"
#include "model.h"
#include "curve.h" // for SegmentType
#include "potential_link_search.h"

#include <fstream>
#include <iostream>
#include <string>
#include <unordered_set>
#include <vector>

namespace VerifyCurves {
// Parallelization parameters
constexpr int NCURVES_PARALLEL_THRESHOLD = 6;
// Barnes–Hut parameters
constexpr double CERTIFICATE_EPSSQ = 1e-24;
constexpr double TARGET_ERROR = 0.1; // This is the Barnes–Hut target error.
constexpr double ERROR_EST_SCALE = 0.05;

/**
 * This constructor constructs a blank certificate.
 */
LinkingNumberCertificate::LinkingNumberCertificate() {}

/**
 * This constructor loads an already-computed linking-number certificate from
 * an input file.
 *
 * @param [in] filename is the name of the file with an existing
 *     linking-number certificate.
 *
 * @post the certificate_matrix_ matrix is created with the contents of the
 *     input file.
 */
LinkingNumberCertificate::LinkingNumberCertificate(
    const std::string &filename) {
  ImportFromFile(filename);
}

/// This calls the functions in potential_link_search.h to generate curve
/// bounding boxes and perform a potential link search and generate a list of
/// potential links suitable for discretization. This means every loop pair is
/// duplicated, e.g. (3,4) is also stored as (4,3), so that both curves 3 and
/// 4 will be checked against each other in the discretization passes.
///
/// @param [inout] model is the input curve collection.
///
/// @returns the output list of potential links suitable for discretization.
///
/// @post The curves in model have their curve bounding boxes precomputed and
///     cached.
std::vector<std::vector<int>> PotentialLinkSearch(const Model &model) {
  PrecomputeCurveBoundingBoxes(model.GetCurves());
  return GetPotentialLinksListPerCurve(model.GetCurves());
}

/// This calls the functions in discretize.h and potential_link_search.h to
/// discretize the curves and get a new, unique (each potential link pair only
/// appears once), list of potential links suitable for linking-number
/// computation.
///
/// @param [in] model is the input curve collection.
/// @param [in] potential_links_list is the input vector-of-vectors list of
///     potential links from PotentialLinkSearch.
/// @param [out] out_discretized_curves is the collection of discretized curves.
/// @param [out] out_potential_links_unique_per_curve is the output list of
///     unique potential links suitable for linking-number computation.
///
/// @post The discretized curves have their curve bounding boxes precomputed and
///     cached.
void DiscretizeAndGetNewPotentialLinks(
    const Model &model,
    const std::vector<std::vector<int>> &potential_links_list,
    std::vector<Curve> &out_discretized_curves,
    std::vector<std::vector<int>> &out_potential_links_unique_per_curve) {
  out_discretized_curves = Discretize(model.GetCurves(), potential_links_list);
  out_potential_links_unique_per_curve =
      GetPotentialLinksUniqueListPerCurve(out_discretized_curves);
}

/// This calls the functions in compute_link.h to compute the linking numbers
/// using Barnes–Hut.
///
/// @param [in] discretized_curves is the collection of discretized curves.
/// @param [in] potential_links_unique_per_curve is the list of unique potential
///     links
/// @param [in] barnes_hut_init_beta is the initial Barnes–Hut beta parameter
///     in the first run, before error estimation.
/// @param [in] barnes_hut_beta_limit is the highest allowed beta to be used in
///     Barnes–Hut.
/// @param [out] out_result is a sparse lower triangular (row < col) matrix that
///     will be initialized to the right size, with its values set to the
///     linking numbers.
void ComputeLinkingNumbersBarnesHut(
    const std::vector<Curve> &discretized_curves,
    const std::vector<std::vector<int>> &potential_links_unique_per_curve,
    double barnes_hut_init_beta, double barnes_hut_beta_limit,
    Eigen::SparseMatrix<int> &out_result) {
  const int ncurves = discretized_curves.size();
  // Precompute segment trees and moments, and count the number of potentially
  // linked loop pairs.
  std::vector<Eigen::ModifiedKdBVH<double, 3, SegmentIntEntry>> segment_trees;
  segment_trees = MakeSegmentTrees(discretized_curves);
  std::vector<std::vector<Moment>> segment_tree_moments(ncurves);
  // This, nplinkssize, is the total number of potentially linked loop pairs.
  size_t nplinkssize = 0;
#pragma omp parallel for reduction(+ : nplinkssize)
  for (int i = 0; i < ncurves; ++i) {
    std::vector<Moment> &tree_moments = segment_tree_moments[i];
    const Eigen::ModifiedKdBVH<double, 3, SegmentIntEntry> &segment_tree =
        segment_trees[i];
    tree_moments = std::vector<Moment>(segment_tree.getBoxArraySize() +
                                       discretized_curves[i].get_q().size());
    int root = segment_tree.getRootIndex();
    // Recursively build tree.
    BuildTreeMoments(segment_tree, discretized_curves[i].GetMidpoints(),
                     discretized_curves[i].GetTangents(), tree_moments, root);
    nplinkssize += potential_links_unique_per_curve[i].size();
  }
  // Preallocate the output matrix's triplet list.
  std::vector<Eigen::Triplet<int>> out_triplet_list(nplinkssize);
  size_t current_output_index = 0;

  // Note: We print the off-integer errors as a guide for the user. This is NOT
  // the error metric we used in our paper [Qu et al. 2021]. For our paper, we
  // just compared each computed value against a ground-truth sparse matrix.
  double off_integer_error = 0;
  double max_off_integer_error = 0;
  // When the number of curves is small, we parallelize the tree.
  bool use_parallel_tree_eval = ncurves < NCURVES_PARALLEL_THRESHOLD;
  // Otherwise, we parallelize iterating the loop pairs
#ifdef _MSC_VER
#pragma omp parallel for reduction(+:off_integer_error) \
    if (!use_parallel_tree_eval)
  for (int index = 0; index < ncurves; ++index) {
#else
#pragma omp parallel for reduction(+:off_integer_error) \
    reduction(max:max_off_integer_error) if (!use_parallel_tree_eval)
  for (int index = 0; index < ncurves; ++index) {
#endif
    int i = index;
    // Load balancing index swap, because we use static scheduling with OpenMP.
    if (i % 2 == 0) {
      i = ncurves - 2 - i + (ncurves % 2);
    }
    // Iterate the potentially linked curves.
    for (int j_ind = 0; j_ind < potential_links_unique_per_curve[i].size();
         ++j_ind) {
      int j = potential_links_unique_per_curve[i][j_ind];
      double error_est = 0;
      double linking_number = EvalLinkingNumberBarnesHut(
          segment_trees[i], segment_trees[j], segment_tree_moments[i],
          segment_tree_moments[j], CERTIFICATE_EPSSQ,
          barnes_hut_init_beta * barnes_hut_init_beta, error_est,
          use_parallel_tree_eval);

      error_est *= ERROR_EST_SCALE;
      double init_beta = barnes_hut_init_beta;
      double target_beta =
          std::min(std::pow(error_est / TARGET_ERROR, 0.25) * init_beta,
                   barnes_hut_beta_limit);
      if (target_beta > init_beta) {
        if (use_parallel_tree_eval) {
          // Print the error estimate only if it's single-threaded, to not spam
          // the user.
          std::cout << "Error estimate of " << error_est
                    << " exceeded the target of " << TARGET_ERROR
                    << "; retrying with beta = " << target_beta << "."
                    << std::endl;
        }
        error_est = 0;
        linking_number = EvalLinkingNumberBarnesHut(
            segment_trees[i], segment_trees[j], segment_tree_moments[i],
            segment_tree_moments[j], CERTIFICATE_EPSSQ,
            target_beta * target_beta, error_est, use_parallel_tree_eval);
      }
      double nearest_int = std::round(linking_number);
      double local_off_integer_error = std::abs(linking_number - nearest_int);
      off_integer_error += local_off_integer_error;
#ifndef _MSC_VER
      max_off_integer_error =
          std::max(max_off_integer_error, local_off_integer_error);
#endif
      // Push the value into the table if it's nonzero.
      // If it is NaN, we output an error message but proceed quietly.
      if (!(std::abs(nearest_int) < 0.25)) {
        size_t out_ind;
#ifdef _MSC_VER
#pragma omp critical
        {
          out_ind = current_output_index++;
        }
#else
#pragma omp atomic capture
        out_ind = current_output_index++;
#endif
        out_triplet_list[out_ind] = Eigen::Triplet<int>(j, i, nearest_int);
        if (std::isnan(nearest_int)) {
#pragma omp critical
          {
            std::cerr << "Linking number (" << i << "," << j << ") is NaN."
                      << std::endl;
          }
        }
      }
    }
  }
  if (nplinkssize > 0) {
    off_integer_error /= nplinkssize;
    std::cout << "Linking-Number Computation, Mean Off-Integer Error: "
              << off_integer_error << std::endl;

#ifndef _MSC_VER
    std::cout << "Linking-Number Computation, Max  Off-Integer Error: "
              << max_off_integer_error << std::endl;
#endif
  }
  out_triplet_list.resize(current_output_index);
  out_result = Eigen::SparseMatrix<int>(ncurves, ncurves);
  out_result.setFromTriplets(out_triplet_list.begin(), out_triplet_list.end());
}

/// This calls the functions in compute_link.h to compute the linking numbers
/// using Direct summation.
///
/// @param [in] discretized_curves is the collection of discretized curves.
/// @param [in] potential_links_unique_per_curve is the list of unique potential
///     links
///
/// @param [out] out_result is a sparse lower triangular (row < col) matrix that
///     will be initialized to the right size, with its values set to the
///     linking numbers.
void ComputeLinkingNumbersDirectSum(
    std::vector<Curve> &discretized_curves,
    std::vector<std::vector<int>> &potential_links_unique_per_curve,
    Eigen::SparseMatrix<int> &out_result) {
  const int ncurves = discretized_curves.size();

  // Count the number of potentially linked loop pairs.
  size_t nplinkssize = 0;
#pragma omp parallel for reduction(+ : nplinkssize) if (ncurves > 10000)
  for (int i = 0; i < ncurves; ++i) {
    nplinkssize += potential_links_unique_per_curve[i].size();
  }
  // Preallocate the output matrix's triplet list.
  std::vector<Eigen::Triplet<int>> out_triplet_list(nplinkssize);
  size_t current_output_index = 0;

  // Note: We print the off-integer errors as a guide for the user. This is NOT
  // the error metric we used in our paper [Qu et al. 2021]. For our paper, we
  // just compared each computed value against a ground-truth sparse matrix.
  double off_integer_error = 0;
  double max_off_integer_error = 0;
  // When the number of curves is small, we parallelize the direct sum.
  bool use_parallel_direct_sum = ncurves < NCURVES_PARALLEL_THRESHOLD;
  // Otherwise, we parallelize iterating the loop pairs.
#ifdef _MSC_VER
#pragma omp parallel for reduction(+:off_integer_error) if (!use_parallel_direct_sum)
  for (int index = 0; index < ncurves; ++index) {
#else
#pragma omp parallel for reduction(+:off_integer_error) \
    reduction(max:max_off_integer_error) if (!use_parallel_direct_sum)
  for (int index = 0; index < ncurves; ++index) {
#endif
    int i = index;
    // Load balancing index swap, because we use static scheduling with OpenMP.
    if (i % 2 == 0) {
      i = ncurves - 2 - i + (ncurves % 2);
    }
    // Iterate the potentially linked curves.
    for (int j_ind = 0; j_ind < potential_links_unique_per_curve[i].size();
         ++j_ind) {
      int j = potential_links_unique_per_curve[i][j_ind];
      double linking_number = std::numeric_limits<double>::quiet_NaN();
      if (!use_parallel_direct_sum) {
        // Single-threaded evaluation, useful for, e.g. chainmail and gloves,
        // because we are parallelizing on the loop pair level.
        linking_number = ComputeLinkingNumberDirectFast(discretized_curves[i],
                                                        discretized_curves[j]);
      } else {
        // Multithreaded evaluation, useful for, e.g. ribbons and torii, because
        // there're only a few, or one, loop pairs.
        linking_number = ComputeLinkingNumberDirectOneArctanPerPair(
            discretized_curves[i], discretized_curves[j], true /* parallel */);
      }
      double nearest_int = std::round(linking_number);
      double local_off_integer_error = std::abs(linking_number - nearest_int);
      off_integer_error += local_off_integer_error;
#ifndef _MSC_VER
      max_off_integer_error =
          std::max(max_off_integer_error, local_off_integer_error);
#endif
      // Push the value into the table if it's nonzero.
      // If it is NaN, we output an error message but proceed quietly.
      if (!(std::abs(nearest_int) < 0.25)) {
        size_t out_ind;
#ifdef _MSC_VER
#pragma omp critical
        {
          out_ind = current_output_index++;
        }
#else
#pragma omp atomic capture
        out_ind = current_output_index++;
#endif
        out_triplet_list[out_ind] = Eigen::Triplet<int>(j, i, nearest_int);
        if (std::isnan(nearest_int)) {
#pragma omp critical
          {
            std::cerr << "Linking number (" << i << "," << j << ") is NaN."
                      << std::endl;
          }
        }
      }
    }
  }

  if (nplinkssize > 0) {
    off_integer_error /= nplinkssize;
    std::cout << "Linking-Number Computation, Mean Off-Integer Error: "
              << off_integer_error << std::endl;
#ifndef _MSC_VER
    std::cout << "Linking-Number Computation, Max  Off-Integer Error: "
              << max_off_integer_error << std::endl;
#endif
  }
  out_triplet_list.resize(current_output_index);
  out_result = Eigen::SparseMatrix<int>(ncurves, ncurves);
  out_result.setFromTriplets(out_triplet_list.begin(), out_triplet_list.end());
}

/**
 * This method computes the linking-number certificate from an input model. By
 * default, it will use Barnes–Hut (Section 3.4.4, Appendix B.2, and the
 * Supplemental docment, of [Qu and James 2021]) to compute the linking
 * number. The force_direct_sum option forces it to use Direct Summation
 * (Section 3.4.2 of [Qu and James 2021]) instead of Barnes–Hut.
 *
 * @param [in] model is the input collection of curves.
 * @param [in] force_direct_sum determines whether to use Direct Summation
 *     instead of Barnes–Hut.
 * @param [in] barnes_hut_init_beta is the initial Barnes–Hut beta parameter
 *     in the first run, before error estimation.
 * @param [in] barnes_hut_beta_limit is the highest allowed beta to be used in
 *     Barnes–Hut.
 *
 * @post the certificate_matrix_ matrix is created with the linking numbers
 *     between curves.
 */
void LinkingNumberCertificate::ComputeFromModel(const Model &model,
                                                bool force_direct_sum,
                                                double barnes_hut_init_beta,
                                                double barnes_hut_beta_limit) {
  std::cout << "Starting Potential Link Search." << std::endl;
  std::vector<std::vector<int>> potential_links = PotentialLinkSearch(model);
  std::vector<Curve> discretized_curves;
  std::vector<std::vector<int>> potential_links_unique_per_curve;

  if(model.get_segment_type() != SegmentType::Polyline){
  std::cout << "Starting Discretization." << std::endl;
  DiscretizeAndGetNewPotentialLinks(model, potential_links, discretized_curves,
                                    potential_links_unique_per_curve);
                                    }
  else{
    std::cout << "Skipping Discretization for Polylines." << std::endl;
    discretized_curves = model.GetCurves();
    potential_links_unique_per_curve = potential_links;
  }
  if (force_direct_sum) {
    std::cout << "Starting Direct Summation." << std::endl;
    ComputeLinkingNumbersDirectSum(discretized_curves,
                                   potential_links_unique_per_curve,
                                   certificate_matrix_);
  } else {
    std::cout << "Starting Barnes–Hut." << std::endl;
    ComputeLinkingNumbersBarnesHut(
        discretized_curves, potential_links_unique_per_curve,
        barnes_hut_init_beta, barnes_hut_beta_limit, certificate_matrix_);
  }
  // Prune any zero entries.
  certificate_matrix_.prune(1);
}

/// This helper subroutine reads an Eigen::SparseMatrix<int> from a text file of
/// triplets. The first line must be the number of curves, and every row after
/// is treated as a (col, row, value) triplet.
///
/// @param [in] filename is the name of the input file.
/// @param [out] out_matrix is the output sparse matrix.
void ReadMatrix(const std::string &filename,
                Eigen::SparseMatrix<int> &out_matrix) {
  typedef Eigen::Triplet<int> Trip;
  std::vector<Trip> triplet_list;
  std::ifstream file(filename);
  std::string line;
  int nrows;
  if (std::getline(file, line)) {
    nrows = std::atoi(line.c_str());
    while (std::getline(file, line)) {
      std::stringstream linestream(line);
      std::string value;
      int row, col, val;
      if (std::getline(linestream, value, ',')) {
        col = std::atoi(value.c_str());
        if (std::getline(linestream, value, ',')) {
          row = std::atoi(value.c_str());
          if (std::getline(linestream, value, ',')) {
            val = std::atoi(value.c_str());
            triplet_list.push_back(Trip(row, col, val));
          }
        }
      }
    }
  }
  file.close();
  out_matrix.resize(nrows, nrows);
  out_matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());
}

/// This helper subroutine saves an Eigen::SparseMatrix<int> into a text file of
/// triplets. The first line is the number of curves, and every row after
/// is treated as a (col, row, value) triplet. If omit_matrix_size is set, this
/// subroutine will omit saving the first row and only save the triplets.
///
/// @param [in] filename is the name of the output file.
/// @param [in] certificate_matrix is the sparse matrix to save.
/// @param [in] omit_matrix_size indicates whether to skip saving the first
///     line, which is the number of rows.
void SaveMatrix(const std::string &filename,
                const Eigen::SparseMatrix<int> &certificate_matrix,
                bool omit_matrix_size) {
  std::ofstream file(filename);
  if (!omit_matrix_size)
    file << certificate_matrix.outerSize() << std::endl;
  for (int k = 0; k < certificate_matrix.outerSize(); ++k) {
    for (Eigen::SparseMatrix<int>::InnerIterator it(certificate_matrix, k); it;
         ++it) {
      file << it.col() << "," << it.row() << "," << it.value() << std::endl;
    }
  }
  file.close();
}

/// This helper function takes a sparse matrix and computes a list of indices
/// that either point to a row with nonzero entries or a column with nonzero
/// entries, indicating a curve participated in nonzero linkages (or violations,
/// in the case of a difference matrix).
///
/// @param [in] certificate_matrix is the sparse matrix.
/// @param [out] out_active_rows_and_cols is the list of indices with nonzero
///     entries in their row or column.
void GetActiveRowsAndCols(const Eigen::SparseMatrix<int> &certificate_matrix,
                          std::unordered_set<int> &out_active_rows_and_cols) {
  for (int k = 0; k < certificate_matrix.outerSize(); ++k) {
    for (Eigen::SparseMatrix<int>::InnerIterator it(certificate_matrix, k); it;
         ++it) {
      out_active_rows_and_cols.insert(it.row());
      out_active_rows_and_cols.insert(it.col());
    }
  }
}

/**
 * This method, ExportToFile, exports the linking-number certificate to a
 * certificate text file. The first line consists of the square matrix's size,
 * and every subsequent line is a triplet of (col index, row index,
 * linking-number value). We swapped the row and column order in the triplets
 * so that the outputs are naturally in ascending order.
 *
 * @param [in] filename is the name of the output certificate text file.
 */
void LinkingNumberCertificate::ExportToFile(const std::string &filename) const {
  SaveMatrix(filename, certificate_matrix_, false /* omit_matrix_size */);
}

/**
 * This method, ImportFromFile, imports an already-computed linking-number
 * certificate from an input certificate text file. The first line consists of
 * the square matrix's size, and every subsequent line is a triplet of (col
 * index, row index, linking-number value). We swapped the row and column
 * order in the triplets so that the outputs are naturally in ascending order.
 *
 * @param [in] filename is the name of the input certificate text file.
 * @post the certificate_matrix_ matrix is created with the contents of the
 *     input file.
 */
void LinkingNumberCertificate::ImportFromFile(const std::string &filename) {
  ReadMatrix(filename, certificate_matrix_);
}

/**
 * This method, ExportDiffsToFiles, compares against another
 * LinkingNumberCertificate and exports lists of:
 *
 * (1) loop pairs that are different into out_diff_links_filename, and
 * (2) curves that participate in differing loop pairs into
 *     out_diff_curves_filename.
 *
 * @param [in] out_diff_links_filename is the filename for exporting the
 *     different loop pairs. These loop pairs are exported as triplets of (col
 *     index, row index, linkage difference), with col < row. Unlike
 *     ExportToFile, we do not output the matrix size in the first row.
 * @param [in] out_diff_curves_filename is the filename for exporting the
 *     curves that participate in differing loop pairs. One curve index is
 *     stored per row.
 * @param [in] other is the other LinkingNumberCertificate to compare against.
 */
void LinkingNumberCertificate::ExportDiffsToFiles(
    const std::string &out_diff_links_filename,
    const std::string &out_diff_curves_filename,
    const LinkingNumberCertificate &other) const {
  if (certificate_matrix_.rows() != other.certificate_matrix_.rows()) {
    std::cerr << "Unable to compare: The two models have different numbers of "
                 "curves: "
              << certificate_matrix_.rows() << " vs. "
              << other.certificate_matrix_.rows() << "." << std::endl;
    std::ofstream diff_links_file(out_diff_links_filename);
    diff_links_file << "Unable to compare: The two models have different "
                       "numbers of curves: "
                    << certificate_matrix_.rows() << " vs. "
                    << other.certificate_matrix_.rows() << "." << std::endl;
    diff_links_file.close();
    std::ofstream diff_curves_file(out_diff_curves_filename);
    diff_curves_file << "Unable to compare: The two models have different "
                        "numbers of curves: "
                     << certificate_matrix_.rows() << " vs. "
                     << other.certificate_matrix_.rows() << "." << std::endl;
    diff_curves_file.close();
    return;
  }
  Eigen::SparseMatrix<int> difference =
      certificate_matrix_ - other.certificate_matrix_;
  difference.prune(1);
  if (out_diff_links_filename.compare("[none]") != 0) {
    SaveMatrix(out_diff_links_filename, difference,
               true /* omit_matrix_size */);
  }
  if (out_diff_curves_filename.compare("[none]") != 0) {
    std::unordered_set<int> active_rows_and_cols;
    GetActiveRowsAndCols(certificate_matrix_, active_rows_and_cols);
    std::ofstream diff_curves_file(out_diff_curves_filename);
    for (int i = 0; i < difference.outerSize(); ++i) {
      if (active_rows_and_cols.find(i) != active_rows_and_cols.end()) {
        diff_curves_file << i << std::endl;
      }
    }
    diff_curves_file.close();
  }
}

} // namespace VerifyCurves