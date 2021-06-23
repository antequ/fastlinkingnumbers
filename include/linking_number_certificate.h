/**
 * Copyright (C) 2021 Ante Qu <antequ@cs.stanford.edu>.
 */

#pragma once
#include "model.h"

#include <Eigen/Sparse>

namespace VerifyCurves {

/**
 * A LinkingNumberCertificate stores a linking-number certificate as a sparse,
 * lower-triangular (row > column) matrix. This class handles certificate
 * computation from an input curve-collection Model, using [Qu and James 2021].
 * It also handles I/O for saving and loading a certificate.
 *
 * Reference:
 *
 * - [Qu and James 2021]         Ante Qu and Doug L. James. 2021. Fast Linking
 *                               Numbers for Topology Verification of Loopy
 *                               Structures. ACM Trans. Graph. 40, 4, Article
 *                               106 (August 2021), 19 pages.
 */

class LinkingNumberCertificate {
public:
  /**
   * The default constructor simply constructs a blank certificate.
   */
  LinkingNumberCertificate();

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
  LinkingNumberCertificate(const std::string &filename);

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
  void ComputeFromModel(const Model &model, bool force_direct_sum = false,
                        double barnes_hut_init_beta = 2.0,
                        double barnes_hut_beta_limit = 10.0);

  /**
   * This method, ExportToFile, exports the linking-number certificate to a
   * certificate text file. The first line consists of the square matrix's size,
   * and every subsequent line is a triplet of (col index, row index,
   * linking-number value). We swapped the row and column order in the triplets
   * so that the outputs are naturally in ascending order.
   *
   * @param [in] filename is the name of the output certificate text file.
   */
  void ExportToFile(const std::string &filename) const;

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
  void ImportFromFile(const std::string &filename);

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
  void ExportDiffsToFiles(const std::string &out_diff_links_filename,
                          const std::string &out_diff_curves_filename,
                          const LinkingNumberCertificate &other) const;

  /**
   * These methods and operators check for equality between
   * LinkingNumberCertificates.
   */
  bool Equals(const LinkingNumberCertificate &other) const {
    return (
        (certificate_matrix_.rows() == other.certificate_matrix_.rows()) &&
        ((certificate_matrix_ - other.certificate_matrix_).squaredNorm() == 0));
  };

  bool operator==(const LinkingNumberCertificate &other) const {
    return Equals(other);
  }
  bool operator!=(const LinkingNumberCertificate &other) const {
    return !(Equals(other));
  }

private:
  // Certificate Data. We maintain that the row > column for every entry.
  Eigen::SparseMatrix<int> certificate_matrix_;
};

} // namespace VerifyCurves