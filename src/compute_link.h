/**
 * Copyright (C) 2021 Ante Qu <antequ@cs.stanford.edu>.
 */

#pragma once

#include "BVH.h"
#include "curve.h"

#include <unordered_map>
#include <utility>
#include <vector>

namespace VerifyCurves {

/**
 * This module computes pairwise linking numbers from collections of discretized
 * curves, based on Section 3.4 of [Qu and James 2021]. The first two functions,
 * ComputeLinkingNumberDirectOneArctanPerPair and
 * ComputeLinkingNumberDirectFast, use direct summation (Section 3.4.2 of [Qu
 * and James 2021]). The rest of this module prepares and evaluates the 2-tree
 * Barnes–Hut method (Section 3.4.4, Appendix B.2, and the Supplemental docment,
 * of [Qu and James 2021]).
 *
 * References:
 *
 * - [Qu and James 2021]         Ante Qu and Doug L. James. 2021. Fast Linking
 *                               Numbers for Topology Verification of Loopy
 *                               Structures. ACM Trans. Graph. 40, 4, Article
 *                               106 (August 2021), 19 pages.
 * - [Arai 2013]                 Zin Arai. 2013. A rigorous numerical algorithm
 *                               for computing the linking number of links.
 *                               Nonlinear Theory and its Applications.
 * - [Bertolazzi et al. 2019]    Enrico Bertolazzi, Riccardo Ghiloni, and Ruben
 *                               Specogna. 2019. Efficient computation of
 *                               linking number with certification. arXiv
 *                               preprint.
 * - [Burscher and Pingali 2011] Martin Burscher and Keshav Pingali. 2011. An
 *                               efficient CUDA implementation of the tree-based
 *                               Barnes Hut N-body algorithm. In GPU Computing
 *                               Gems Emerald Edition.
 */

/**
 * This function, ComputeLinkingNumberDirectOneArctanPerPair, computes the
 * linking number by directly summing the Gauss linking integral, using the
 * exact finite-segment-pair expression from [Arai 2013]. Within each segment
 * pair, it combines the two arctan calls into one arctan. This performs
 * slightly better than ComputeLinkingNumberDirectFast when parallelization is
 * enabled. We therefore call this on models with few curves, such as the
 * double-helix ribbons and the torii, where there may be only one
 * linking-number evaluation. This means we benefit the most from parallelizing
 * within the linking-number evaluation rather than at the loop pair level.
 *
 * @param [in] curve1, the first curve, which should be a discretized polyline.
 * @param [in] curve2, the second curve, which should be another discretized
 *     polyline
 * @param [in] parallel, determines whether to parallelize this evaluation.
 *
 * @returns a double, representing the linking number between the two curves.
 */
double ComputeLinkingNumberDirectOneArctanPerPair(const Curve &curve1,
                                                  const Curve &curve2,
                                                  bool parallel = false);

/**
 * This function, ComputeLinkingNumberDirectFast, computes the linking number by
 * directly summing the Gauss linking integral, using the exact
 * finite-segment-pair expression from [Arai 2013], combined with the arctan
 * addition optimizations from [Bertolazzi et al. 2019]. In particular, for each
 * segment on curve 1, we repeatedly use the arctan addition formula while
 * iterating segments on curve 2, so that only one arctan is called for each
 * segment on curve 1. This performs slightly better than
 * ComputeLinkingNumberDirectOneArctanPerPair when single-threaded. We therefore
 * call this on models with many curves, such as the glove and the chainmail,
 * because we can parallelize across loop pairs rather than within each
 * linking-number evaluation.
 *
 * @param [in] curve1, the first curve, which should be a discretized polyline.
 * @param [in] curve2, the second curve, which should be another discretized
 *     polyline
 * @param [in] parallel, determines whether to parallelize this evaluation.
 *
 * @returns a double, representing the linking number between the two curves.
 */
double ComputeLinkingNumberDirectFast(const Curve &curve1, const Curve &curve2,
                                      bool parallel = false);

/// SegmentIntAddr stores the (curve index, segment index) for a segment.
typedef std::pair<int, int> SegmentIntAddr;
/// SegmentIntEntry stores a (SegmentIntAddr, Bounding Box) for a segment.
typedef std::pair<SegmentIntAddr, Box3d> SegmentIntEntry;

/// These are the second-order moments we store for CPU Barnes-Hut.
struct Moment {
  /// The quadrupole moments are separated along the first dimension, that is,
  /// quadrupole0 stores C_Q₀ᵢⱼ, quadrupole1 stores C_Q₁ᵢⱼ, and quadrupole2
  /// stores C_Q₂ᵢⱼ.
  Eigen::Matrix3d quadrupole0 = Eigen::Matrix3d::Zero();
  Eigen::Matrix3d quadrupole1 = Eigen::Matrix3d::Zero();
  Eigen::Matrix3d quadrupole2 = Eigen::Matrix3d::Zero();
  /// This is C_D, the dipole moment.
  Eigen::Matrix3d dipole = Eigen::Matrix3d::Zero();
  /// This is cₘ, the monopole moment, which for a line segment is in the
  /// direction of its tangent, and for a node just the sum of all its child
  /// monopole moments.
  Eigen::Vector3d tangent = Eigen::Vector3d::Zero();
  /// This is r̃, the location of the moment, which is the center of the node or
  /// midpoint of the line segment.
  Eigen::Vector3d midpoint = Eigen::Vector3d::Zero();
  double tannorm = 0;
  double dipolenorm = 0;
  double quadrupolenorm = 0;
};

/**
 * This function, MakeSegmentTrees, builds a tree of segments for each curve in
 * preparation for CPU Barnes-Hut. It is called once to generate segment trees
 * for all curves. These segment trees are needed for MakeTreeMoments, and both
 * are needed as a precomputation before evaluating every potentially-linked
 * loop pair with Barnes-Hut.
 *
 * @param [in] curves, a list of curves already discretized into polylines.
 *
 * @returns a vector of segment trees, one tree for each curve, where the leaf
 *     elements are segment addresses and bounding boxes.
 */
std::vector<Eigen::ModifiedKdBVH<double, 3, SegmentIntEntry>>
MakeSegmentTrees(const std::vector<Curve> &curves);

/**
 * This subroutine, BuildTreeMoments, builds a tree of moments for a single
 * curve in preparation for CPU Barnes-Hut. It is called once for each curve as
 * a precomputation before evaluating Barnes-Hut, and it requires the segment
 * trees from MakeSegmentTrees.
 *
 * Input:
 * @param [in] segment_tree, a KdBVH of segments in the curve.
 * @param [in] curve_midpoints, a list of segment midpoints.
 * @param [in] curve_tangents, a list of segment length vectors.
 * @param [in] node, the index to the root node of the segment tree.
 *
 * Output:
 * @param [out] tree_moments, a vector that lays out the moments of each node in
 *     segment_tree, concatenated with the moments of each leaf in segment_tree.
 */
void BuildTreeMoments(
    const Eigen::ModifiedKdBVH<double, 3, SegmentIntEntry> &segment_tree,
    const std::vector<Point3> &curve_midpoints,
    const std::vector<Point3> &curve_tangents,
    std::vector<Moment> &tree_moments, int node);

/**
 * This function, EvalLinkingNumberBarnesHut, evaluates the 2-tree Barnes-Hut
 * algorithm described in Appendix B.2 of [Qu and James 2021], with the
 * additional terms from the supplemental document. It takes a segment tree and
 * the tree moments built for each curve, and returns a linking number.
 *
 * @param [in] segment_tree1, a KdBVH of segments in curve 1.
 * @param [in] segment_tree2, a KdBVH of segments in curve 2.
 * @param [in] tree_moments1, a vector of moments of segment_tree1 from
 *     BuildTreeMoments.
 * @param [in] tree_moments2, a vector of moments of segment_tree2 from
 *     BuildTreeMoments.
 * @param [in] epssq, the square of machine epsilon used to prevent
 *     singularities when the distance between nodes is zero. This variable name
 *     comes from [Burtscher and Pingali 2011].
 * @param [in] itolsq, the square of beta; this variable name also comes from
 *     [Burtscher and Pingali 2011].
 * @param [in] parallel, determines whether to parallelize this evaluation.
 *
 * @param [out] evals_errs_total, the next-order error estimate is accumulated
 *     into this value.
 * @returns a double, representing the linking number between the two curves.
 */
double EvalLinkingNumberBarnesHut(
    const Eigen::ModifiedKdBVH<double, 3, SegmentIntEntry> &segment_tree1,
    const Eigen::ModifiedKdBVH<double, 3, SegmentIntEntry> &segment_tree2,
    const std::vector<Moment> &tree_moments1,
    const std::vector<Moment> &tree_moments2, const double epssq,
    const double itolsq, double &eval_errs_total, bool parallel = true);

} // namespace VerifyCurves