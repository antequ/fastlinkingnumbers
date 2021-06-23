/**
 * Copyright (C) 2021 Ante Qu <antequ@cs.stanford.edu>.
 */

#include "compute_link.h"

#include <iostream>
#include <vector>

// Predefine 1/sqrt(3), 1/(2 π), 1/(4 π), and 4 π
#define INV_RT3                                                                \
  (0.5773502691896257645091487805019574556476017512701268760186023264)
#define INV_TWO_PI                                                             \
  (0.1591549430918953357688837633725143620344596457404564487476673440)
#define INV_FOUR_PI                                                            \
  (0.0795774715459476678844418816862571810172298228702282243738336720)
#define FOUR_PI                                                                \
  (12.566370614359172953850573533118011536788677597500423283899778369)

namespace Eigen {
template <typename T>
inline VerifyCurves::Box3d
bounding_box(const std::pair<T, VerifyCurves::Box3d> &spline) {
  return spline.second;
}
} // namespace Eigen

namespace VerifyCurves {

// These are direct summation helper routines. While they are as fast as I could
// get them, they are still much slower than what should be possible with
// vectorization, and initial profiling showed some undesired bottlenecks, too.

// These two compute the triple scalar product a ⋅ (b × c). determ acts on row
// vectors; determ_col acts on column vectors.
inline Eigen::RowVectorXd determ(const Eigen::Matrix3Xd &a,
                                 const Eigen::Matrix3Xd &b,
                                 const Eigen::Matrix3Xd &c) {
  return a.row(0).cwiseProduct(b.row(1).cwiseProduct(c.row(2)) -
                               b.row(2).cwiseProduct(c.row(1))) +
         a.row(2).cwiseProduct(b.row(0).cwiseProduct(c.row(1)) -
                               b.row(1).cwiseProduct(c.row(0))) +
         a.row(1).cwiseProduct(b.row(2).cwiseProduct(c.row(0)) -
                               b.row(0).cwiseProduct(c.row(2)));
}

inline Eigen::VectorXd determ_col(const Eigen::MatrixX3d &a,
                                  const Eigen::MatrixX3d &b,
                                  const Eigen::MatrixX3d &c) {
  return a.col(0).cwiseProduct(b.col(1).cwiseProduct(c.col(2)) -
                               b.col(2).cwiseProduct(c.col(1))) +
         a.col(2).cwiseProduct(b.col(0).cwiseProduct(c.col(1)) -
                               b.col(1).cwiseProduct(c.col(0))) +
         a.col(1).cwiseProduct(b.col(2).cwiseProduct(c.col(0)) -
                               b.col(0).cwiseProduct(c.col(2)));
}

// These two compute the solid angles needed for the denominator
// of the [Arai 2013] expression. a,b,c,d are inputs, tmp1, tmp2,
// and tmp3 are scratch space, and div1 and div2 are outputs.
// solid_angle_triangle_denoms acts on row vectors;
// solid_angle_triangle_denoms_col acts on column vectors.
inline void solid_angle_triangle_denoms(
    const Eigen::Matrix3Xd &a, const Eigen::Matrix3Xd &b,
    const Eigen::Matrix3Xd &c, const Eigen::Matrix3Xd &d,
    Eigen::RowVectorXd &tmp1, Eigen::RowVectorXd &tmp2,
    Eigen::RowVectorXd &tmp3, Eigen::RowVectorXd &div1,
    Eigen::RowVectorXd &div2) {
  tmp1 = a.colwise().norm();
  tmp2 = b.colwise().norm();
  tmp3 = c.colwise().norm();
  // solid angle triangle(alpha, beta, gamma)
  div2 = c.cwiseProduct(a).colwise().sum();
  div1 = tmp1.cwiseProduct(tmp2).cwiseProduct(tmp3) +
         a.cwiseProduct(b).colwise().sum().cwiseProduct(tmp3) +
         div2.cwiseProduct(tmp2) +
         c.cwiseProduct(b).colwise().sum().cwiseProduct(tmp1);
  // solid angle triangle(gamma, delta, alpha)
  tmp2 = d.colwise().norm();
  div2 = tmp1.cwiseProduct(tmp2).cwiseProduct(tmp3) +
         a.cwiseProduct(d).colwise().sum().cwiseProduct(tmp3) +
         div2.cwiseProduct(tmp2) +
         c.cwiseProduct(d).colwise().sum().cwiseProduct(tmp1);
}

inline void solid_angle_triangle_denoms_col(
    const Eigen::MatrixX3d &a, const Eigen::MatrixX3d &b,
    const Eigen::MatrixX3d &c, const Eigen::MatrixX3d &d, Eigen::VectorXd &tmp1,
    Eigen::VectorXd &tmp2, Eigen::VectorXd &tmp3, Eigen::VectorXd &div1,
    Eigen::VectorXd &div2) {
  tmp1 = a.rowwise().norm();
  tmp2 = b.rowwise().norm();
  tmp3 = c.rowwise().norm();
  // solid angle triangle(alpha, beta, gamma)
  div2 = c.cwiseProduct(a).rowwise().sum();
  div1 = tmp1.cwiseProduct(tmp2).cwiseProduct(tmp3) +
         a.cwiseProduct(b).rowwise().sum().cwiseProduct(tmp3) +
         div2.cwiseProduct(tmp2) +
         c.cwiseProduct(b).rowwise().sum().cwiseProduct(tmp1);
  // solid angle triangle(gamma, delta, alpha)
  tmp2 = d.rowwise().norm();
  div2 = tmp1.cwiseProduct(tmp2).cwiseProduct(tmp3) +
         a.cwiseProduct(d).rowwise().sum().cwiseProduct(tmp3) +
         div2.cwiseProduct(tmp2) +
         c.cwiseProduct(d).rowwise().sum().cwiseProduct(tmp1);
}

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
                                                  bool parallel) {
  const int n1 = curve1.get_q().size();
  const int n2 = curve2.get_q().size();
  std::vector<Point3> startPoints1 = curve1.GetStartPoints();
  std::vector<Point3> loopedStartPoints2 = curve2.GetLoopedStartPoints();
  const Eigen::Map<Eigen::Matrix3Xd> g0s(startPoints1[0].data(), 3, n1);
  const Eigen::Map<Eigen::Matrix3Xd> g1s(loopedStartPoints2[0].data(), 3, n2);
  const Eigen::Map<Eigen::Matrix3Xd> g1p1s(loopedStartPoints2[1].data(), 3, n2);
  double linking_number = 0;
  Eigen::Matrix3Xd alpha = Eigen::Matrix3Xd(3, n2);
  Eigen::Matrix3Xd beta = Eigen::Matrix3Xd(3, n2);
  Eigen::Matrix3Xd gamma = Eigen::Matrix3Xd(3, n2);
  Eigen::Matrix3Xd delta = Eigen::Matrix3Xd(3, n2);
  Eigen::RowVectorXd tmp1 = Eigen::RowVectorXd(n2);
  Eigen::RowVectorXd tmp2 = Eigen::RowVectorXd(n2);
  Eigen::RowVectorXd tmp3 = Eigen::RowVectorXd(n2);
  Eigen::RowVectorXd div1 = Eigen::RowVectorXd(n2);
  Eigen::RowVectorXd div2 = Eigen::RowVectorXd(n2);
  Eigen::RowVectorXd determinant = Eigen::RowVectorXd(n2);
  const Eigen::RowVectorXd ones = Eigen::RowVectorXd::Ones(n2);
#pragma omp parallel for reduction(+:linking_number) if(parallel) \
firstprivate(alpha, beta, gamma, delta, tmp1, tmp2, tmp3, div1, div2, \
determinant)
  for (int i = 0; i < n1; ++i) {
    alpha = g1s.colwise() - g0s.col(i);
    beta = g1s.colwise() - g0s.col((i + 1) % n1);
    gamma = g1p1s.colwise() - g0s.col((i + 1) % n1);
    delta = g1p1s.colwise() - g0s.col(i);
    solid_angle_triangle_denoms(alpha, beta, gamma, delta, tmp1, tmp2, tmp3,
                                div1, div2);
    determinant = determ(alpha, beta, gamma);
    // Combining two arctans into one arctan:
    // x'
    tmp1 = div1.cwiseProduct(div2) - determinant.cwiseAbs2();
    // y'
    tmp2 = determinant.cwiseProduct(div1 + div2);
    // 2 * atan2(y, x)
    tmp3 = tmp2.binaryExpr(
        tmp1, [](double a, double b) { return 2.0 * std::atan2(a, b); });
    // check if signs are equal: if new y is a differnet sign than old y, add
    // 4 π * sign(old y)
    tmp1 = tmp2.cwiseProduct(determinant);
    tmp3 =
        tmp3 + determinant.cwiseSign().binaryExpr(tmp1, [](double a, double b) {
          return (b < 0) ? (FOUR_PI)*a : 0.0;
        });
    linking_number += tmp3.sum();
  }
  linking_number *= INV_FOUR_PI;
  return linking_number;
}

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
                                      bool parallel) {
  const int n1 = curve1.get_q().size();
  const int n2 = curve2.get_q().size();
  std::vector<Point3> startPoints1 = curve1.GetStartPoints();
  std::vector<Point3> loopedStartPoints2 = curve2.GetLoopedStartPoints();
  const Eigen::Map<Eigen::Matrix3Xd> rg0s(startPoints1[0].data(), 3, n1);
  const Eigen::Map<Eigen::Matrix3Xd> rg1s(loopedStartPoints2[0].data(), 3, n2);
  const Eigen::Map<Eigen::Matrix3Xd> rg1p1s(loopedStartPoints2[1].data(), 3,
                                            n2);
  Eigen::MatrixX3d g0s = rg0s.transpose();
  Eigen::MatrixX3d g1s = rg1s.transpose();
  Eigen::MatrixX3d g1p1s = rg1p1s.transpose();
  Eigen::MatrixX3d alpha = Eigen::MatrixX3d(n2, 3);
  Eigen::MatrixX3d beta = Eigen::MatrixX3d(n2, 3);
  Eigen::MatrixX3d gamma = Eigen::MatrixX3d(n2, 3);
  Eigen::MatrixX3d delta = Eigen::MatrixX3d(n2, 3);
  Eigen::VectorXd tmp1 = Eigen::VectorXd(n2);
  Eigen::VectorXd tmp2 = Eigen::VectorXd(n2);
  Eigen::VectorXd tmp3 = Eigen::VectorXd(n2);
  Eigen::VectorXd div1 = Eigen::VectorXd(n2);
  Eigen::VectorXd div2 = Eigen::VectorXd(n2);
  Eigen::VectorXd determinant = Eigen::VectorXd(n2);
  Eigen::VectorXd S = -Eigen::VectorXd::Ones(n2);
  Eigen::VectorXd xs = Eigen::VectorXd::Ones(n2);
  Eigen::VectorXd ys = Eigen::VectorXd::Zero(n2);
  // offsets is the number of negative-x-axis crossings; it is always an
  // integer.
  Eigen::VectorXd offsets = Eigen::VectorXd::Zero(n2);
  double linking_number = 0;

  for (int i = 0; i < n1; ++i) {
    alpha = g1s.rowwise() - g0s.row(i);
    beta = g1s.rowwise() - g0s.row((i + 1) % n1);
    gamma = g1p1s.rowwise() - g0s.row((i + 1) % n1);
    delta = g1p1s.rowwise() - g0s.row(i);
    solid_angle_triangle_denoms_col(alpha, beta, gamma, delta, tmp1, tmp2, tmp3,
                                    div1, div2);
    determinant = determ_col(alpha, beta, gamma);
    // Combining two arctans into one arctan:
    // x' for the two
    tmp1 = div1.cwiseProduct(div2) - determinant.cwiseAbs2();
    // y' for the two
    tmp2 = (determinant.cwiseProduct(div1 + div2));
    // new y" = x1 y2 + x2 y1
    div2 = (xs.cwiseProduct(tmp2) + tmp1.cwiseProduct(ys));
    // new x" = x1 x2 - y1 y2
    div1 = (xs.cwiseProduct(tmp1) - ys.cwiseProduct(tmp2));
    // use max(|x''|, |y''|) to scale; otherwise the values get too large.
    tmp1 = div1.cwiseAbs().cwiseMax(div2.cwiseAbs());
    div1 = div1.cwiseQuotient(tmp1);
    div2 = div2.cwiseQuotient(tmp1);
    tmp1 = tmp2.cwiseProduct(determinant);
    offsets += (determinant.cwiseSign().binaryExpr(
        tmp1, [](double a, double b) { return (b < 0) ? a : 0.0; }));
    tmp3 = div1.binaryExpr(div2, [](double a, double b) {
      return ((b > 0.0) || (b == 0.0 && a < 0.0)) ? 1.0 : -1.0;
    });
    tmp1 = tmp3.cwiseProduct(S);
    tmp2 = tmp2.cwiseProduct(S);
    // cross detection expression
    tmp1 = tmp1.binaryExpr(tmp2, [](double a, double b) {
      return ((a < 0) && (b >= 0)) ? 1.0 : 0.0;
    });
    offsets += S.cwiseProduct(tmp1);
    // set xs, ys, and S
    using std::swap;
    swap(xs, div1);
    swap(ys, div2);
    swap(S, tmp3);
  }

  linking_number +=
      ys.binaryExpr(xs,
                    [](double a, double b) {
                      return (2.0 * INV_FOUR_PI) * std::atan2(a, b);
                    })
          .sum();
  linking_number += offsets.sum();
  return linking_number;
}

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
MakeSegmentTrees(const std::vector<Curve> &curves) {
  const int ncurves = curves.size();
  std::vector<std::vector<SegmentIntEntry>> entry_list(
      ncurves); // (entry_list_size);
  std::vector<Eigen::ModifiedKdBVH<double, 3, SegmentIntEntry>> trees(ncurves);
// This parallelization doesn't speed anything up.
#pragma omp parallel for
  for (int yi = 0; yi < curves.size(); ++yi) {
    entry_list[yi] = std::vector<SegmentIntEntry>(curves[yi].get_q().size());
    for (int j = 0; j < curves[yi].get_q().size(); ++j) {
      const Point3 &current_point = curves[yi].get_q()[j];
      const Point3 &next_point =
          curves[yi].get_q()[(j + 1) % curves[yi].get_q().size()];
      entry_list[yi][j] = SegmentIntEntry(
          SegmentIntAddr(yi, j), curves[yi].get_segment_bounding_box(j));
    }
    trees[yi] = Eigen::ModifiedKdBVH<double, 3, SegmentIntEntry>(
        entry_list[yi].begin(), entry_list[yi].end());
  }
  return trees;
}

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
    std::vector<Moment> &tree_moments, int node) {
  const Eigen::AlignedBox3d &bbox = segment_tree.getVolume(node);
  Eigen::Vector3d node_midpoint = 0.5 * (bbox.max() + bbox.min());
  Moment &moment = tree_moments[node];
  moment.midpoint = node_midpoint;
  std::pair<int, int> indices = segment_tree.getChildrenIndices(node);
  const int moment_index_offset = segment_tree.getBoxArraySize();
  for (size_t i = 0; i < 2; ++i) {
    int index = (i == 0) ? indices.first : indices.second;
    Moment &child_moment = tree_moments[index];
    if (index < moment_index_offset) {
      BuildTreeMoments(segment_tree, curve_midpoints, curve_tangents,
                       tree_moments, index);
    } else {
      // This initializes leaf nodes, as per Appendix B.2 of [Qu and James 2021]
      const int spl_index =
          segment_tree.getObject(index - moment_index_offset).first.second;
      child_moment.midpoint = curve_midpoints[spl_index];
      Eigen::Vector3d c1 = curve_tangents[spl_index];
      child_moment.tangent = c1;
      child_moment.tannorm = c1.norm();
      Eigen::Matrix3d c1outerc1 = (c1 * c1.transpose());
      child_moment.quadrupole0 = (1.0 / 12.0) * c1[0] * c1outerc1;
      child_moment.quadrupole1 = (1.0 / 12.0) * c1[1] * c1outerc1;
      child_moment.quadrupole2 = (1.0 / 12.0) * c1[2] * c1outerc1;
      child_moment.quadrupolenorm = (1.0 / 36.0) * child_moment.tannorm *
                                    child_moment.tannorm * child_moment.tannorm;
    }
    // This code shifts the child node on to the parent node, as per Eqs.
    // (21-23) of [Qu and James 2021] Appendix B.2.
    Eigen::Vector3d cc = child_moment.tangent;
    moment.tangent += cc;
    Eigen::Vector3d rc = child_moment.midpoint - node_midpoint;
    moment.dipole += child_moment.dipole + cc * rc.transpose();
    Eigen::Matrix3d rcouterrc = (rc * rc.transpose());
    Eigen::Matrix3d tmp = rc * child_moment.dipole.row(0);
    moment.quadrupole0 +=
        child_moment.quadrupole0 + tmp.transpose() + tmp + cc[0] * rcouterrc;
    tmp = rc * child_moment.dipole.row(1);
    moment.quadrupole1 +=
        child_moment.quadrupole1 + tmp.transpose() + tmp + cc[1] * rcouterrc;
    tmp = rc * child_moment.dipole.row(2);
    moment.quadrupole2 +=
        child_moment.quadrupole2 + tmp.transpose() + tmp + cc[2] * rcouterrc;
  }
  // These norms aid in error estimation.
  moment.tannorm = moment.tangent.norm();
  moment.dipolenorm = moment.dipole.norm() * INV_RT3;
  moment.quadrupolenorm =
      (1.0 / 3.0) * std::sqrt(moment.quadrupole0.squaredNorm() +
                              moment.quadrupole1.squaredNorm() +
                              moment.quadrupole2.squaredNorm());
}

/// This helper function accumulates dipole and quadrupole terms of the
/// Barnes–Hut far-field expansion.
///
/// @param [in] moment1 contains the moments at this node in curve 1.
/// @param [in] moment2 contains the moments at this node in curve 2.
/// @param [in] dr is the displacement vector, r̃₂ - r̃₁, between the node centers
/// @param [in] rdist2 is the reciprocal of the squared norm of dr, |r̃₂ - r̃₁|⁻².
/// @param [in] b1radius is the circumradius (or semidiagonal length) of the
///     bounding box of node 1.
/// @param [in] b2radius is the circumradius (or semidiagonal length) of the
///     bounding box of node 2.
/// @param [inout] eval_error accumulates the next-order error estimate of this
///     evaluation, scaled by 4 π |r|³.
///
/// @returns the linking-number contribution from this pair of nodes, scaled by
///     4 π |r|³.
///
/// @note Important: the outputs assume they will be multiplied by a factor of
///    1/(4 π |r|³), and so do not include this factor.
inline double
EvalFarFieldHOCorrection(const Moment &moment1, const Moment &moment2,
                         const Eigen::Vector3d &dr, const double rdist2,
                         const double b1radius, const double b2radius,
                         double &eval_error) {
  // Refer to Eq. (29) and (30) of the Supplemental document for these
  // expressions. The supplemental document is here:
  // https://graphics.stanford.edu/papers/fastlinkingnumbers/assets/Qu2021BarnesHutExpandedTermsReference.pdf

  // DIPOLE TERMS (Eq. (29))
  // Confusingly, these names (C1,C2,C3,D1,D2,D3) are different from those in
  // the paper. This is how they map (code variable name on the left, paper name
  // on the right):
  //   C1 = c_M₁, a 3-value vector.
  //   C2 = εᵢⱼₖ C_D₁ⱼₖ, a 3-value vector.
  //   C3 = C_D₁, a 3 × 3 matrix.
  //   D1 = c_M₂, a 3-value vector.
  //   D2 = εᵢⱼₖ C_D₂ⱼₖ, a 3-value vector.
  //   D3 = c_D₂, a 3 × 3 matrix.

  // Furthermore, in the code, dr = r̃₂ - r̃₁.

  double result = 0;
  const Eigen::Matrix3d &D3 = moment2.dipole;
  const Eigen::Matrix3d &C3 = moment1.dipole;
  const Eigen::Vector3d &D1 = moment2.tangent;
  const Eigen::Vector3d &C1 = moment1.tangent;
  // We compute D2 and C2. The prime (') is a derivative with respect to the
  // parameter (so it's a tangent vector).
  //   D2 = ∫ r₂'× r₂,
  //   D3 = ∫ r₂' r₂ᵀ, and so
  //   D2 = [D3₂₃ - D3₃₂, D3₃₁ - D3₁₃, D3₁₂ - D3₂₁]
  // The same computation is done for C2 and r₁.
  const Eigen::Vector3d D2 = Eigen::Vector3d(
      D3(1, 2) - D3(2, 1), D3(2, 0) - D3(0, 2), D3(0, 1) - D3(1, 0));
  const Eigen::Vector3d C2 = Eigen::Vector3d(
      C3(1, 2) - C3(2, 1), C3(2, 0) - C3(0, 2), C3(0, 1) - C3(1, 0));
  // This is the first term on the RHS of Eq. (29)
  result -= C1.dot(D2) + C2.dot(D1);
  // And this is the second term on the RHS of Eq. (29).
  result -= 3.0 * rdist2 *
            ((D3 * dr).cross(C1).dot(dr) + (C3 * dr).cross(D1).dot(dr));
  // QUADRUPOLE TERMS (Eq. (30))
  {
    // We first compute T_Q_second, which are the terms in the second row of Eq.
    // (30), negated.

    // Let's compute T_Q_second_first = -εᵢⱼₖ rᵢ(c_M₁ⱼ C_Q₂ₖₗₗ - c_M₂ⱼ C_Q₁ₖₗₗ)
    // from Eq. (30). As for the intermediate variable names, p refers to r₁ and
    // q refers to r₂. So pprimepp is C_Q₁ₖₗₗ, and qprimeqq is C_Q₂ₖₗₗ.
    Eigen::Vector3d pprimepp = {moment1.quadrupole0.trace(),
                                moment1.quadrupole1.trace(),
                                moment1.quadrupole2.trace()};
    Eigen::Vector3d qprimeqq = {moment2.quadrupole0.trace(),
                                moment2.quadrupole1.trace(),
                                moment2.quadrupole2.trace()};
    double T_Q_second_first = -dr.dot(pprimepp.cross(D1) + C1.cross(qprimeqq));
    // Now let's compute T_Q_second_second = 2 εᵢⱼₖ rᵢ C_D₁ⱼₗ C_D₂ₖₗ from Eq.
    // (30).
    double product01 = C3.row(0).dot(D3.row(1));
    double product02 = C3.row(0).dot(D3.row(2));
    double product10 = C3.row(1).dot(D3.row(0));
    double product12 = C3.row(1).dot(D3.row(2));
    double product20 = C3.row(2).dot(D3.row(0));
    double product21 = C3.row(2).dot(D3.row(1));
    double T_Q_second_second = 2.0 * (dr(0) * (product12 - product21) +
                                      dr(1) * (product20 - product02) +
                                      dr(2) * (product01 - product10));
    double T_Q_second = T_Q_second_first + T_Q_second_second;

    // Now let's compute T_Q_third, which is 0.5 times the third row of Eq.
    // (30), negated.
    // First let's compute εᵢⱼₖ c_M₂ᵢ C_Q₁ⱼₖₗ rₗ.
    Eigen::Vector3d pprimeXppTr = {
        (moment1.quadrupole1.row(2) - moment1.quadrupole2.row(1)) * dr,
        (moment1.quadrupole2.row(0) - moment1.quadrupole0.row(2)) * dr,
        (moment1.quadrupole0.row(1) - moment1.quadrupole1.row(0)) * dr};
    double T_Q_third = D1.dot(pprimeXppTr);
    // And now -εᵢⱼₖ c_M₁ᵢ C_Q₂ⱼₖₗ rₗ.
    Eigen::Vector3d qprimeXqqTr = {
        (moment2.quadrupole1.row(2) - moment2.quadrupole2.row(1)) * dr,
        (moment2.quadrupole2.row(0) - moment2.quadrupole0.row(2)) * dr,
        (moment2.quadrupole0.row(1) - moment2.quadrupole1.row(0)) * dr};
    T_Q_third += -C1.dot(qprimeXqqTr);
    // Now let's add εᵢⱼₖ rₗ C_D₁ᵢₗ C_D₂ⱼₖ - εᵢⱼₖ rₗ C_D₂ᵢₗ C_D₁ⱼₖ.
    T_Q_third += dr.dot(C3.transpose() * D2 - D3.transpose() * C2);

    // Now let's compute T_Q_first, the first line of Eq. (30).
    // First let's compute -εᵢⱼₖ rᵢ rₗ rₘ c_M₂ⱼ C_Q₁ₖₗₘ.
    Eigen::Vector3d tensor_product_pprimepp_rr = {
        dr.dot(moment1.quadrupole0 * dr), dr.dot(moment1.quadrupole1 * dr),
        dr.dot(moment1.quadrupole2 * dr)};
    double T_Q_firstpp = -dr.dot(D1.cross(tensor_product_pprimepp_rr));
    // And here's εᵢⱼₖ rᵢ rₗ rₘ c_M₁ⱼ C_Q₂ₖₗₘ.
    Eigen::Vector3d tensor_product_qprimeqq_rr = {
        dr.dot(moment2.quadrupole0 * dr), dr.dot(moment2.quadrupole1 * dr),
        dr.dot(moment2.quadrupole2 * dr)};
    double T_Q_firstqq = dr.dot(C1.cross(tensor_product_qprimeqq_rr));
    // Finally, we should add -2 εᵢⱼₖ rᵢ rₗ rₘ C_D₁ⱼₗ C_D₂ₖₘ
    double T_Q_firstpq = 2.0 * (C3 * dr).dot((dr.cross(D3 * dr)));
    double T_Q_first = T_Q_firstpp + T_Q_firstqq + T_Q_firstpq;

    // To combine the terms, we can subtract all of this times |r|⁻². This gets
    // us t_Q.
    result -= 1.5 * rdist2 *
              ((T_Q_second + 2.0 * T_Q_third) + 5.0 * rdist2 * T_Q_first);
  }
  // This is the error estimate from Eq. (14) of [Qu and James 2021].
  eval_error =
      rdist2 * (moment1.quadrupolenorm *
                    (b1radius * moment2.tannorm + 3 * moment2.dipolenorm) +
                moment2.quadrupolenorm *
                    (b2radius * moment1.tannorm + 3 * moment1.dipolenorm));
  return result;
}

/// This is a recursive helper to evaluate the linking numbers using Barnes–Hut,
/// by traversing the segment tree depth-first.
///
/// @param [in] segment_tree1 is the segment tree of curve 1.
/// @param [in] segment_tree2 is the segment tree of curve 2.
/// @param [in] tree_moments1 contains the moments for segment_tree1.
/// @param [in] tree_moments2 contains the moments for segment_tree2.
/// @param [in] node1 is the root node of the subtree in segment_tree1 we are
///     traversing from.
/// @param [in] node2 is the root node of the subtree in segment_tree2 we are
///     traversing from.
/// @param [in] epssq is the square of machine epsilon used to prevent
///     singularities when the distance between nodes is zero. This variable
///     name comes from [Burtscher and Pingali 2011].
/// @param [in] itolsq is the square of beta; this variable name also comes from
///     [Burtscher and Pingali 2011].
/// @param [in] parallel determines whether to parallelize this evaluation.
/// @param [inout] eval_error accumulates the next-order error estimate for this
///     subtree.
///
/// @returns the linking number accumulated from this pair of subtrees.
double EvalLinkingNumberRecurse(
    const Eigen::ModifiedKdBVH<double, 3, SegmentIntEntry> &segment_tree1,
    const Eigen::ModifiedKdBVH<double, 3, SegmentIntEntry> &segment_tree2,
    const std::vector<Moment> &tree_moments1,
    const std::vector<Moment> &tree_moments2, int node1, int node2,
    double &eval_error, double epssq, double itolsq) {
  const int moment_index_offset1 = segment_tree1.getBoxArraySize();
  const int moment_index_offset2 = segment_tree2.getBoxArraySize();
  const int obj1 = node1 - moment_index_offset1;
  const int obj2 = node2 - moment_index_offset2;
  const Eigen::AlignedBox3d &bbox1 = (obj1 < 0)
                                         ? segment_tree1.getVolume(node1)
                                         : segment_tree1.getObject(obj1).second;
  const Eigen::AlignedBox3d &bbox2 = (obj2 < 0)
                                         ? segment_tree2.getVolume(node2)
                                         : segment_tree2.getObject(obj2).second;
  double radius1sq = 0.25 * (bbox1.max() - bbox1.min()).squaredNorm();
  double radius2sq = 0.25 * (bbox2.max() - bbox2.min()).squaredNorm();
  double radius1 = std::sqrt(radius1sq);
  double radius2 = std::sqrt(radius2sq);
  using std::sqrt;
  double r1plusr2sq = (radius1 + radius2) * (radius1 + radius2);
  const Moment &moment1 = tree_moments1[node1];
  const Moment &moment2 = tree_moments2[node2];
  // Notationwise, we say dr = r = r̃₂ - r̃₁.
  Eigen::Vector3d dr = moment2.midpoint - moment1.midpoint;
  double dist_sq = dr.squaredNorm() + epssq;
  double dist_threshold = r1plusr2sq * itolsq + epssq;

  if (obj1 >= 0 && obj2 >= 0) {
    // Leaf nodes: Use the one-arctan expression for two finite segments.
    double determinant = dr.dot(moment2.tangent.cross(moment1.tangent));
    Eigen::Vector3d alpha = dr + 0.5 * (-moment2.tangent + moment1.tangent);
    Eigen::Vector3d beta = dr + 0.5 * (-moment2.tangent - moment1.tangent);
    Eigen::Vector3d gamma = dr + 0.5 * (moment2.tangent - moment1.tangent);
    double la = alpha.norm();
    double lb = beta.norm();
    double lc = gamma.norm();
    double ac = alpha.dot(gamma);
    double div1 =
        la * lb * lc + alpha.dot(beta) * lc + ac * lb + beta.dot(gamma) * la;
    Eigen::Vector3d delta = dr + 0.5 * (moment2.tangent + moment1.tangent);
    double ld = delta.norm();
    double div2 =
        la * ld * lc + alpha.dot(delta) * lc + ac * ld + delta.dot(gamma) * la;
    // x'
    la = div1 * div2 - determinant * determinant;
    using std::atan2;
    using std::copysign;
    // ac is s'' = y' * determinant. We also want it negative when det = 0 and
    // div1, div2 are both negative.
    ac = (determinant == 0) ? copysign(1.0, div1) + copysign(1.0, div2)
                            : (div1 + div2); // = y' * determinant
    // y'
    lb = determinant * ac;
    // sign of the determinant
    lc = (determinant == 0) ? -copysign(1.0, div1) : copysign(1.0, determinant);
    return INV_TWO_PI * atan2(lb, la) +
           (((ac < 0) || ((ac == 0) && (determinant < 0))) ? lc : 0.0);
  }

  if (dist_sq > dist_threshold) {
    // Far field nodes: Evaluate using moments.
    // Start with the monopole term r ⋅ (c_M₂ × c_M₁).
    double result = dr.dot(moment2.tangent.cross(moment1.tangent));
    // Declare rdist2 = |r|⁻², and multipler = 1/(4 π |r|³)
    double rdist2 = 1.0 / dist_sq;
    double multiplier = rdist2 / sqrt(dist_sq) * INV_FOUR_PI;
    // Now let's add the dipole and quadrupole terms and compute the error
    // estimate:
    {
      double eval_error_without_multiplier;
      result +=
          EvalFarFieldHOCorrection(moment1, moment2, dr, rdist2, radius1,
                                   radius2, eval_error_without_multiplier);
      eval_error_without_multiplier *= multiplier;

      eval_error += eval_error_without_multiplier;
    }
    // Return the multiplier times the three terms added up.
    return multiplier * result;
  } else {
    // Near-field nodes: Subdivide the larger node and traverse their children.
    int subdivider = obj1;
    int keep = obj2;
    int s_idx = node1;
    int k_idx = node2;
    const Eigen::ModifiedKdBVH<double, 3, SegmentIntEntry> *s_tree =
        &segment_tree1;
    const Eigen::ModifiedKdBVH<double, 3, SegmentIntEntry> *k_tree =
        &segment_tree2;
    const std::vector<Moment> *s_moments = &tree_moments1;
    const std::vector<Moment> *k_moments = &tree_moments2;

    if (obj1 < 0 && obj2 < 0) {
      // Subdivide the larger box.
      if (radius2sq > radius1sq) {
        std::swap(subdivider, keep);
        std::swap(s_idx, k_idx);
        std::swap(s_tree, k_tree);
        std::swap(s_moments, k_moments);
      }
    } else if (obj2 < 0) {
      // Subdivide obj2 because it can be subdivided.
      std::swap(subdivider, keep);
      std::swap(s_idx, k_idx);
      std::swap(s_tree, k_tree);
      std::swap(s_moments, k_moments);
    }

    // Subdivide and recurse on the children.
    std::pair<int, int> indices = s_tree->getChildrenIndices(s_idx);
    return EvalLinkingNumberRecurse(*k_tree, *s_tree, *k_moments, *s_moments,
                                    k_idx, indices.first, eval_error, epssq,
                                    itolsq) +
           EvalLinkingNumberRecurse(*k_tree, *s_tree, *k_moments, *s_moments,
                                    k_idx, indices.second, eval_error, epssq,
                                    itolsq);
  }
}

/// This is a helper to evaluate the linking numbers using Barnes–Hut, by
/// traversing the segment tree with a breadth-first pass. Any far-field and
/// leaf node pairs are evaluated and accumulated. For near-field node pairs, we
/// add the two children of the larger node paired with the node from the other
/// tree into a list of nodes for the next pass.
///
/// @param [in] segment_tree1 is the segment tree of curve 1.
/// @param [in] segment_tree2 is the segment tree of curve 2.
/// @param [in] tree_moments1 contains the moments for segment_tree1.
/// @param [in] tree_moments2 contains the moments for segment_tree2.
/// @param [in] nodes_1 is a list of current nodes in segment_tree1 to evaluate.
/// @param [in] nodes_2 is a list of current nodes in segment_tree2 to evaluate.
/// @param [in] epssq is the square of machine epsilon used to prevent
///     singularities when the distance between nodes is zero. This variable
///     name comes from [Burtscher and Pingali 2011].
/// @param [in] itolsq is the square of beta; this variable name also comes from
///     [Burtscher and Pingali 2011].
/// @param [in] parallel determines whether to parallelize this evaluation.
/// @param [inout] eval_error accumulates the next-order error estimate for from
///     this pass.
/// @param [out] nodes_1_next is the resulting list of nodes in segment_tree1 to
///     evaluate in the next pass.
/// @param [out] nodes_2_next is the resulting list of nodes in segment_tree2 to
///     evaluate in the next pass.
///
/// @returns the linking number accumulated from any far-field or leaf node
///     pairs in this pass.
double EvalLinkingNumberOneIter(
    const Eigen::ModifiedKdBVH<double, 3, SegmentIntEntry> &segment_tree1,
    const Eigen::ModifiedKdBVH<double, 3, SegmentIntEntry> &segment_tree2,
    const std::vector<Moment> &tree_moments1,
    const std::vector<Moment> &tree_moments2, const std::vector<int> &nodes_1,
    const std::vector<int> &nodes_2, std::vector<int> &nodes_1_next,
    std::vector<int> &nodes_2_next, double &eval_error, bool parallel,
    double epssq, double itolsq) {
  const int nisize = nodes_1.size();
  if (nodes_2.size() != nisize) {
    std::cout << "EvalLinkingNumberOneIter sizes mismatch." << std::endl;
    throw std::runtime_error("EvalLinkingNumberOneIter error");
  }
  if (nisize == 0) {
    return 0.0;
  }
  nodes_1_next.resize(nisize * 2);
  nodes_2_next.resize(nisize * 2);
  int next_inds = 0;
  const int moment_index_offset1 = segment_tree1.getBoxArraySize();
  const int moment_index_offset2 = segment_tree2.getBoxArraySize();
  double result = 0.;
  double local_eval_error = 0;

// Note that even though Appendix B.2 of [Qu and James 2021] says 1000 as the
// minimum size for parallelization, it is fine to just use 100 here.
#pragma omp parallel for reduction(+:result) reduction(+:local_eval_error) \
    if(nisize >= 100 && parallel)
  for (int ni = 0; ni < nisize; ++ni) {
    int node1 = nodes_1[ni];
    int node2 = nodes_2[ni];
    const int obj1 = node1 - moment_index_offset1;
    const int obj2 = node2 - moment_index_offset2;
    const Eigen::AlignedBox3d &bbox1 =
        (obj1 < 0) ? segment_tree1.getVolume(node1)
                   : segment_tree1.getObject(obj1).second;
    const Eigen::AlignedBox3d &bbox2 =
        (obj2 < 0) ? segment_tree2.getVolume(node2)
                   : segment_tree2.getObject(obj2).second;
    double radius1sq = 0.25 * (bbox1.max() - bbox1.min()).squaredNorm();
    double radius2sq = 0.25 * (bbox2.max() - bbox2.min()).squaredNorm();
    double radius1 = std::sqrt(radius1sq);
    double radius2 = std::sqrt(radius2sq);
    using std::sqrt;
    double r1plusr2sq = (radius1 + radius2) * (radius1 + radius2);
    const Moment &moment1 = tree_moments1[node1];
    const Moment &moment2 = tree_moments2[node2];

    // Notationwise, we say dr = r = r̃₂ - r̃₁.
    Eigen::Vector3d dr = moment2.midpoint - moment1.midpoint;
    double dist_sq = dr.squaredNorm() + epssq;
    double dist_threshold = r1plusr2sq * itolsq + epssq;
    if (obj1 >= 0 && obj2 >= 0) {
      // Leaf nodes: Use the one-arctan expression for two finite segments.
      double determinant = dr.dot(moment2.tangent.cross(moment1.tangent));
      Eigen::Vector3d alpha = dr + 0.5 * (-moment2.tangent + moment1.tangent);
      Eigen::Vector3d beta = dr + 0.5 * (-moment2.tangent - moment1.tangent);
      Eigen::Vector3d gamma = dr + 0.5 * (moment2.tangent - moment1.tangent);
      double la = alpha.norm();
      double lb = beta.norm();
      double lc = gamma.norm();
      double ac = alpha.dot(gamma);
      double div1 =
          la * lb * lc + alpha.dot(beta) * lc + ac * lb + beta.dot(gamma) * la;
      Eigen::Vector3d delta = dr + 0.5 * (moment2.tangent + moment1.tangent);
      double ld = delta.norm();
      double div2 = la * ld * lc + alpha.dot(delta) * lc + ac * ld +
                    delta.dot(gamma) * la;
      // x'
      using std::copysign;
      la = div1 * div2 - determinant * determinant;
      // ac is s'' = y' * determinant. We also want it negative when det = 0 and
      // div1, div2 are both negative.
      ac = (determinant == 0) ? copysign(1.0, div1) + copysign(1.0, div2)
                              : (div1 + div2); // = y' * determinant
      // y'
      lb = determinant * ac;
      lc = (determinant == 0) ? -copysign(1.0, div1)
                              : copysign(1.0, determinant);
      using std::atan2;
      result += INV_TWO_PI * atan2(lb, la) +
                (((ac < 0) || ((ac == 0) && (determinant < 0))) ? lc : 0.0);
    } else if (dist_sq > dist_threshold) {
      // Far field nodes: Evaluate using moments.
      // Start with the monopole term r ⋅ (c_M₂ × c_M₁).
      double local_result = dr.dot(moment2.tangent.cross(moment1.tangent));
      // Declare rdist2 = |r|⁻², and multipler = 1/(4 π |r|³)
      double rdist2 = 1.0 / dist_sq;
      double multiplier = rdist2 / sqrt(dist_sq) * INV_FOUR_PI;
      // Now let's add the dipole and quadrupole terms and compute the error
      // estimate:
      {
        double eval_error_without_multiplier;
        local_result +=
            EvalFarFieldHOCorrection(moment1, moment2, dr, rdist2, radius1,
                                     radius2, eval_error_without_multiplier);
        eval_error_without_multiplier *= multiplier;
        local_eval_error += eval_error_without_multiplier;
      }
      // Accumulate the multiplier times the three terms added up.
      result += multiplier * local_result;
    } else {
      // Near-field nodes: Subdivide the larger node and traverse their
      // children.
      bool subdivide_obj1 = true;

      if (obj1 < 0 && obj2 < 0) {
        // Subdivide the larger box.
        if (radius2sq > radius1sq) {
          subdivide_obj1 = false;
        }
      } else if (obj2 < 0) {
        // Subdivide obj2 because it can be subdivided.
        subdivide_obj1 = false;
      }

      // Subdivide and add the children to the next list for next pass.
      int half_ind = -1;
#ifdef _MSC_VER
#pragma omp critical
      {
        half_ind = next_inds++;
      }
#else
#pragma omp atomic capture
      half_ind = next_inds++;
#endif
      if (subdivide_obj1) {
        std::pair<int, int> indices = segment_tree1.getChildrenIndices(node1);
        nodes_1_next[2 * half_ind] = indices.first;
        nodes_1_next[2 * half_ind + 1] = indices.second;
        nodes_2_next[2 * half_ind] = node2;
        nodes_2_next[2 * half_ind + 1] = node2;
      } else {
        std::pair<int, int> indices = segment_tree2.getChildrenIndices(node2);
        nodes_1_next[2 * half_ind] = node1;
        nodes_1_next[2 * half_ind + 1] = node1;
        nodes_2_next[2 * half_ind] = indices.first;
        nodes_2_next[2 * half_ind + 1] = indices.second;
      }
    }
  }
  nodes_1_next.resize(2 * next_inds);
  nodes_2_next.resize(2 * next_inds);
  eval_error = local_eval_error;
  return result;
}

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
 * @param [in] epssq_in, the square of machine epsilon used to prevent
 *     singularities when the distance between nodes is zero. This variable name
 *     comes from [Burtscher and Pingali 2011].
 * @param [in] itolsq_in, the square of beta; this variable name also comes from
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
    const std::vector<Moment> &tree_moments2, const double epssq_in,
    const double itolsq_in, double &eval_errs_total,
    bool parallel /* = true */) {
  eval_errs_total = 0;
  double eval_error;
  if (parallel) {
    // To parallelize this while maintaining load balancing, we take a
    // breadth-first pass at a time and maintain a vector of node pairs to check
    // in the next pass. The node pairs are stored in nodes1 and nodes2. We keep
    // track of the number of pairs in each pass, as curr_size.
    // EvalLinkingNumberOneIter parallelizes the evaluation when curr_size is
    // greater than a lower threshold (set to 100, see
    // EvalLinkingNumberOneIter.). We also set an upper bound, because it gets
    // expensive to heap-allocate so many vectors if we use too much memory.
    // This is range is described in [Qu and James 2021] as "1000 to 500,000
    // node pairs" (Appendix B.2), but in the code we tweaked these values
    // slightly.
    std::vector<int> nodes1;
    std::vector<int> nodes2;
    std::vector<int> nextnodes1;
    std::vector<int> nextnodes2;
    double result = 0.;
    nextnodes1.push_back(segment_tree1.getRootIndex());
    nextnodes2.push_back(segment_tree2.getRootIndex());
    int curr_size = 1;
    // Set the upper bound to about 75k; if there are too many, then memory
    // swapping becomes an issue. While Appendix B.2 in [Qu and James 2021] says
    // 500k, 75k might be a better bound for smaller machines. (My machine, used
    // in [Qu and James 2021], has 64GB of RAM.)
    while (curr_size > 0 && curr_size < 75000) {
      // Set the current nodes to the nextnodes from last iteration.
      std::swap(nextnodes1, nodes1);
      std::swap(nextnodes2, nodes2);
      result += EvalLinkingNumberOneIter(
          segment_tree1, segment_tree2, tree_moments1, tree_moments2, nodes1,
          nodes2, nextnodes1, nextnodes2, eval_error, parallel, epssq_in,
          itolsq_in);
      curr_size = nextnodes1.size();
      eval_errs_total += eval_error;
    }
    eval_error = 0;
    // After we're done with the breadth-first parallellism, we just call the
    // recursive (depth-first) version in parallel, using each of the 75k+ node
    // pairs (nextnodes1[i], nextnodes2[i]) as the root node. Dynamic scheduling
    // helps slightly here.
#pragma omp parallel for reduction(+ : result, eval_error) schedule(dynamic)
    for (int i = 0; i < curr_size; ++i) {
      result += EvalLinkingNumberRecurse(
          segment_tree1, segment_tree2, tree_moments1, tree_moments2,
          nextnodes1[i], nextnodes2[i], eval_error, epssq_in, itolsq_in);
    }
    eval_errs_total += eval_error;
    return result;
  } else {
    // When this is single-threaded, just use the recursive (depth-first)
    // version, which is both memory- and compute- efficient.
    double eval_error = 0;
    float result = EvalLinkingNumberRecurse(
        segment_tree1, segment_tree2, tree_moments1, tree_moments2,
        segment_tree1.getRootIndex(), segment_tree2.getRootIndex(), eval_error,
        epssq_in, itolsq_in);
    eval_errs_total += eval_error;
    return result;
  }
}

} // namespace VerifyCurves