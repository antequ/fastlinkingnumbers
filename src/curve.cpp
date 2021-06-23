/**
 * Copyright (C) 2021 Ante Qu <antequ@cs.stanford.edu>.
 * Copyright (C) 2018 Jonathan Leaf
 * Some snippets from Jonathan Leaf, Xinru (Lucy) Hua, and Jonathan Kaldor;
 * all modified by Ante Qu.
 */

#include "curve.h"

#include <iostream>
#include <vector>

namespace VerifyCurves {

/// These helper functions retrieve the spline basis vector B(t) given the
/// parameter t.
///
/// @param [in] t, a real value within [0,1] of the point on the spline.
/// @param [out] B, the spline basis vector at t.
inline void get_b(double t, double *B) {
  /**
   * This subroutine was originally supplied by Jonathan Leaf
   * and modified by Ante Qu.
   */
  B[0] = ((-(1.0 / 6) * t + (1.0 / 2)) * t - (1.0 / 2)) * t + (1.0 / 6);
  B[1] = ((1.0 / 2) * t - 1.0) * (t * t) + (2.0 / 3);
  B[2] = ((-(1.0 / 2) * t + (1.0 / 2)) * t + (1.0 / 2)) * t + (1.0 / 6);
  B[3] = (1.0 / 6) * (t * t) * t;
}

inline void get_catmull_rom_b(double t, double *B) {
  B[0] = ((-0.5 * t + 1.0) * t - 0.5) * t;
  B[1] = (1.5 * t - 2.5) * (t * t) + 1.0;
  B[2] = ((-1.5 * t + 2.0) * t + 0.5) * t;
  B[3] = (0.5 * t - 0.5) * (t * t);
}

inline void get_bezier_b(double t, double *B) {
  B[0] = ((-1.0 * t + 3.0) * t - 3.0) * t + 1.0;
  B[1] = ((3.0 * t - 6.0) * t + 3.0) * t;
  B[2] = (-3.0 * t + 3.0) * (t * t);
  B[3] = 1.0 * t * t * t;
}

/**
 * The default constructor generates an empty, unlooped, B-Spline curve.
 */
Curve::Curve() {}

/**
 * Translate rigidly translates all points of the curve.
 *
 * @param [in] d is the translation vector.
 *
 * @post the underlying control points q_[i] are each replaced with q_[i] + d.
 */
void Curve::Translate(const Point3 &d) {
  for (auto it = q_.begin(); it != q_.end(); ++it) {
    (*it) += d;
  }
}

/**
 * Rotate rigidly rotates all points of the curve about the origin.
 *
 * @param [in] R is the rotation matrix
 *
 * @post the underlying control points q_[i] are each replaced with R q_[i].
 */
void Curve::Rotate(const Eigen::Matrix3d &R) {
  for (int i = 0; i < q_.size(); i++) {
    q_[i] = R * q_[i];
  }
}

/**
 * Scale rigidly scales all points of the curve about the origin.
 *
 * @param [in] s is the scaling vector
 *
 * @post the underlying control points q_[i] are each replaced with (q_[i].x
 *     s.x, q_[i].y s.y, q_[i].z s.z).
 */
void Curve::Scale(const Point3 &s) {
  for (auto it = q_.begin(); it != q_.end(); ++it) {
    (*it) = (*it).cwiseProduct(s);
  }
}

/**
 * Only for closed-loop curves: GetMidpoints(), GetTangents(),
 * GetStartPoints() return vectors of midpoints, "tangents" (which are
 * finishpoints minus startpoints), and startpoints, respectively.
 * GetLoopedStartPoints() returns a vector of startpoints, plus the last
 * endpoint.
 *
 * @returns a vector of the desired points.
 */
std::vector<Point3> Curve::GetMidpoints() const {
  if (!looped_)
    throw std::runtime_error(
        "Curve::GetMidpoints() does not support unlooped curves.");
  std::vector<Point3> midpoints(q_.size());
  for (int i = 0; i < q_.size(); ++i) {
    midpoints[i] = get_point(i, 0.5);
  }
  return midpoints;
}

std::vector<Point3> Curve::GetTangents() const {
  if (!looped_)
    throw std::runtime_error(
        "Curve::GetTangents() does not support unlooped curves.");

  std::vector<Point3> startpoints = GetStartPoints();
  const size_t n = startpoints.size();
  std::vector<Point3> tangents = std::vector<Point3>(n);
  for (int i = 0; i < n - 1; ++i) {
    tangents[i] = startpoints[i + 1] - startpoints[i];
  }
  tangents[n - 1] = startpoints[0] - startpoints[n - 1];
  return tangents;
}

std::vector<Point3> Curve::GetStartPoints() const {
  if (!looped_)
    throw std::runtime_error(
        "Curve::GetStartPoints() does not support unlooped curves.");
  std::vector<Point3> startpoints(q_.size());
  for (int i = 0; i < q_.size(); ++i) {
    startpoints[i] = get_point(i, 0.0);
  }
  return startpoints;
}

std::vector<Point3> Curve::GetLoopedStartPoints() const {
  if (!looped_)
    throw std::runtime_error(
        "Curve::GetLoopedStartPoints() does not support unlooped curve.");
  std::vector<Point3> startpoints = GetStartPoints();
  startpoints.push_back(startpoints[0]);
  return startpoints;
}

/**
 * This method grabs a specific point on a specific spline segment, using the
 * underlying parameter t.
 *
 * @param [in] seg_index the index of the specific segment, under 0-based
 *     indexing.
 * @param [in] t the parameter, a real value within [0,1] of the point on the
 *     segment.
 *
 * @returns the desired point.
 */
Point3 Curve::get_point(int seg_index, double t) const {
  if (t < 0.0 || t > 1.0 || seg_index < 0 || seg_index >= num_segments()) {
    std::cerr << "Curve::get_point() has an invalid i or t: i: " << seg_index
              << ", t: " << t << ", n: " << q_.size() << "." << std::endl;
    throw std::runtime_error("Curve::get_point() has an invalid i or t.");
  }
  if (segment_type_ == SegmentType::Polyline) {
    const Point3 &start = q_[seg_index];
    const Point3 &end = q_[(seg_index + 1) % q_.size()];
    return (1.0 - t) * start + t * end;
  } else {
    // std::vector<const Point3&> q_4s = std::vector<Point3>(4);
    const Point3 &q_0 = q_[seg_index];
    const Point3 &q_1 = q_[(seg_index + 1) % q_.size()];
    const Point3 &q_2 = q_[(seg_index + 2) % q_.size()];
    const Point3 &q_3 = q_[(seg_index + 3) % q_.size()];
    double b[4] = {0};
    if (segment_type_ == SegmentType::BSpline) {
      get_b(t, b);
    } else if (segment_type_ == SegmentType::CatmullRom) {
      get_catmull_rom_b(t, b);
    } else if (segment_type_ == SegmentType::Bezier) {
      get_bezier_b(t, b);
    } else {
      throw std::runtime_error("Invalid segment type in Curve::get_point().");
    }
    return b[0] * q_0 + b[1] * q_1 + b[2] * q_2 + b[3] * q_3;
  }
}

/// This helper function computes the minimum and maximum of a spline
/// coordinate for a subspline, given its control point values (q_1 to q_4) and
/// endpoint parameter values (t0 and t1).
///
/// @param [in] seg_type, the spline segment type,
/// @param [in] q_1–q_4, the coordinate at control points 1 through 4,
/// @param [in] t0, the t parameter for the low end of the subspline,
/// @param [in] t1, the t parameter for the high end of the subspline,
///
/// @param [out] min, the minimum coordinate for this subspline,
/// @param [out] max, the maximum coordinate for this subspline.
inline void compute_box_1d(SegmentType seg_type, double q_1, double q_2,
                           double q_3, double q_4, double t0, double t1,
                           double &min, double &max) {
  /**
   * Parts of this subroutine were originally supplied by Jonathan Kaldor and
   * modified by Xinru Hua and Ante Qu.
   */
  double x0, x1, x2, x3, t;
  if (seg_type == SegmentType::BSpline) {
    // t^3 term
    x0 = -(1.0 / 6.0) * q_1 + 0.5 * q_2 - 0.5 * q_3 + (1.0 / 6.0) * q_4;
    // t^2 term
    x1 = 0.5 * q_1 - q_2 + 0.5 * q_3;
    // t^1 term
    x2 = -0.5 * q_1 + 0.5 * q_3;
    // constant term
    x3 = (1.0 / 6.0) * q_1 + (2.0 / 3.0) * q_2 + (1.0 / 6.0) * q_3;
  } else if (seg_type == SegmentType::CatmullRom) {
    // t^3 term
    x0 = -0.5 * q_1 + 1.5 * q_2 - 1.5 * q_3 + 0.5 * q_4;
    // t^2 term
    x1 = 1.0 * q_1 - 2.5 * q_2 + 2.0 * q_3 - 0.5 * q_4;
    // t^1 term
    x2 = -0.5 * q_1 + 0.5 * q_3;
    // constant term
    x3 = 1.0 * q_2;
  } else {
    // Bezier curve
    // t^3 term
    x0 = -1.0 * q_1 + 3.0 * q_2 - 3.0 * q_3 + 1.0 * q_4;
    // t^2 term
    x1 = 3.0 * q_1 - 6.0 * q_2 + 3.0 * q_3;
    // t^1 term
    x2 = -3.0 * q_1 + 3.0 * q_2;
    // constant term
    x3 = 1.0 * q_1;
  }
  // one endpoint
  min = ((x0 * t0 + x1) * t0 + x2) * t0 + x3;
  // other endpoint
  max = ((x0 * t1 + x1) * t1 + x2) * t1 + x3;
  if (max < min) {
    std::swap(min, max);
  }
  double deltax = 4 * x1 * x1 - 12 * x2 * x0;
  if (deltax == 0) { // one root
    t = -2 * x1 / (6 * x0);
    if (t >= t0 && t <= t1) {
      double v = ((x0 * t + x1) * t + x2) * t + x3;
      min = std::min(min, v);
      max = std::max(max, v);
    }
  } else { // d > 0 - two real roots
    double radical;

    if (x1 >= 0) {
      radical = -2 * x1 - sqrt(deltax);
      t = radical / (6 * x0);
      if (t >= t0 && t <= t1) {
        double v = ((x0 * t + x1) * t + x2) * t + x3;
        min = std::min(min, v);
        max = std::max(max, v);
      }
    } else {
      radical = -2 * x1 + sqrt(deltax);
      t = radical / (6 * x0);
      if (t >= t0 && t <= t1) {
        double v = ((x0 * t + x1) * t + x2) * t + x3;
        min = std::min(min, v);
        max = std::max(max, v);
      }
    }

    t = 2 * x2 / radical;
    if (t >= t0 && t <= t1) {
      double v = ((x0 * t + x1) * t + x2) * t + x3;
      min = std::min(min, v);
      max = std::max(max, v);
    }
  }
}

/**
 * This method computes the bounding box of a subsegment of a spline segment,
 * with padding if desired.
 *
 * @param [in] seg_index the index of the specific segment, under 0-based
 *     indexing.
 * @param [in] t_lo the parameter of the low endpoint, a real value within
 *     [0,1] of the endpoint on the segment.
 * @param [in] t_hi the parameter of the high endpoint, a real value within
 *     [0,1] of the endpoint on the segment. t_hi must >= t_lo.
 * @param [in] pad the box is dilated by this distance in every direction.
 *
 * @returns the bounding box.
 */
Box3d Curve::get_segment_bounding_box(int seg_index, double t_lo, double t_hi,
                                      double pad) const {
  Point3 min, max;
  if (seg_index < 0 || seg_index >= num_segments() || t_lo < 0 || t_hi > 1 ||
      t_lo > t_hi) {
    std::cerr << "Curve::get_segment_bounding_box() has an invalid i or t: i: "
              << seg_index << ", tlo: " << t_lo << ", thi: " << t_hi
              << ", n: " << q_.size() << "." << std::endl;
    throw std::runtime_error(
        "Curve::get_segment_bounding_box() has an invalid i or t.");
  }
  if (segment_type_ == SegmentType::Polyline) {
    const Point3 lo_pt = get_point(seg_index, t_lo);
    const Point3 hi_pt = get_point(seg_index, t_hi);
    min = lo_pt.cwiseMin(hi_pt);
    max = lo_pt.cwiseMax(hi_pt);
  } else {
    const Point3 &q_0 = q_[seg_index];
    const Point3 &q_1 = q_[(seg_index + 1) % q_.size()];
    const Point3 &q_2 = q_[(seg_index + 2) % q_.size()];
    const Point3 &q_3 = q_[(seg_index + 3) % q_.size()];
    compute_box_1d(segment_type_, q_0.x(), q_1.x(), q_2.x(), q_3.x(), t_lo,
                   t_hi, min.x(), max.x());
    compute_box_1d(segment_type_, q_0.y(), q_1.y(), q_2.y(), q_3.y(), t_lo,
                   t_hi, min.y(), max.y());
    compute_box_1d(segment_type_, q_0.z(), q_1.z(), q_2.z(), q_3.z(), t_lo,
                   t_hi, min.z(), max.z());
  }
  min.array() -= pad;
  max.array() += pad;
  return Box3d(min, max);
}

/**
 * This method precomputes the bounding box of the entire curve and caches it.
 *
 * @post stores the curve bounding box in bounding_box_.
 */
void Curve::ComputeCurveBoundingBox() const {
  int ulimit = num_segments();
  if (ulimit == 0) {
    bounding_box_ = Box3d();
  } else {
    // Don't start with a null box. Start with the first box.
    bounding_box_ = get_segment_bounding_box(0);
    for (int i = 1; i < ulimit; ++i) {
      bounding_box_.extend(get_segment_bounding_box(i));
    }
  }
  bounding_box_initialized_ = true;
}

} // namespace VerifyCurves