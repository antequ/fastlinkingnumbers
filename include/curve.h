/**
 * Copyright (C) 2021 Ante Qu <antequ@cs.stanford.edu>.
 */

#pragma once

#include <Eigen/StdVector>
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Vector3d)
#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#define PI                                                                     \
  (3.1415926535897932384626433832795028841971693993751058209749445923078164)

namespace VerifyCurves {
/**
 * A Curve stores a series of control points, either for a cubic spline or a
 * polyline curve. It supports B-Splines, uniform Catmull–Rom splines, Bezier
 * curves, and polylines. This class also comes with helper methods to
 * retrieve midpoints, specific points on spline segments, and bounding boxes
 * for splines and subsplines.
 *
 * It also has a helper method to compute and cache a bounding box for the
 * whole curve, as a mutable value. ComputeCurveBoundingBox() must be called to
 * precompute the bounding box of the entire curve before
 * get_last_curve_bounding_box() can be called. Both are const methods that can
 * be called on const Curves, as they do not change the underlying curve.
 *
 * In addition to the control points, this class also stores its segment type
 * (segment_type_) as well as whether it should be interpreted as a closed loop
 * (looped_).
 */

enum class SegmentType { BSpline, CatmullRom, Bezier, Polyline };

typedef Eigen::Vector3d Point3;
typedef Eigen::AlignedBox<double, 3> Box3d;

// Stores the data for a single curve.
class Curve {

public:
  /**
   * The default constructor generates an empty, unlooped, B-Spline curve.
   */
  Curve();

  /**
   * This method returns the number of degrees of freedom of the curve.
   *
   * @returns the number of coordinates stored in the control points in q_.
   */
  int num_dofs() const { return q_.size() * 3; }

  /**
   * This method returns the number of spline segments this curve contains.
   *
   * @returns the number of spline segments the control points in q_ represent.
   */
  int num_segments() const {
    return (looped_)
               ? q_.size()
               : ((segment_type_ == SegmentType::Polyline) ? (q_.size() - 1)
                                                           : (q_.size() - 3));
  }

  /**
   * Translate rigidly translates all points of the curve.
   *
   * @param [in] d is the translation vector.
   *
   * @post the underlying control points q_[i] are each replaced with q_[i] + d.
   */
  void Translate(const Point3 &d);

  /**
   * Rotate rigidly rotates all points of the curve about the origin.
   *
   * @param [in] R is the rotation matrix
   *
   * @post the underlying control points q_[i] are each replaced with R q_[i].
   */
  void Rotate(const Eigen::Matrix3d &R);

  /**
   * Scale rigidly scales all points of the curve about the origin.
   *
   * @param [in] s is the scaling vector
   *
   * @post the underlying control points q_[i] are each replaced with (q_[i].x
   *     s.x, q_[i].y s.y, q_[i].z s.z).
   */
  void Scale(const Point3 &s);

  /**
   * Only for closed-loop curves: GetMidpoints(), GetTangents(),
   * GetStartPoints() return vectors of midpoints, "tangents" (which are
   * finishpoints minus startpoints), and startpoints, respectively.
   * GetLoopedStartPoints() returns a vector of startpoints, plus the last
   * endpoint.
   *
   * @returns a vector of the desired points.
   */
  std::vector<Point3> GetMidpoints() const;
  std::vector<Point3> GetTangents() const;
  std::vector<Point3> GetStartPoints() const;
  std::vector<Point3> GetLoopedStartPoints() const;

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
  Point3 get_point(int seg_index, double t) const;

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
  Box3d get_segment_bounding_box(int seg_index, double t_lo = 0.0,
                                 double t_hi = 1.0, double pad = 0.0) const;

  /**
   * This method precomputes the bounding box of the entire curve and caches it.
   *
   * @post stores the curve bounding box in bounding_box_.
   */
  void ComputeCurveBoundingBox() const;

  /**
   * Retrieves the bounding box of the entire curve. ComputeCurveBoundingBox()
   * must have been called before this to get a valid bounding box.
   *
   * @returns the cached bounding box of the entire curve.
   */
  Box3d get_last_curve_bounding_box() const {
    if (!bounding_box_initialized_)
      throw std::runtime_error("Bounding box uninitialized. Use "
                               "ComputeCurveBoundingBox() to initialize it.");
    return bounding_box_;
  }

  /**
   * These two methods get the list of control points, either as a modifiable or
   * const reference.
   */
  std::vector<Point3> &get_q() { return q_; }
  const std::vector<Point3> &get_q() const { return q_; }

  /**
   * These two methods get and set whether the curve is a closed loop.
   */
  bool is_looped() const { return looped_; }
  void set_looped(bool looped) { looped_ = looped; }

  /**
   * These two methods get and set the segment spline type. This value indicates
   * whether the control points are B-Splines, uniform Catmull–Rom splines,
   * Bezier curves, or polylines.
   */
  SegmentType get_segment_type() const { return segment_type_; }
  void set_segment_type(SegmentType type) { segment_type_ = type; }

private:
  // q_ is the list of control points.
  std::vector<Point3> q_;

  // This value indicates whether the control points are B-Splines, uniform
  // Catmull–Rom splines, Bezier curves, or polylines.
  SegmentType segment_type_ = SegmentType::BSpline;

  // This indicates whether this curve is a closed loop.
  bool looped_ = false;

  // These values are used to cache the curve bounding box.
  mutable Box3d bounding_box_;
  mutable bool bounding_box_initialized_ = false;
};

} // namespace VerifyCurves