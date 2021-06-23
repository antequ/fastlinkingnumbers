/**
 * Copyright (C) 2021 Ante Qu <antequ@cs.stanford.edu>.
 */

#include "discretize.h"
#include "BVH.h" // This is Eigen's BVH module.
#include "curve.h"

#include <iostream>
#include <unordered_set>
#include <utility>
#include <vector>

#define SQRT_OF_THREE (1.73205080757)

namespace Eigen {
template <typename T>
inline VerifyCurves::Box3d
bounding_box(const std::pair<T, VerifyCurves::Box3d> &spline) {
  return spline.second;
}
} // namespace Eigen

namespace VerifyCurves {
// SplineAddress is (curve index, segment index, t-low, t-high) where t ranges
// from 0 to 1 for a spline segment.
typedef std::tuple<int, int, double, double> SplineAddress;

// SplineEntry is (SplineAddress, Spline Bounding Box), which is useful for BVH
// traversal.
typedef std::pair<SplineAddress, Box3d> SplineEntry;

// Hash specialization for std::hash based on
// https://www.techiedelight.com/use-std-pair-key-std-unordered_map-cpp/
struct spline_address_hash {
  std::size_t operator()(const SplineAddress &addr) const {

    return std::hash<int>()(std::get<0>(addr)) ^
           std::hash<int>()(std::get<1>(addr)) ^
           std::hash<double>()(std::get<2>(addr)) ^
           std::hash<double>()(std::get<3>(addr));
  }
};

// This is an Intersector to use with Eigen's BVH module. We run this for a
// segment tree of a single curve against the tree of another curve. We traverse
// further when two boxes overlap, and at the leaf level, we mark segments from
// the first curve, into "results", as segments that need subdivision after this
// pass.
struct BoxBoxIntersectorForSegments {
  typedef SplineEntry EntryType;
  BoxBoxIntersectorForSegments() : npairs(0), results_list() {}
  typedef double Scalar;
  bool intersectVolumeVolume(const Box3d &r1, const Box3d &r2) {
    return r1.intersects(r2);
  }
  bool intersectVolumeObject(const Box3d &r, const EntryType &s) {
    return r.intersects(s.second);
  }
  bool intersectObjectVolume(const EntryType &s, const Box3d &r) {
    return r.intersects(s.second);
  }
  // If v2 overlaps with v1, then we accumulate thev1 entry into a result list
  // for this curve, indicating that it needs to be subdivided.
  bool intersectObjectObject(const EntryType &v1, const EntryType &v2) {
    // They must be from different curves. (Should already be true.)
    if ((std::get<0>(v1.first) != std::get<0>(v2.first)) &&
        v1.second.intersects(v2.second)) {
      ++npairs;
      results_list.insert(v1.first);
    }
    return false;
  }
  void reset() {
    npairs = 0;
    results_list.clear();
  }
  // results_list marks the spline segments that need subdivision because they
  //     overlap with another curve's splines.
  std::unordered_set<SplineAddress, spline_address_hash> results_list;
  // npairs keeps track of the number of overlaps that segments of this curve
  //     have with segments of other curves. Note that this is bigger than
  //     results_list.size() because a single segment can overlap with multiple
  //     other segments.
  int npairs;
};

/**
 * This function, Discretize, discretizes curves into line segments.
 *
 * @param [in] curves, a list of curves.
 * @param [in] potential_links_list, a list of potential links between curves
 *     (curves that have overlapping bounding boxes). potential_links_list must
 *     be represented as a vector of vectors, with the first vector indexed by
 *     curve id, and the second vector storing a list of ALL curves the first
 *     curve overlaps with. Each potentially linked looped pair is therefore
 *     represented twice.
 *
 * @returns It returns a list of discretized curves, with their new curve
 *     bounding boxes already precomputed and cached.
 */
std::vector<Curve>
Discretize(const std::vector<Curve> &curves,
           const std::vector<std::vector<int>> &potential_links_list) {
  const int ncurves = curves.size();
  const double eps = std::numeric_limits<double>::epsilon();

  // sum up number of entries
  int entrycount = 0;
  std::vector<std::vector<SplineEntry>> entry_trees(ncurves);
#ifdef _MSC_VER
  std::vector<double> min_box_diags(ncurves,
                                    std::numeric_limits<double>::infinity());
#else
  double min_box_diag = std::numeric_limits<double>::infinity();
#endif
  int min_curve_index, min_box_index;
  double avg_coordinate_size = 0.0;
  using std::abs;
  using std::min;
#ifdef _MSC_VER
#pragma omp parallel for reduction(+ : entrycount, avg_coordinate_size)
  for (int i = 0; i < ncurves; ++i) {
    double min_box_diag = min_box_diags[i];
#else
#pragma omp parallel for reduction(+:entrycount, avg_coordinate_size) \
    reduction(min:min_box_diag)
  for (int i = 0; i < ncurves; ++i) {
#endif
    entry_trees[i] = std::vector<SplineEntry>(curves[i].num_segments());
    for (int j = 0; j < curves[i].num_segments(); ++j) {
      const Box3d &bbox = curves[i].get_segment_bounding_box(j);
      entry_trees[i][j] = SplineEntry(SplineAddress(i, j, 0.0, 1.0), bbox);
      avg_coordinate_size += curves[i].get_q()[j].norm();
      double bboxsize = bbox.sizes().norm();
      if (!(bboxsize >= min_box_diag)) {
        min_box_diag = bboxsize;
        min_curve_index = i;
        min_box_index = j;
      }
    }
    entrycount += curves[i].num_segments();
  }
  avg_coordinate_size /= entrycount * SQRT_OF_THREE;
  const double box_diag_threshold = eps * avg_coordinate_size;
  // use inversion in case min_box_diag is nan.
#ifdef _MSC_VER
  for (int i = 0; i < ncurves; ++i) {
    double min_box_diag = min_box_diags[i];
#endif
    if (!(min_box_diag >= box_diag_threshold)) {
      std::cerr << "Discretize: Smallest input segment, (" << min_curve_index
                << "," << min_box_index << ") has a size of " << min_box_diag
                << ", which is smaller than eps * avg coordinate. Consider "
                   "combining them into one vertex."
                << std::endl;
      throw std::runtime_error(
          "Discretize: Smallest input segment is smaller "
          "than eps * avg coordinate. Consider combining.");
    }
#ifdef _MSC_VER
  }
#endif

  BoxBoxIntersectorForSegments intersector;

  std::vector<std::unordered_set<SplineAddress, spline_address_hash>>
      collided_splines(ncurves);

  std::vector<std::vector<SplineEntry>> collided_splines_trees(ncurves);
  std::vector<std::vector<SplineAddress>> uncollided_splines(ncurves);

  std::vector<Eigen::ModifiedKdBVH<double, 3, SplineEntry>> trees(ncurves);
  constexpr int MAX_RECURSION = 55;
  int npairs = 0;
  for (int iteration = 0; iteration < MAX_RECURSION; ++iteration) {
// Build trees.
#pragma omp parallel for
    for (int i = 0; i < ncurves; ++i) {
      trees[i] = Eigen::ModifiedKdBVH<double, 3, SplineEntry>(
          entry_trees[i].begin(), entry_trees[i].end());
    }
    npairs = 0;
// intersect curves
#pragma omp parallel for reduction(+ : npairs) private(intersector)
    for (int i = 0; i < ncurves; ++i) {
      double min_box_diag = std::numeric_limits<double>::infinity();
      intersector.reset();
      for (int j : potential_links_list[i]) {
        Eigen::BVIntersect(trees[i], trees[j], intersector);
      }
      for (const SplineAddress &result : intersector.results_list) {
        int x1, x2;
        double xtl, xth;
        x1 = std::get<0>(result);  // curve index
        x2 = std::get<1>(result);  // segment index
        xtl = std::get<2>(result); // t of start point
        xth = std::get<3>(result); // t of end point
        double xmid = 0.5 * (xtl + xth);
        Box3d bbox1 = curves[x1].get_segment_bounding_box(x2, xtl, xmid);
        Box3d bbox2 = curves[x1].get_segment_bounding_box(x2, xmid, xth);
        collided_splines_trees[i].push_back(
            SplineEntry(SplineAddress(x1, x2, xtl, xmid), bbox1));
        collided_splines_trees[i].push_back(
            SplineEntry(SplineAddress(x1, x2, xmid, xth), bbox2));
        collided_splines[i].insert(result);
        min_box_diag =
            min(min_box_diag, min(bbox1.sizes().norm(), bbox2.sizes().norm()));
      }
      // append to uncollided splines list
      for (SplineEntry &entry : entry_trees[i]) {
        if (collided_splines[i].find(entry.first) ==
            collided_splines[i].end()) {
          uncollided_splines[i].push_back(entry.first);
        }
      }
      npairs += intersector.npairs;
      // Use this expression, which is also triggered if min_box_diag is NaN.
      if (!(min_box_diag >= box_diag_threshold)) {
#pragma omp critical
        {
          std::cout << "Discretize: Smallest segment is smaller than eps * avg "
                       "coordinate. Perhaps curves are too close."
                    << std::endl;
          throw std::runtime_error("Discretize: Smallest segment is smaller "
                                   "than eps * avg coordinate. "
                                   "Perhaps curves are too close.");
        }
      }
    }
    if ((npairs % 2) != 0) {
      std::cout << "Error: npairs is odd: " << npairs << "." << std::endl;
    }
    std::cout << "Iteration " << iteration + 1 << ", " << npairs / 2
              << " overlaps." << std::endl;

    if (npairs == 0) {
      std::cout << "Finished after " << iteration + 1 << " iteration(s)."
                << std::endl;
      break;
    }
    std::swap(entry_trees, collided_splines_trees);
    for (int i = 0; i < ncurves; ++i) {
      collided_splines_trees[i].clear();
      collided_splines[i].clear();
    }
  }
  if (npairs > 0) {
    std::cout << "Max recursion (curves are too close)." << std::endl;
    std::cout << "Number of overlaps: " << npairs / 2 << std::endl;
    throw std::runtime_error(
        "Reached max recursion (curves are too close) in Discretize().");
  }
  Curve segment_curve;
  segment_curve.set_looped(true);
  segment_curve.set_segment_type(SegmentType::Polyline);
  std::vector<Curve> result(ncurves, segment_curve);
  int discretized_count = 0;
#pragma omp parallel for reduction(+ : discretized_count)
  for (int y_ind = 0; y_ind < ncurves; ++y_ind) {
    std::sort(uncollided_splines[y_ind].begin(),
              uncollided_splines[y_ind].end());
    for (const SplineAddress &spline : uncollided_splines[y_ind]) {
      // just push all the minimum endpoints. It's a loop.
      result[y_ind].get_q().push_back(
          curves[y_ind].get_point(std::get<1>(spline), std::get<2>(spline)));
    }
    if (!curves[y_ind].is_looped()) {
      const SplineAddress &spline = uncollided_splines[y_ind].back();
      // if it's not a loop, add the maximum endpoint of the last segment.
      result[y_ind].get_q().push_back(
          curves[y_ind].get_point(std::get<1>(spline), std::get<3>(spline)));
      result[y_ind].set_looped(false);
    }
    result[y_ind].ComputeCurveBoundingBox();
    discretized_count += uncollided_splines[y_ind].size();
  }
  std::cout << "Original segment count: " << entrycount
            << ", discretized count: " << discretized_count << std::endl;
  return result;
}

} // namespace VerifyCurves
