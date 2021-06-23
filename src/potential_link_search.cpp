/**
 * Copyright (C) 2021 Ante Qu <antequ@cs.stanford.edu>.
 */

#include "potential_link_search.h"
#include "BVH.h"
#include "curve.h"

#include <parallel/algorithm>
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

template <typename T> struct BoxBoxIntersectorCurve {
  typedef std::pair<T, Box3d> EntryType;
  BoxBoxIntersectorCurve(int ncurves)
      : npairs(0), results(), ncurvesminusone(ncurves - 1) {}
  typedef double Scalar;
  // returns true if product of volumes intersects the query
  bool intersectVolumeVolume(const Box3d &r1, const Box3d &r2) {
    return r1.intersects(r2);
  }
  // returns true if the volume-object product intersects the query
  bool intersectVolumeObject(const Box3d &r, const EntryType &s) {
    return r.intersects(s.second);
  }
  // returns true if the volume-object product intersects the query
  bool intersectObjectVolume(const EntryType &s, const Box3d &r) {
    return r.intersects(s.second);
  }
  // returns true if the search should terminate immediately
  bool intersectObjectObject(const EntryType &v1, const EntryType &v2) {
    // They must be from different curves.
    if ((v1.first < v2.first) && v1.second.intersects(v2.second)) {
      ++npairs;
      results.emplace_back(v1.first, v2.first);
    }
    return false;
  }
  void reset() {
    npairs = 0;
    results.clear();
  }
  std::vector<std::pair<T, T>> results;
  int npairs;
  const int ncurvesminusone;
};

/**
 * This function, PrecomputeCurveBoundingBoxes, precomputes the bounding box of
 * each curve by calling Curve::ComputeCurveBoundingBox().
 *
 * @param [inout] curves, a list of curves.
 */
void PrecomputeCurveBoundingBoxes(const std::vector<Curve> &curves) {
#pragma omp parallel for
  for (int i = 0; i < curves.size(); ++i) {
    curves[i].ComputeCurveBoundingBox();
  }
}

typedef std::pair<int, Box3d> CurveEntry;

/**
 * This function, GetPotentialLinksUniqueList, computes a list of potentially
 * linked loop pairs from the curves by intersecting their bounding boxes in a
 * KdBVH. The out list consists of all unique overlaps.
 *
 * @param [in] curves, a list of curves with precomputed bounding boxes.
 *
 * @returns a list of pairs of curves, each pair representing one overlap. The
 *     pairs are ordered in that their lower index is the first element, and
 *     upper index is the second, and they are sorted.
 */
std::vector<std::pair<int, int>>
GetPotentialLinksUniqueList(const std::vector<Curve> &curves) {
  const int ncurves = curves.size();
  constexpr int bvh_threshold = 50;
  constexpr int parallel_threshold = 100000;
  if (ncurves >= bvh_threshold) {
    std::vector<CurveEntry> tree_list(ncurves);
    for (int i = 0; i < ncurves; ++i) {
      tree_list[i] = {i, curves[i].get_last_curve_bounding_box()};
    }
    Eigen::ModifiedKdBVH<double, 3, CurveEntry> tree(tree_list.begin(),
                                                     tree_list.end());
    BoxBoxIntersectorCurve<int> curve_intersector(ncurves);
    Eigen::BVIntersect(tree, tree, curve_intersector);
    if (curve_intersector.results.size() > parallel_threshold) {
      __gnu_parallel::sort(curve_intersector.results.begin(),
                           curve_intersector.results.end());
    } else {
      std::sort(curve_intersector.results.begin(),
                curve_intersector.results.end());
    }
    return curve_intersector.results;
  } else {
    std::vector<std::pair<int, int>> results;
    for (int i = 0; i < ncurves; ++i) {
      for (int j = i + 1; j < ncurves; ++j) {
        if (curves[i].get_last_curve_bounding_box().intersects(
                curves[j].get_last_curve_bounding_box())) {
          results.emplace_back(i, j);
        }
      }
    }
    return results;
  }
}

/**
 * These helper functions take the list of pairs and break them down into a list
 * per curve. GetPotentialLinksUniqueListPerCurve does not repeat pairs and only
 * puts each pair in the list for the lower-indexed curve.
 * GetPotentialLinksListPerCurve repeats each pair by putting it into the lists
 * of both curves that overlap.
 *
 * @param [in] unique_list, the list of potentially linked loop pairs from
 *     GetPotentialLinksUniqueList.
 * @param [in] ncurves, the number of curves in the model.
 *
 * @returns a vector of vectors of potential links; the out vector indexes over
 *     curves, and the inner vector lists all curves that overlap with this
 *     curve.
 */
std::vector<std::vector<int>> GetPotentialLinksListPerCurve(
    const std::vector<std::pair<int, int>> &unique_list, const int ncurves) {
  std::vector<std::vector<int>> results(ncurves);
  for (int i = 0; i < unique_list.size(); ++i) {
    results[unique_list[i].first].push_back(unique_list[i].second);
    results[unique_list[i].second].push_back(unique_list[i].first);
  }
#pragma omp parallel for
  for (int i = 0; i < ncurves; ++i) {
    std::sort(results[i].begin(), results[i].end());
  }
  return results;
}

std::vector<std::vector<int>> GetPotentialLinksUniqueListPerCurve(
    const std::vector<std::pair<int, int>> &unique_list, const int ncurves) {
  std::vector<std::vector<int>> results(ncurves);
  for (int i = 0; i < unique_list.size(); ++i) {
    results[unique_list[i].first].push_back(unique_list[i].second);
  }
  return results;
}

} // namespace VerifyCurves
