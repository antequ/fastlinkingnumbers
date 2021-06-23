/**
 * Copyright (C) 2021 Ante Qu <antequ@cs.stanford.edu>.
 */

#pragma once

#include "curve.h"
#include <vector>

namespace VerifyCurves {
/**
 * This module computes a list of potentially linked loop pairs from a list of
 * curves, based on Section 3.2 of [Qu and James 2021]:
 *
 * - [Qu and James 2021] Ante Qu and Doug L. James. 2021. Fast Linking Numbers
 *                       for Topology Verification of Loopy Structures. ACM
 *                       Trans. Graph. 40, 4, Article 106 (August 2021), 19
 *                       pages.
 */

/**
 * This function, PrecomputeCurveBoundingBoxes, precomputes the bounding box of
 * each curve by calling Curve::ComputeCurveBoundingBox().
 *
 * @param [inout] curves, a list of curves.
 */
void PrecomputeCurveBoundingBoxes(const std::vector<Curve> &curves);

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
GetPotentialLinksUniqueList(const std::vector<Curve> &curves);

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
std::vector<std::vector<int>> GetPotentialLinksUniqueListPerCurve(
    const std::vector<std::pair<int, int>> &unique_list, const int ncurves);

std::vector<std::vector<int>> GetPotentialLinksListPerCurve(
    const std::vector<std::pair<int, int>> &unique_list, const int ncurves);

/**
 * These wrapper functions combine the above two operations into one, taking the
 * list of curves as input and outputting the vector of vectors of potential
 * links. GetPotentialLinksUniqueListPerCurve does not repeat pairs and only
 * puts each pair in the list for the lower-indexed curve.
 * GetPotentialLinksListPerCurve repeats each pair by putting it into the lists
 * of both curves that overlap.
 *
 * @param [in] curves, a list of curves with precomputed bounding boxes.
 *
 * @returns a vector of vectors of potential links; the out vector indexes over
 *     curves, and the inner vector lists all curves that overlap with this
 *     curve.
 */
inline std::vector<std::vector<int>>
GetPotentialLinksListPerCurve(const std::vector<Curve> &curves) {
  return GetPotentialLinksListPerCurve(GetPotentialLinksUniqueList(curves),
                                       curves.size());
}

inline std::vector<std::vector<int>>
GetPotentialLinksUniqueListPerCurve(const std::vector<Curve> &curves) {
  return GetPotentialLinksUniqueListPerCurve(
      GetPotentialLinksUniqueList(curves), curves.size());
}

} // namespace VerifyCurves