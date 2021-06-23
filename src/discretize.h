/**
 * Copyright (C) 2021 Ante Qu <antequ@cs.stanford.edu>.
 */

#pragma once

#include "curve.h"
#include <vector>

namespace VerifyCurves {

/**
 * This module computes pairwise linking numbers from collections of
 * discretized curves, based on Section 3.3 of [Qu and James 2021]:
 *
 * - [Qu and James 2021] Ante Qu and Doug L. James. 2021. Fast Linking Numbers
 *                       for Topology Verification of Loopy Structures. ACM
 *                       Trans. Graph. 40, 4, Article 106 (August 2021), 19
 *                       pages.
 */

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
           const std::vector<std::vector<int>> &potential_links_list);

} // namespace VerifyCurves