#pragma once
/**
 * Copyright 2020 Ante Qu.  All rights reserved.
 */

// In this file, every instance of "Yarn" simply refers to a "Curve". This code is distributed as-is.
// "helper_math.h" is from Nvidia, please find it here for example:
// https://github.com/NVIDIA/cuda-samples/blob/master/Common/helper_math.h

#include "helper_math.h"
#include <vector>

struct YarnTree {
  // device pointers
  float4 *midpts;
  float4 *tangents;
  int *sortl;
  int *childl;
  // relevant tree metadata
  float radius;
  int nnodes;
};

struct YarnDipoleTree {
  // device pointers
  float4 *midpts;
  float4 *tangents;
  float4 *dipole1s;
  float4 *dipole3s;
  int *sortl;
  int *childl;
  int* countl;
  int* startl;
  // relevant tree metadata
  float radius;
  int nnodes;
};

// Barnes Hut, using only monopole terms
float ComputeLinkingNumberBH(const YarnTree &source, int nbodiess,
                             const YarnTree &target, int nbodiest,
                             float epssq = 0.0025f, float itolsq = 4.f);

// Barnes Hut, using dipole terms
float ComputeLinkingNumberBHDipole(const YarnDipoleTree &source, int nbodiess,
                                   const YarnDipoleTree &target, int nbodiest,
                                   const int blocks, float epssq = 0.0025f,
                                   float itolsq = 4.f, float distratiosq = 4.8f);

// Experimental function for computing the writhe instead of the linking number.
float ComputeWritheBHDipole(const YarnDipoleTree& source, int nbodiess,
                            const int blocks, float epssq = 0.0025f,
                            float itolsq = 4.f, float distratiosq = 4.8f,
                            bool use_sorted_order = false);

// Call prepareyarn first to compute the proper w values.
void PrepareYarn(std::vector<float4> &midpts, std::vector<float4> &tangents);

void PrepareYarn(std::vector<float4> &midpts,
                 std::vector<float4> &tangents,
                 std::vector<float4> &dipole1s,
                 std::vector<float4> &dipole3s);

// Call PrepareTree to build the trees. This uses the GPU and
// is not thread-safe; only call this once at a time.
YarnTree PrepareTree(const std::vector<float4> &midpts,
                     const std::vector<float4> &tangents);


// Call this before preparing the dipole trees.
void AllocateYarnRadiiAndBottoms(float** radii, int** bottoms, int* blocks, const int nyarns);

YarnDipoleTree PrepareDipoleTree(const std::vector<float4> &midpts,
                     const std::vector<float4> &tangents,
                     const std::vector<float4> &dipole1s,
                     const std::vector<float4> &dipole3s,
                     float* radius_dev, int* bottom_dev,
                     const int blocks, const int min_threshold);

// call this after PrepareDipoleTree.
void GrabYarnRadiiAndBottoms(std::vector<YarnDipoleTree>& yarn_dipole_trees,
    std::vector<int>& bottoms, float* radii_dev, int* bottoms_dev);

// call this after PrepareDipoleTree and after GrabYarnRadiiAndBottoms
void ComputeKernel3(YarnDipoleTree& yarn_dipole_tree, const int nbodies, const int bottom, const int blocks);


/*
// Example usage:


DEFINE_bool(verify_no_delta, true, "Topology verification: "
    "Disable the delta check (set delta=0) and just use beta.");
DEFINE_double(verify_beta, 4.0, "Topology verification: "
    "Set the Barnes-Hut Beta value.");
DEFINE_double(verify_delta, 2.2, "Topology verification: "
    "Set the FMM correction delta value.");

std::vector<float4*> PrepareGpuPointsAligned(const std::vector<Yarn>& yarns) {
  std::vector<std::vector<float>> yarn_pts(yarns.size());
  // #pragma omp parallel for // this is faster without parallelization
  for(int i = 0; i < yarns.size(); ++i) {
    yarn_pts[i] = std::vector<float>(4 * yarns[i].q.size() + 4);
    for(int j = 0; j < yarns[i].q.size(); ++j) {
      std::memcpy(&(yarn_pts[i][4*j]), &(yarns[i].q[j]), 3 * sizeof(float));
    }
    if (yarns[i].q.size() > 0) std::memcpy(&(yarn_pts[i][4*yarns[i].q.size()]), &(yarns[i].q[0]), 3 * sizeof(float));
  }
  return CopyYarnPointsAligned(yarn_pts);
}

template <typename T>
void FreeGpuPoints(std::vector<T*>& to_free) {
  for (int i = 0; i < to_free.size(); ++i) {
    cudaFree(to_free[i]);
    to_free[i] = nullptr;
  }
  to_free.clear();
}
std::vector<YarnDipoleTree> PrepareGpuBHDipoleTrees(const std::vector<Yarn>& yarns, int& blocks, const int min_threshold, bool optimize_for_two) {
  int nyarns = yarns.size();
  // If there are just two yarns, just make one of the trees.
  std::vector<YarnDipoleTree> yarn_trees(nyarns);
  std::vector<int> bottoms(nyarns);
  float* radii_dev;
  int* bottoms_dev;
  
  AllocateYarnRadiiAndBottoms(&radii_dev, &bottoms_dev, &blocks, nyarns);
  int small = 0;
  int big = 1;
  if (optimize_for_two && (yarns.size() == 2) && (yarns[small].q.size() > yarns[big].q.size())) {
    std::swap(small, big);
  }
  // This prep can now be parallel. It's still faster single-threaded.
  // #pragma omp parallel for
  for(int i = 0; i < nyarns; ++i) {
    int tmp_min_threshold = min_threshold;
    // If there are just two yarns, make a tree for the smaller yarn.
    if (optimize_for_two && (yarns.size() == 2) && (i == big)) {
      tmp_min_threshold = yarns[i].q.size() + 1;
    }  
    std::vector<Point3> midpoints = yarns[i].getMidpoints();
    std::vector<Point3> tangents = yarns[i].getTangents();
    const int n = midpoints.size();
    if (n != tangents.size()) {
      std::cout << "PrepareGpuBHDipoleTrees: Yarn " << i << " midpoints and tangents aren't the same size." << std::endl;
      throw std::runtime_error("PrepareGphBHDipoleTrees midpt tangent size issue.");
    }
    std::vector<float4> midpts_f4(n);
    std::vector<float4> tangents_f4(n);
    std::vector<float4> dipole1s_f4(n);
    std::vector<float4> dipole3s_f4(n);
    for(int j = 0; j < n; ++j) {
      std::memcpy(&(midpts_f4[j]), &(midpoints[j]), 3 * sizeof(float));
      std::memcpy(&(tangents_f4[j]), &(tangents[j]), 3 * sizeof(float));
    }
    PrepareYarn(midpts_f4, tangents_f4, dipole1s_f4, dipole3s_f4);
    yarn_trees[i] = PrepareDipoleTree(midpts_f4, tangents_f4, dipole1s_f4, dipole3s_f4, radii_dev+i, bottoms_dev+i, blocks, tmp_min_threshold);
  }
  GrabYarnRadiiAndBottoms(yarn_trees, bottoms, radii_dev, bottoms_dev);
  // #pragma omp parallel for
  for (int i = 0; i < nyarns; ++i) {
    // If there are just two yarns, make a tree for the smaller yarn.
    const int n = yarns[i].q.size();
    ComputeKernel3(yarn_trees[i], n, bottoms[i], blocks);
  }
  return yarn_trees;
}

void FreeGpuYarnDipoleTrees(std::vector<YarnDipoleTree>& to_free) {
  for (int i = 0; i < to_free.size(); ++i) {
    if (to_free[i].midpts != nullptr) cudaFree(to_free[i].midpts);
    if (to_free[i].tangents != nullptr) cudaFree(to_free[i].tangents);
    if (to_free[i].dipole1s != nullptr) cudaFree(to_free[i].dipole1s);
    if (to_free[i].dipole3s != nullptr) cudaFree(to_free[i].dipole3s);
    if (to_free[i].sortl != nullptr) cudaFree(to_free[i].sortl);
    if (to_free[i].childl != nullptr) cudaFree(to_free[i].childl);
  }
  to_free.clear();
}

std::vector<Yarn> discretized_yarns = Discretize();
// warm up the GPU before running timing experiments.
std::vector<float4*> yarn_gpu_pts = PrepareGpuPointsAligned(discretized_yarns);
FreeGpuPoints(yarn_gpu_pts);

// barnes-hut params
constexpr float epssq = 1e-14f;
const float itolsq = FLAGS_verify_beta * FLAGS_verify_beta;
const float itolsq_dipole = FLAGS_verify_beta * FLAGS_verify_beta; // normal accuracy // beta^2
//constexpr float itolsq_dipole = 55.f; // high accuracy
const float distratiosq = FLAGS_verify_no_delta ? (FLAGS_verify_beta ) *(FLAGS_verify_beta )  : FLAGS_verify_delta * FLAGS_verify_delta; // normal accuracy // delta

begin = std::chrono::steady_clock::now();
int blocks;
std::vector<YarnDipoleTree> yarn_bh_trees = PrepareGpuBHDipoleTrees(discretized_yarns, blocks, 0, true);
end = std::chrono::steady_clock::now();
int64_t prep_time = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
std::cout << "GPU BH Dip Prep time: " << (1e-3f * prep_time) << " ms." << std::endl;

begin = std::chrono::steady_clock::now();
int skip = 0;
int computed = 0;
float error = 0.f;
int err_count = 0;
float max_error = 0.f;
//#pragma omp parallel for
for(int i = 0; i < discretized_yarns.size(); ++i) {
  skip += discretized_yarns.size() - (i+1) - yarn_overlap_unique_per_yarn[i].size();
  computed += yarn_overlap_unique_per_yarn[i].size();
  int64_t nsi = discretized_yarns[i].NumSegments();
  for (int j : yarn_overlap_unique_per_yarn[i]) {
    // std::cout << "Overlap: (" << i << "," << j << "), performing computation." << std::endl;
    float linking_number = std::numeric_limits<float>::quiet_NaN();

    int big = j;
    int small = i;
    int nb = discretized_yarns[j].q.size();
    int ns = discretized_yarns[i].q.size();
    if (nb < ns) {
      std::swap(big, small);
      std::swap(nb, ns);
    }
    float beta2 = itolsq_dipole;
    float delta2 = distratiosq;
    int64_t nsj = discretized_yarns[j].NumSegments();
    if (nsi * nsj >= 100000000ULL) {
      std::cout << "Using beta=16 on (" << i << ", " << j << ")" << std::endl;
      beta2 = 16.0f * 16.0f;
      delta2 = beta2;
    }
    linking_number = ComputeLinkingNumberBHDipole(
        yarn_bh_trees[small], ns,
        yarn_bh_trees[big], nb,
        blocks,
        epssq, beta2,
        delta2);
    float nearest_int = std::roundf(linking_number);
    if (FLAGS_use_input_links) nearest_int = in_matrix.coeff(i,j);
    max_error = std::max(max_error, (float) std::abs(nearest_int - linking_number));
    error += std::abs(nearest_int - linking_number);
    ++err_count;
    if (!FLAGS_squelch_links && !(std::abs(linking_number - in_matrix.coeff(i, j)) < 0.25f)) {  
      //#pragma omp critical
      {
        std::cout << "-------------- (" << i << "," << j << "): computed link: " << linking_number << ", expected: " << in_matrix.coeff(i,j) << std::endl;
      }
    }
    if (!(std::abs(linking_number) < 0.25f)) {
      out_triplet_list.push_back(Eigen::Triplet<int>(i, j, std::roundf(linking_number)));
    }
  }
}
std::cout << "L0 max error: " << max_error << " among " << err_count << " computes." << std::endl;
std::cout << "L1 avg error: " << (error / err_count) << " among " << err_count << " computes. <---- BH GPU ERROR: " << (error/err_count) << std::endl;
std::cout << "Skipped: " << skip << ", Computed: " << computed << std::endl;;
end = std::chrono::steady_clock::now();
time = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
std::cout << "GPU BH Dip Eval time: " << (1e-3f * time) << " ms.             <---- BH GPU TIME: " << (1e-3f * (prep_time + time)) << std::endl;
FreeGpuYarnDipoleTrees(yarn_bh_trees);

*/
