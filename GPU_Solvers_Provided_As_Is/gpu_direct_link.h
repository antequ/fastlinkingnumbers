/**
 * Copyright 2020 Ante Qu.  All rights reserved.
 */

// In this file, every instance of "Yarn" simply refers to a "Curve". This code is distributed as-is.
// "helper_math.h" is from Nvidia, please find it here for example:
// https://github.com/NVIDIA/cuda-samples/blob/master/Common/helper_math.h

#pragma once
#include <vector>
#include "helper_math.h"

// Loads memory from the host onto the GPU.
std::vector<float*> CopyYarnPoints(const std::vector<std::vector<float>>& yarn_points);

// All functions below take device pointers as inputs.

// These are unoptimized versions provided for clarity and readability.
float ComputeLink(const float* yarn1_points, const float* yarn2_points, int n1, int n2);
float ComputeLinkBoth(const float* yarn1_points, const float* yarn2_points, int n1, int n2);
std::vector<float4*> CopyYarnPointsAligned(const std::vector<std::vector<float>>& yarn_points);
float ComputeLinkAligned(const float4* yarn1_points, const float4* yarn2_points, int n1, int n2);

// Use this optimized version for small examples -- it can saturate the GPU and give much faster results.
// This is just parallelizing across both for loops.
float ComputeLinkBothAligned(const float4* yarn1_points, const float4* yarn2_points, int n1, int n2);

// Use this optimized version for large examples -- e.g. >20,000 ctrlpts per yarn -- this uses less memory.
// This parallelizes across one for loop and takes advantage of tiling and read-only cache.
float ComputeLinkAlignedSm(const float4* yarn1_points, const float4* yarn2_points, int n1, int n2);

// this performs the same as the above read-only cache implementations.
typedef unsigned long long cudaTextureObject_t;
std::vector<cudaTextureObject_t> CopyYarnPointsAlignedTextures(const std::vector<std::vector<float>>& yarn_points);
float ComputeLinkBothAlignedTexture(cudaTextureObject_t yarn1_points, cudaTextureObject_t yarn2_points, int n1, int n2);


/* 

// Example usage:
    
constexpr uint64_t kDeThreshold = 1400000ULL;

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

std::vector<Yarn> discretized_yarns = Discretize();

// Warm up the GPU before timing it.
std::vector<float4*> yarn_gpu_pts = PrepareGpuPointsAligned(discretized_yarns);
FreeGpuPoints(yarn_gpu_pts);
begin = std::chrono::steady_clock::now();
yarn_gpu_pts = PrepareGpuPointsAligned(discretized_yarns);
end = std::chrono::steady_clock::now();
int64_t prep_time = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
std::cout << "GPU DE Preparation time: " << (1e-3f * prep_time) << " ms." << std::endl;
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
    for (int j : yarn_overlap_unique_per_yarn[i]) {
    // std::cout << "Overlap: (" << i << "," << j << "), performing computation." << std::endl;
    float linking_number = std::numeric_limits<float>::quiet_NaN();
    if (static_cast<uint64_t>(discretized_yarns[i].q.size()) *
        static_cast<uint64_t>(discretized_yarns[j].q.size()) >
        kDeThreshold) {
        // std::cout << "Using 1D aligned shared memory implementation instead of 2D for " << i << ", " << j << std::endl;
        int big = i;
        int small = j;
        int nb = discretized_yarns[i].q.size();
        int ns = discretized_yarns[j].q.size();
        if (nb < ns) {
        std::swap(big, small);
        std::swap(nb, ns);
        }
        linking_number = ComputeLinkAlignedSm(
            yarn_gpu_pts[big], yarn_gpu_pts[small], nb, ns);
    } else {
        linking_number = ComputeLinkBothAligned(
            yarn_gpu_pts[i], yarn_gpu_pts[j],
            discretized_yarns[i].q.size(),
            discretized_yarns[j].q.size());
    }
    float nearest_int = std::roundf(linking_number);
    
    if (FLAGS_use_input_links) nearest_int = in_matrix.coeff(i,j);
    max_error = std::max(max_error, (float) std::abs(nearest_int - linking_number));
    error += std::abs(nearest_int - linking_number);
    ++err_count;
    if (!FLAGS_squelch_links && !(std::abs(linking_number - in_matrix.coeff(i, j)) < 0.25f)) {  
        //#pragma omp critical
        {
        std::cout << "-------------- (" << i << "," << j << "): link = " << linking_number << "." << std::endl;
        }
    }
    if (!use_gpu_bh_dipole_link && !(std::abs(linking_number) < 0.25f)) {
        out_triplet_list.push_back(Eigen::Triplet<int>(i, j, std::roundf(linking_number)));
    }
    }
}
std::cout << "L0 max error: " << max_error << " among " << err_count << " computes." << std::endl;
std::cout << "L1 avg error: " << (error / err_count) << " among " << err_count << " computes." << std::endl;
std::cout << "Skipped: " << skip << ", Computed: " << computed << std::endl;;
end = std::chrono::steady_clock::now();
time = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
std::cout << "GPU DE Evaluation time: " << (1e-3f * time) << " ms.          <---- DE GPU TIME: " << (1e-3f * (time + prep_time)) << std::endl;
FreeGpuPoints(yarn_gpu_pts);

*/