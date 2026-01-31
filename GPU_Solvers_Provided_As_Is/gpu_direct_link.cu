/**
 * Copyright 2020 Ante Qu.  All rights reserved.
 */

#include "gpu_direct_link.h"
#include "helper_math.h"
#include <vector>
#include <thrust/device_vector.h>
#include <thrust/reduce.h>

// In this file, every instance of "Yarn" simply refers to a "Curve". This code is distributed as-is.
// "helper_math.h" is from Nvidia, please find it here for example:
// https://github.com/NVIDIA/cuda-samples/blob/master/Common/helper_math.h

// Some optimizations include tiling (preload "tiles" in parallel in shared memory) and taking advantage of texture memory ("read-only cache").

// Also 128 and 4 works pretty well.
constexpr int blocksize = 64;
constexpr int blockfactor = 8;

//1D grid of 1D blocks
__device__ int getGlobalIdx_1D_1D()
{
	return blockIdx.x *blockDim.x + threadIdx.x;
}

__host__ std::vector<float*> CopyYarnPoints(const std::vector<std::vector<float>>& yarn_points) {
  std::vector<float*> yarn_point_pointers(yarn_points.size());
  for (int i = 0; i < yarn_points.size(); ++i) {
    // loop it around
    cudaMalloc((void**)&(yarn_point_pointers[i]), sizeof(float) * (yarn_points[i].size() + 3));
    cudaMemcpy((void*) yarn_point_pointers[i], (void*) yarn_points[i].data(), sizeof(float) * yarn_points[i].size(), cudaMemcpyHostToDevice);
    cudaMemcpy((void*) (yarn_point_pointers[i] + yarn_points[i].size() ), (void*) yarn_points[i].data(), sizeof(float) * 3, cudaMemcpyHostToDevice);
  }
  return yarn_point_pointers;
}

__host__ std::vector<float4*> CopyYarnPointsAligned(const std::vector<std::vector<float>>& yarn_points) {
  std::vector<float4*> yarn_point_pointers(yarn_points.size());
  for (int i = 0; i < yarn_points.size(); ++i) {
    // round up to the nearest block size for the shared memory implementation.
    int nearest_blocksize = yarn_points[i].size() + (4 * blocksize * blockfactor - 1);
    nearest_blocksize -= (nearest_blocksize % (4 * blocksize * blockfactor));
    // add one to accomodate the shared memory implementation.
    nearest_blocksize += 4;
    // loop it around
    cudaMalloc((void**)&(yarn_point_pointers[i]), sizeof(float) * nearest_blocksize);
    cudaMemcpy((void*) yarn_point_pointers[i], (void*) yarn_points[i].data(), sizeof(float) * yarn_points[i].size(), cudaMemcpyHostToDevice);
    // cudaMemcpy((void*) (yarn_point_pointers[i] + (yarn_points[i].size()/4) ), (void*) yarn_points[i].data(), sizeof(float4), cudaMemcpyHostToDevice);
  }
  return yarn_point_pointers;
}

__host__ std::vector<cudaTextureObject_t> CopyYarnPointsAlignedTextures(const std::vector<std::vector<float>>& yarn_points) {
  std::vector<cudaTextureObject_t> yarn_point_pointers(yarn_points.size());
  for (int i = 0; i < yarn_points.size(); ++i) {
    // round up to the nearest block size for the shared memory implementation.
    int nearest_blocksize = yarn_points[i].size() + 4 + (4 * blocksize * blockfactor - 1);
    nearest_blocksize -= (nearest_blocksize % (4 * blocksize * blockfactor));
    nearest_blocksize += 4;

    float4* device_pointer = nullptr;

    // loop it around
    cudaMalloc((void**)&(device_pointer), sizeof(float) * (nearest_blocksize));
    cudaMemcpy((void*) device_pointer, (void*) yarn_points[i].data(), sizeof(float) * yarn_points[i].size(), cudaMemcpyHostToDevice);
    cudaMemcpy((void*) (device_pointer + (yarn_points[i].size()/4) ), (void*) yarn_points[i].data(), sizeof(float4), cudaMemcpyHostToDevice);

    // specify texture
    struct cudaResourceDesc resDesc;
    resDesc.resType = cudaResourceTypeLinear;
    resDesc.res.linear.devPtr = device_pointer;
    resDesc.res.linear.desc = cudaCreateChannelDesc<float4>();
    resDesc.res.linear.sizeInBytes = sizeof(float) * nearest_blocksize;

    // specify texture object parameters
    struct cudaTextureDesc texDesc;
    std::memset(&texDesc, 0, sizeof(struct cudaTextureDesc));
    texDesc.sRGB = 0;

    
    // create texture object: we only have to do this once! Blah blah blah this is still a pain.
    cudaCreateTextureObject(&(yarn_point_pointers[i]), &resDesc, &texDesc, NULL);
  }
  return yarn_point_pointers;
}


#define TWO_PI (6.28318530718f)
__global__ void computeLinksAgainstYarn2(float* linking_number, const float* __restrict__ yarn1_points, const float* __restrict__ yarn2_points, int n1, int n2){
  int y1ind = getGlobalIdx_1D_1D(); //yarn 1 index
  if (y1ind >= n1) return;
  float3 y1 = *(reinterpret_cast<const float3*>(& (yarn1_points[3 * y1ind])));
  //int next_ind = 3 * ((y1ind + 1) % n1);
  float3 y1p1 = *(reinterpret_cast<const float3*>(& (yarn1_points[3 * y1ind + 3])));
  float3 y2, y2p1, alpha, beta, gamma, delta;
  float determinant;
  float div1, div2;
  float al, bl, cl, dl, ac;
  float link_sum = 0;
  y2 = *(reinterpret_cast<const float3*>(yarn2_points));
  for (int j = 1; j <= n2; ++j) {
    // next_ind = 3 * ((j + 1) % n2);
    y2p1 = *(reinterpret_cast<const float3*>(& (yarn2_points[3 * j])));
    alpha = y2 - y1;
    beta = y2 - y1p1;
    gamma = y2p1 - y1p1;
    delta = y2p1 - y1;
    determinant = dot(alpha, cross(beta, gamma));
    al = length(alpha);
    bl = length(beta);
    cl = length(gamma);
    dl = length(delta);
    ac = dot(alpha, gamma);
    div1 = al*bl*cl + dot(alpha, beta) * cl + dot(beta, gamma) * al + ac * bl;
    div2 = al*dl*cl + dot(alpha, delta) * cl + dot(delta, gamma) * al + ac * dl;
    link_sum += atan2f(determinant, div1) + atan2f(determinant, div2);
    y2 = y2p1;
  }
  linking_number[y1ind] = link_sum / TWO_PI;
}

__global__ void computeLinksAgainstYarnBoth(float* linking_number, const float* __restrict__ yarn1_points, const float* __restrict__ yarn2_points, int n1, int n2){
  int ind = getGlobalIdx_1D_1D();
  if (ind >= n1 * n2) return;
  int y1ind = ind / n2; //yarn 1 index
  int j = (ind % n2);
  float3 y1 = *(reinterpret_cast<const float3*>(& (yarn1_points[3 * y1ind])));
  //int next_ind = 3 * ((y1ind + 1) % n1);
  float3 y1p1 = *(reinterpret_cast<const float3*>(& (yarn1_points[3 * y1ind + 3])));
  float3 y2, y2p1, alpha, beta, gamma, delta;
  float determinant;
  float div1, div2;
  float al, bl, cl, dl, ac;
  y2 = *(reinterpret_cast<const float3*>(& (yarn2_points[3 * j])));
  //next_ind = 3 * ((j + 1) % n2);
  y2p1 = *(reinterpret_cast<const float3*>(& (yarn2_points[3 * j + 3])));
  alpha = y2 - y1;
  beta = y2 - y1p1;
  gamma = y2p1 - y1p1;
  delta = y2p1 - y1;
  determinant = dot(alpha, cross(beta, gamma));
  al = length(alpha);
  bl = length(beta);
  cl = length(gamma);
  dl = length(delta);
  ac = dot(alpha, gamma);
  div1 = al*bl*cl + dot(alpha, beta) * cl + dot(beta, gamma) * al + ac * bl;
  div2 = al*dl*cl + dot(alpha, delta) * cl + dot(delta, gamma) * al + ac * dl;
  linking_number[ind] = atan2f(determinant, div1) + atan2f(determinant, div2);

}

__global__ void computeLinksAgainstYarnBothFast(float* linking_number, const float* __restrict__ yarn1_points, const float* __restrict__ yarn2_points, int n1, int n2){
  int ind = getGlobalIdx_1D_1D();
  if (ind >= n1 * n2) return;
  int y1ind = ind / n2; //yarn 1 index
  int j = ind % n2;
  float3 y1 = *(reinterpret_cast<const float3*>(& (yarn1_points[3 * y1ind])));
  //int next_ind = 3 * ((y1ind + 1) % n1);
  float3 y1p1 = *(reinterpret_cast<const float3*>(& (yarn1_points[3 * y1ind + 3])));
  float3 y2, y2p1, alpha, beta, gamma, delta;
  float determinant;
  float div1, div2;
  float al, bl, cl, dl, ac;
  y2 = *(reinterpret_cast<const float3*>(& (yarn2_points[3 * j])));
  //next_ind = 3 * ((j + 1) % n2);
  y2p1 = *(reinterpret_cast<const float3*>(& (yarn2_points[3 * j + 3])));
  alpha = y2 - y1;
  beta = y2 - y1p1;
  gamma = y2p1 - y1p1;
  delta = y2p1 - y1;
  determinant = dot(alpha, cross(beta, gamma));
  al = length(alpha);
  bl = length(beta);
  cl = length(gamma);
  dl = length(delta);
  ac = dot(alpha, gamma);
  div1 = al*bl*cl + dot(alpha, beta) * cl + dot(beta, gamma) * al + ac * bl;
  div2 = al*dl*cl + dot(alpha, delta) * cl + dot(delta, gamma) * al + ac * dl;
  // x'
  al = div1 * div2 - determinant * determinant;
  
  // ac is s s'' = y' * determinant. We also want it negative when det = 0 and div1, div2 are both negative.
  ac = (determinant == 0) ? copysignf(1.f, div1) + copysignf(1.f, div2) : (div1 + div2); // = y' * distth
  // y'
  bl = determinant * ac;
  cl = (determinant == 0) ? -copysignf(TWO_PI, div1) : copysignf(TWO_PI, determinant);
  
  linking_number[ind] = atan2f(bl, al) + (((ac < 0) || ((ac == 0) && (determinant < 0))) ? cl : 0);

}


__global__ void computeLinksAgainstYarnBothFastAligned(float* linking_number, const float4* __restrict__ yarn1_points, const float4* __restrict__ yarn2_points, int n1, int n2){
  int ind = getGlobalIdx_1D_1D();
  if (ind >= n1 * n2) return;
  int y1ind = ind / n2; //yarn 1 index
  int j = ind % n2;
  float3 y1 = make_float3(yarn1_points[y1ind]);
  //int next_ind = (y1ind + 1) % n1;
  float3 y1p1 = make_float3(yarn1_points[y1ind+1]);
  float3 y2, y2p1, alpha, beta, gamma, delta;
  float determinant;
  float div1, div2;
  float al, bl, cl, dl, ac;
  y2 = make_float3(yarn2_points[j]);
  //next_ind =  ((j + 1) % n2);
  y2p1 = make_float3(yarn2_points[j+1]);
  alpha = y2 - y1;
  beta = y2 - y1p1;
  gamma = y2p1 - y1p1;
  delta = y2p1 - y1;
  determinant = dot(alpha, cross(beta, gamma));
  al = length(alpha);
  bl = length(beta);
  cl = length(gamma);
  dl = length(delta);
  ac = dot(alpha, gamma);
  div1 = al*bl*cl + dot(alpha, beta) * cl + dot(beta, gamma) * al + ac * bl;
  div2 = al*dl*cl + dot(alpha, delta) * cl + dot(delta, gamma) * al + ac * dl;
  // x'
  al = div1 * div2 - determinant * determinant;
  // ac is s s'' = y' * determinant. We also want it negative when det = 0 and div1, div2 are both negative.
  ac = (determinant == 0) ? copysignf(1.f, div1) + copysignf(1.f, div2) : (div1 + div2); // = y' * distth
  // y'
  bl = determinant * ac;
  cl = (determinant == 0) ? -copysignf(TWO_PI, div1) : copysignf(TWO_PI, determinant);

  linking_number[ind] = atan2f(bl, al) + (((ac < 0) || ((ac == 0) && (determinant < 0))) ? cl : 0);

}


__global__ void computeLinksAgainstYarnBothFastAlignedTexture(float* linking_number, cudaTextureObject_t yarn1_points, cudaTextureObject_t yarn2_points, int n1, int n2){
  int ind = getGlobalIdx_1D_1D();
  if (ind >= n1 * n2) return;
  int y1ind = ind / n2; //yarn 1 index
  int j = ind % n2;
  float3 y1 = make_float3(tex1Dfetch<float4>(yarn1_points, y1ind));
  //int next_ind = (y1ind + 1) % n1;
  float3 y1p1 = make_float3(tex1Dfetch<float4>(yarn1_points,y1ind+1));
  float3 y2, y2p1, alpha, beta, gamma, delta;
  float determinant;
  float div1, div2;
  float al, bl, cl, dl, ac;
  y2 = make_float3(tex1Dfetch<float4>(yarn2_points, j));
  //next_ind =  ((j + 1) % n2);
  y2p1 = make_float3(tex1Dfetch<float4>(yarn2_points, j+1));
  alpha = y2 - y1;
  beta = y2 - y1p1;
  gamma = y2p1 - y1p1;
  delta = y2p1 - y1;
  determinant = dot(alpha, cross(beta, gamma));
  al = length(alpha);
  bl = length(beta);
  cl = length(gamma);
  dl = length(delta);
  ac = dot(alpha, gamma);
  div1 = al*bl*cl + dot(alpha, beta) * cl + dot(beta, gamma) * al + ac * bl;
  div2 = al*dl*cl + dot(alpha, delta) * cl + dot(delta, gamma) * al + ac * dl;
  // x'
  al = div1 * div2 - determinant * determinant;
  // ac is s s'' = y' * determinant. We also want it negative when det = 0 and div1, div2 are both negative.
  ac = (determinant == 0) ? copysignf(1.f, div1) + copysignf(1.f, div2) : (div1 + div2); // = y' * distth
  // y'
  bl = determinant * ac;
  cl = (determinant == 0) ? -copysignf(TWO_PI, div1) : copysignf(TWO_PI, determinant);

  linking_number[ind] = atan2f(bl, al) + (((ac < 0) || ((ac == 0) && (determinant < 0))) ? cl : 0);

}

__global__ void computeLinksAgainstYarn2Fast(float* linking_number, const float* __restrict__ yarn1_points, const float* __restrict__ yarn2_points, int n1, int n2){
  int y1ind = getGlobalIdx_1D_1D(); //yarn 1 index
  if (y1ind >= n1) return;
  float3 y1 = *(reinterpret_cast<const float3*>(& (yarn1_points[3 * y1ind])));
  //int next_ind = 3 * ((y1ind + 1) % n1);
  float3 y1p1 = *(reinterpret_cast<const float3*>(& (yarn1_points[3 * y1ind + 3])));
  float3 y2, y2p1, alpha, beta, gamma, delta;
  float determinant;
  float div1, div2;
  float al, bl, cl, dl, ac;
  float S = -1.f;
  float xs = 1.f;
  float ys = 0;
  float offsets = 0;
  y2 = *(reinterpret_cast<const float3*>(yarn2_points));
  for (int j = 1; j <= n2; ++j) {
    //next_ind = 3 * ((j + 1) % n2);
    y2p1 = *(reinterpret_cast<const float3*>(& (yarn2_points[3 * j])));
    alpha = y2 - y1;
    beta = y2 - y1p1;
    gamma = y2p1 - y1p1;
    delta = y2p1 - y1;
    determinant = dot(alpha, cross(beta, gamma));
    al = length(alpha);
    bl = length(beta);
    cl = length(gamma);
    dl = length(delta);
    ac = dot(alpha, gamma);
    div1 = al*bl*cl + dot(alpha, beta) * cl + dot(beta, gamma) * al + ac * bl;
    div2 = al*dl*cl + dot(alpha, delta) * cl + dot(delta, gamma) * al + ac * dl;
    // xprime for the two
    al = div1 * div2 - determinant * determinant;
    // yprime for the two
    bl = determinant * (div1 + div2);
    // ac is s s'' = y' * determinant. We also want it negative when det = 0 and div1, div2 are both negative.
    ac = (determinant == 0) ? copysignf(1.f, div1) + copysignf(1.f, div2) : (div1 + div2);
    cl = (determinant == 0) ? -copysignf(1.f, div1) : copysignf(1.f, determinant);
    // add the offset for the two angles together
    offsets += (((ac < 0) || ((ac == 0) && (determinant < 0))) ? cl : 0);

    // new y" = x1 y2 + x2 y1
    div2 = xs * bl + ys * al;
    // new x" = x1 x2 - y1 y2
    div1 = xs * al - ys * bl;
    // use max(|x"|, |y"|) to scale
    cl = fmaxf(fabsf(div1), fabsf(div2));
    div1 /= cl;
    div2 /= cl;
    // check if signs are equal: if new y is a different sign than old y, add sign(old y) to offset
    // al = bl * determinant;
    
    // this is sign(x",y")
    cl = ((div2 > 0.f) || (div2 == 0.f && div1 < 0.f)) ? 1.f : -1.f;
    
    // sy sy'
    bl = (((bl > 0) || (bl == 0 && al < 0)) ? 1.f : -1.f ) * S;
    // Cross detection expression: (sy sy" < 0) && (sy sy' >= 0)
    al = ((cl * S < 0) && (bl >= 0)) ? 1.f : 0.f;
    offsets += S * al;
    xs = div1;
    ys = div2;
    S = cl;
    y2 = y2p1;
  }
  linking_number[y1ind] = atan2(ys, xs) / TWO_PI + offsets;
}

__global__ void computeLinksAgainstYarn2FastAligned(float* linking_number, const float4* __restrict__ yarn1_points, const float4* __restrict__ yarn2_points, int n1, int n2){
  int y1ind = getGlobalIdx_1D_1D(); //yarn 1 index
  if (y1ind >= n1) return;
  float3 y1 = make_float3(yarn1_points[y1ind]);
  //int next_ind = 3 * ((y1ind + 1) % n1);
  float3 y1p1 = make_float3(yarn1_points[y1ind + 1]);
  float3 y2, y2p1, alpha, beta, gamma, delta;
  float determinant;
  float div1, div2;
  float al, bl, cl, dl, ac;
  float S = -1.f;
  float xs = 1.f;
  float ys = 0;
  float offsets = 0;
  y2 = make_float3(yarn2_points[0]);
  for (int j = 1; j <= n2; ++j) {
    //next_ind = 3 * ((j + 1) % n2);
    y2p1 = make_float3(yarn2_points[j]);
    alpha = y2 - y1;
    beta = y2 - y1p1;
    gamma = y2p1 - y1p1;
    delta = y2p1 - y1;
    determinant = dot(alpha, cross(beta, gamma));
    al = length(alpha);
    bl = length(beta);
    cl = length(gamma);
    dl = length(delta);
    ac = dot(alpha, gamma);
    div1 = al*bl*cl + dot(alpha, beta) * cl + dot(beta, gamma) * al + ac * bl;
    div2 = al*dl*cl + dot(alpha, delta) * cl + dot(delta, gamma) * al + ac * dl;
    // xprime for the two
    al = div1 * div2 - determinant * determinant;
    // yprime for the two
    bl = determinant * (div1 + div2);
    // ac is s s'' = y' * determinant. We also want it negative when det = 0 and div1, div2 are both negative.
    ac = (determinant == 0) ? copysignf(1.f, div1) + copysignf(1.f, div2) : (div1 + div2);
    cl = (determinant == 0) ? -copysignf(1.f, div1) : copysignf(1.f, determinant);
    // add the offset for the two angles together
    offsets += (((ac < 0) || ((ac == 0) && (determinant < 0))) ? cl : 0);

    // new y" = x1 y2 + x2 y1
    div2 = xs * bl + ys * al;
    // new x" = x1 x2 - y1 y2
    div1 = xs * al - ys * bl;
    // use max(|x"|, |y"|) to scale
    cl = fmaxf(fabsf(div1), fabsf(div2));
    div1 /= cl;
    div2 /= cl;
    // check if signs are equal: if new y is a different sign than old y, add sign(old y) to offset
    // al = bl * determinant;
    
    // this is sign(x",y")
    cl = ((div2 > 0.f) || (div2 == 0.f && div1 < 0.f)) ? 1.f : -1.f;
    
    // sy sy'
    bl = (((bl > 0) || (bl == 0 && al < 0)) ? 1.f : -1.f ) * S;
    // Cross detection expression: (sy sy" < 0) && (sy sy' >= 0)
    al = ((cl * S < 0) && (bl >= 0)) ? 1.f : 0.f;
    offsets += S * al;
    xs = div1;
    ys = div2;
    S = cl;
    y2 = y2p1;
  }
  linking_number[y1ind] = atan2(ys, xs) / TWO_PI + offsets;
}

// this implementation uses shared memory to preload the points quickly (GPU tiling).
__global__ void computeLinksAgainstYarn2FastAlignedSm(float* linking_number, const float4* __restrict__ yarn1_points, const float4* __restrict__ yarn2_points, int n1, int n2){
  const int y1ind = getGlobalIdx_1D_1D(); //yarn 1 index
  // if (y1ind >= n1) y1ind = n1 - 1;
  __shared__ float4 shmem[blocksize * blockfactor];
  const int tid = threadIdx.x;
  float3 y1 = make_float3(yarn1_points[y1ind]);
  //int next_ind = 3 * ((y1ind + 1) % n1);
  float3 y1p1 = make_float3(yarn1_points[y1ind + 1]);
  float3 y2, y2p1, alpha, beta, gamma, delta;
  float determinant;
  float div1, div2;
  float al, bl, cl, dl, ac;
  float S = -1.f;
  float xs = 1.f;
  float ys = 0;
  float offsets = 0;
  y2 = make_float3(yarn2_points[0]);
  for (int j = 1; j <= n2; j+= blocksize * blockfactor) {
    for (int k = 0; k < blocksize * blockfactor; k += blocksize) {
      if (tid + k <= n2) shmem[tid + k] = yarn2_points[j + tid + k];
    }
    __syncthreads();
    for (int k = 0; k < blocksize * blockfactor; ++k) {
      //next_ind = 3 * ((j + 1) % n2);
      if (j + k > n2) break;
      y2p1 = make_float3(shmem[k]);
      alpha = y2 - y1;
      beta = y2 - y1p1;
      gamma = y2p1 - y1p1;
      delta = y2p1 - y1;
      determinant = dot(alpha, cross(beta, gamma));
      al = length(alpha);
      bl = length(beta);
      cl = length(gamma);
      dl = length(delta);
      ac = dot(alpha, gamma);
      div1 = al*bl*cl + dot(alpha, beta) * cl + dot(beta, gamma) * al + ac * bl;
      div2 = al*dl*cl + dot(alpha, delta) * cl + dot(delta, gamma) * al + ac * dl;
      // xprime for the two
      al = div1 * div2 - determinant * determinant;
      // yprime for the two
      bl = determinant * (div1 + div2);
      // ac is s s'' = y' * determinant. We also want it negative when det = 0 and div1, div2 are both negative.
      ac = (determinant == 0) ? copysignf(1.f, div1) + copysignf(1.f, div2) : (div1 + div2);
      cl = (determinant == 0) ? -copysignf(1.f, div1) : copysignf(1.f, determinant);
      // add the offset for the two angles together
      offsets += (((ac < 0) || ((ac == 0) && (determinant < 0))) ? cl : 0);

      // new y" = x1 y2 + x2 y1
      div2 = xs * bl + ys * al;
      // new x" = x1 x2 - y1 y2
      div1 = xs * al - ys * bl;
      // use max(|x"|, |y"|) to scale
      cl = fmaxf(fabsf(div1), fabsf(div2));
      div1 /= cl;
      div2 /= cl;
      // check if signs are equal: if new y is a different sign than old y, add sign(old y) to offset
      // al = bl * determinant;
      
      // this is sign(x",y")
      cl = ((div2 > 0.f) || (div2 == 0.f && div1 < 0.f)) ? 1.f : -1.f;
      
      // sy sy'
      bl = (((bl > 0) || (bl == 0 && al < 0)) ? 1.f : -1.f ) * S;
      // Cross detection expression: (sy sy" < 0) && (sy sy' >= 0)
      al = ((cl * S < 0) && (bl >= 0)) ? 1.f : 0.f;
      offsets += S * al;
      xs = div1;
      ys = div2;
      S = cl;
      y2 = y2p1;
    }
    __syncthreads();
  }
  if (y1ind >= n1) return;
  linking_number[y1ind] = atan2(ys, xs) / TWO_PI + offsets;
}

// these pointers are device pointers
__host__ float ComputeLink(const float* __restrict__ yarn1_points, const float* __restrict__ yarn2_points, int n1, int n2) {
    thrust::device_vector<float> l_ns(n1);
    //KernelConfig kc(n1);
    int gridsize = (n1 + blocksize - 1) / blocksize;
    computeLinksAgainstYarn2Fast<<<gridsize, blocksize>>>(l_ns.data().get(), yarn1_points, yarn2_points, n1, n2);
    cudaDeviceSynchronize();
    return thrust::reduce(l_ns.begin(), l_ns.end());
}

// these pointers are device pointers
__host__ float ComputeLinkBoth(const float* __restrict__ yarn1_points, const float* __restrict__ yarn2_points, int n1, int n2) {
    thrust::device_vector<float> l_ns(n1*n2);
    //KernelConfig kc(n1);
    int gridsize = (n1*n2 + blocksize - 1) / blocksize;
    computeLinksAgainstYarnBothFast<<<gridsize, blocksize>>>(l_ns.data().get(), yarn1_points, yarn2_points, n1, n2);
    cudaDeviceSynchronize();
    return thrust::reduce(l_ns.begin(), l_ns.end()) / TWO_PI;
}

// these pointers are device pointers
__host__ float ComputeLinkAligned(const float4* __restrict__ yarn1_points, const float4* __restrict__ yarn2_points, int n1, int n2) {
    thrust::device_vector<float> l_ns(n1);
    //KernelConfig kc(n1);
    int gridsize = (n1 + blocksize - 1) / blocksize;
    computeLinksAgainstYarn2FastAligned<<<gridsize, blocksize>>>(l_ns.data().get(), yarn1_points, yarn2_points, n1, n2);
    cudaDeviceSynchronize();
    return thrust::reduce(l_ns.begin(), l_ns.end());
}

// these pointers are device pointers
__host__ float ComputeLinkAlignedSm(const float4* __restrict__ yarn1_points, const float4* __restrict__ yarn2_points, int n1, int n2) {
    thrust::device_vector<float> l_ns(n1);
    //KernelConfig kc(n1);
    int gridsize = (n1 + blocksize - 1) / blocksize;
    computeLinksAgainstYarn2FastAlignedSm<<<gridsize, blocksize, blocksize * sizeof(float4)>>>(l_ns.data().get(), yarn1_points, yarn2_points, n1, n2);
    cudaDeviceSynchronize();
    return thrust::reduce(l_ns.begin(), l_ns.end());
}


// these pointers are device pointers
__host__ float ComputeLinkBothAligned(const float4* __restrict__ yarn1_points, const float4* __restrict__ yarn2_points, int n1, int n2) {
    thrust::device_vector<float> l_ns(n1*n2);
    //KernelConfig kc(n1);
    int gridsize = (n1*n2 + blocksize - 1) / blocksize;
    computeLinksAgainstYarnBothFastAligned<<<gridsize, blocksize>>>(l_ns.data().get(), yarn1_points, yarn2_points, n1, n2);
    cudaDeviceSynchronize();
    return thrust::reduce(l_ns.begin(), l_ns.end()) / TWO_PI;
}

// these pointers are device pointers
__host__ float ComputeLinkBothAlignedTexture(cudaTextureObject_t yarn1_points, cudaTextureObject_t yarn2_points, int n1, int n2) {
    thrust::device_vector<float> l_ns(n1*n2);
    //KernelConfig kc(n1);
    int gridsize = (n1*n2 + blocksize - 1) / blocksize;
    computeLinksAgainstYarnBothFastAlignedTexture<<<gridsize, blocksize>>>(l_ns.data().get(), yarn1_points, yarn2_points, n1, n2);
    cudaDeviceSynchronize();
    return thrust::reduce(l_ns.begin(), l_ns.end()) / TWO_PI;
}
