// Ante Qu, 2020.
// This Barnes-Hut linking integral evaluator has been significantly modified from
// the Texas State University's ECL-BH code. I've included their disclaimer
// here:

/*
ECL-BH v4.5: Simulation of the gravitational forces in a star cluster using
the Barnes-Hut n-body algorithm.

Copyright (c) 2010-2020 Texas State University. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

   * Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.
   * Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.
   * Neither the name of Texas State University nor the names of its
     contributors may be used to endorse or promote products derived from
     this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL TEXAS STATE UNIVERSITY BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Authors: Martin Burtscher and Sahar Azimi

URL: The latest version of this code is available at
https://userweb.cs.txstate.edu/~burtscher/research/ECL-BH/.

Publication: This work is described in detail in the following paper.
Martin Burtscher and Keshav Pingali. An Efficient CUDA Implementation of the
Tree-based Barnes Hut n-Body Algorithm. Chapter 6 in GPU Computing Gems
Emerald Edition, pp. 75-92. January 2011.
*/

// Compile:
// nvcc -O1 -arch=sm_35 -ccbin=g++ --compiler-options "-fopenmp -funroll-loops
// -O3" -use_fast_math bh_link_main.cu -o bh_link -I../common

#include <cuda.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

// for linking number input/output processing, float3 dot/cross products, etc
#include "helper_math.h"
#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <vector>
#include "bh_link.h"
#include <math_constants.h>

// tree output structures: child, sort, posMass (-> posTanz), and Tanxy.
// posTanz.w will be the z portion of tangent.
// TODO(antequ 12012020):
// Just convert it to pos and mass float4s instead of the combined.
// we can have an acc scalar float too.

// threads per block
#define THREADS1 1024 /* must be a power of 2 */
#define THREADS2 1024
// make this 288 for 20*288*8<49152 - only 10% slower than 768.
#define THREADS3 768 /* shared-memory limited on some devices */
#define THREADS4 1024
#define THREADS5 1024
#define THREADS6 1024

// block count = factor * #SMs
#define FACTOR1 2
#define FACTOR2 2
// experimentally this works better than FACTOR3=1.
#define FACTOR3 1 /* must all be resident at the same time */
#define FACTOR4 1 /* must all be resident at the same time */
#define FACTOR5 2
#define FACTOR6 2
// values from 64 to 128 were best
#define OVERSATURATE_FACTOR5 16
#define WARPSIZE 32
#define MAXDEPTH 32
// #define DISTRATIO2 (4.8f)
#define DISTRATIO2 (4.8f)
#define INV_FOUR_PI (0.079577471545947667884f)
// ONLY RUN ONE KERNEL AT A TIME.
// set by BoundingBoxKernel, TreeBuildingKernel, and ClearKernel2
//__device__ volatile int bottomd;
__device__ unsigned int blkcntd;
//__device__ volatile float radiusd;

struct TreeMeta {
  int bottom;
  float radius;
};

// set by BoundingBoxKernel, TreeBuildingKernel, and ClearKernel2
// __device__ volatile TreeMeta tree_meta_d;

/******************************************************************************/
/*** initialize memory ********************************************************/
/******************************************************************************/

__global__ void InitializationKernelDipole() { blkcntd = 0; }

/******************************************************************************/
/*** compute center and radius ************************************************/
/******************************************************************************/

__global__ __launch_bounds__(THREADS1, FACTOR1) void BoundingBoxKernelDipole(
    const int nnodesd, const int nbodiesd, int *const __restrict__ startd,
    int *const __restrict__ childd, float4 *const __restrict__ midpointsd,
    float3 *const __restrict__ maxd, float3 *const __restrict__ mind,
    int* bottom_dev, float* radius_dev) {
  int i, j, k, inc;
  float val;
  __shared__ volatile float sminx[THREADS1], smaxx[THREADS1], sminy[THREADS1],
      smaxy[THREADS1], sminz[THREADS1], smaxz[THREADS1];
  float3 min, max;

  // initialize with valid data (in case #bodies < #threads)
  const float4 p0 = midpointsd[0];
  min.x = max.x = p0.x;
  min.y = max.y = p0.y;
  min.z = max.z = p0.z;

  // scan all bodies
  i = threadIdx.x;
  inc = THREADS1 * gridDim.x;
  for (j = i + blockIdx.x * THREADS1; j < nbodiesd; j += inc) {
    const float4 p = midpointsd[j];
    val = p.x;
    min.x = fminf(min.x, val);
    max.x = fmaxf(max.x, val);
    val = p.y;
    min.y = fminf(min.y, val);
    max.y = fmaxf(max.y, val);
    val = p.z;
    min.z = fminf(min.z, val);
    max.z = fmaxf(max.z, val);
  }

  // reduction in shared memory
  sminx[i] = min.x;
  smaxx[i] = max.x;
  sminy[i] = min.y;
  smaxy[i] = max.y;
  sminz[i] = min.z;
  smaxz[i] = max.z;

  for (j = THREADS1 / 2; j > 0; j /= 2) {
    __syncthreads();
    if (i < j) {
      k = i + j;
      sminx[i] = min.x = fminf(min.x, sminx[k]);
      smaxx[i] = max.x = fmaxf(max.x, smaxx[k]);
      sminy[i] = min.y = fminf(min.y, sminy[k]);
      smaxy[i] = max.y = fmaxf(max.y, smaxy[k]);
      sminz[i] = min.z = fminf(min.z, sminz[k]);
      smaxz[i] = max.z = fmaxf(max.z, smaxz[k]);
    }
  }

  // write block result to global memory
  if (i == 0) {
    k = blockIdx.x;
    mind[k] = min;
    maxd[k] = max;
    __threadfence();

    inc = gridDim.x - 1;
    if (inc == atomicInc(&blkcntd, inc)) {
      // I'm the last block, so combine all block results
      for (j = 0; j <= inc; j++) {
        float3 minp = mind[j];
        float3 maxp = maxd[j];
        min.x = fminf(min.x, minp.x);
        max.x = fmaxf(max.x, maxp.x);
        min.y = fminf(min.y, minp.y);
        max.y = fmaxf(max.y, maxp.y);
        min.z = fminf(min.z, minp.z);
        max.z = fmaxf(max.z, maxp.z);
      }

      // compute radius
      val = fmaxf(max.x - min.x, max.y - min.y);
      *radius_dev = fmaxf(val, max.z - min.z) * 0.5f;

      // ADDED TO FIX tree radius to be constant across yarns. 
      // tree_meta_d.radius = 54.f; // see Alien sweater dims.

      // create root node
      k = nnodesd;
      *bottom_dev = k;

      startd[k] = 0;
      float4 p;
      p.x = (min.x + max.x) * 0.5f;
      p.y = (min.y + max.y) * 0.5f;
      p.z = (min.z + max.z) * 0.5f;
      // ADDED TO FIX tree center to be constant across yarns. 
      // p.x = 0.f;
      // p.y = 0.f;
      // p.z = 0.f;
      p.w = -1.0f;
      midpointsd[k] = p;
      k *= 8;
      for (i = 0; i < 8; i++)
        childd[k + i] = -1;

      // stepd++;
    }
  }
}

/******************************************************************************/
/*** build tree ***************************************************************/
/******************************************************************************/

__global__ __launch_bounds__(1024, 1) void ClearKernelDipole1(
    const int nnodesd, const int nbodiesd, int *const __restrict__ childd) {
  int k, inc, top, bottom;

  top = 8 * nnodesd;
  bottom = 8 * nbodiesd;
  inc = blockDim.x * gridDim.x;
  k = (bottom & (-WARPSIZE)) + threadIdx.x +
      blockIdx.x * blockDim.x; // align to warp size
  if (k < bottom)
    k += inc;

  // iterate over all cells assigned to thread
  while (k < top) {
    childd[k] = -1;
    k += inc;
  }
}

__global__ __launch_bounds__(THREADS2, FACTOR2) void TreeBuildingKernelDipole(
    const int nnodesd, const int nbodiesd,
    volatile int *const __restrict__ childd,
    const float4 *const __restrict__ midpointsd,
    int* bottom_dev, float* radius_dev) {
  int i, j, depth, skip, inc;
  float x, y, z, r;
  float dx, dy, dz;
  int ch, n, cell, locked, patch;
  float radius;

  // cache root data
  radius = (*radius_dev) * 0.5f;
  const float4 root = midpointsd[nnodesd];

  skip = 1;
  inc = blockDim.x * gridDim.x;
  i = threadIdx.x + blockIdx.x * blockDim.x;

  // iterate over all bodies assigned to thread
  while (i < nbodiesd) {
    const float4 p = midpointsd[i];
    if (skip != 0) {
      // new body, so start traversing at root
      skip = 0;
      n = nnodesd;
      depth = 1;
      r = radius;
      dx = dy = dz = -r;
      j = 0;
      // determine which child to follow
      if (root.x < p.x) {
        j = 1;
        dx = r;
      }
      if (root.y < p.y) {
        j |= 2;
        dy = r;
      }
      if (root.z < p.z) {
        j |= 4;
        dz = r;
      }
      x = root.x + dx;
      y = root.y + dy;
      z = root.z + dz;
    }

    // follow path to leaf cell
    ch = childd[n * 8 + j];
    while (ch >= nbodiesd) {
      n = ch;
      depth++;
      r *= 0.5f;
      dx = dy = dz = -r;
      j = 0;
      // determine which child to follow
      if (x < p.x) {
        j = 1;
        dx = r;
      }
      if (y < p.y) {
        j |= 2;
        dy = r;
      }
      if (z < p.z) {
        j |= 4;
        dz = r;
      }
      x += dx;
      y += dy;
      z += dz;
      ch = childd[n * 8 + j];
    }

    if (ch != -2) { // skip if child pointer is locked and try again later
      locked = n * 8 + j;
      if (ch == -1) {
        if (-1 == atomicCAS((int *)&childd[locked], -1,
                            i)) { // if null, just insert the new body
          i += inc;               // move on to next body
          skip = 1;
        }
      } else { // there already is a body at this position
        if (ch == atomicCAS((int *)&childd[locked], ch, -2)) { // try to lock
          patch = -1;
          const float4 chp = midpointsd[ch];
          // create new cell(s) and insert the old and new bodies
          do {
            depth++;
            if (depth > MAXDEPTH) {
              printf("ERROR: maximum depth exceeded (bodies are too close "
                     "together)\n");
              asm("trap;");
            }

            cell = atomicSub(bottom_dev, 1) - 1;
            if (cell <= nbodiesd) {
              printf("ERROR: out of cell memory\n");
              asm("trap;");
            }

            if (patch != -1) {
              childd[n * 8 + j] = cell;
            }
            patch = max(patch, cell);

            j = 0;
            if (x < chp.x)
              j = 1;
            if (y < chp.y)
              j |= 2;
            if (z < chp.z)
              j |= 4;
            childd[cell * 8 + j] = ch;

            n = cell;
            r *= 0.5f;
            dx = dy = dz = -r;
            j = 0;
            if (x < p.x) {
              j = 1;
              dx = r;
            }
            if (y < p.y) {
              j |= 2;
              dy = r;
            }
            if (z < p.z) {
              j |= 4;
              dz = r;
            }
            x += dx;
            y += dy;
            z += dz;

            ch = childd[n * 8 + j];
            // repeat until the two bodies are different children
          } while (ch >= 0);
          childd[n * 8 + j] = i;

          i += inc; // move on to next body
          skip = 2;
        }
      }
    }
    __syncthreads(); // optional barrier for performance
    __threadfence();

    if (skip == 2) {
      childd[locked] = patch;
    }
  }
}

__global__ __launch_bounds__(1024, 1) void ClearKernelDipole2(
    const int nnodesd, int *const __restrict__ startd,
    float4 *const __restrict__ midpointsd, int* bottom_dev) {
  int k, inc, bottom;

  bottom = *bottom_dev;
  inc = blockDim.x * gridDim.x;
  k = (bottom & (-WARPSIZE)) + threadIdx.x +
      blockIdx.x * blockDim.x; // align to warp size
  if (k < bottom)
    k += inc;

  // iterate over all cells assigned to thread
  while (k < nnodesd) {
    midpointsd[k].w = -1.0f;
    startd[k] = -1;
    k += inc;
  }
}

/******************************************************************************/
/*** compute center of mass ***************************************************/
/******************************************************************************/
// just do float4 addition of tangent vector, and use w of tangent for COM comp.
// propagate max for the w component of midpoint vector.
// set mass equal to tangent. <-- just make it equal to w component.
// no point in caching the rest of the tangent.
// let w component of tangent be initialized to -1 for the children
// 12/2/2020: just setting COM to be the center of max and min, and setting
// w to be the max radius.

// DIPOLE: use D1 and D3. x,y,z correspond to rows 1 and 3.
// Row 2 will be D3.w, Tangent.w, D1.w.
// Combining: Parent Dipole = Sum(child_dipole + child_monopole (R_c - R_p)^T)
// Eval: D2 (vector) = <D1.w-D3.y, D3.x-D1.z, D1.y-D3.w>
__global__ __launch_bounds__(THREADS3, FACTOR3) void SummarizationKernelDipole(
    const int nnodesd, const int nbodiesd, const int bottom,
    volatile int *const __restrict__ countd,
    const int *const __restrict__ childd, float4 *const __restrict__ midpointsd,
    float4 *const __restrict__ tangentsd, float4 *const __restrict__ dipole1sd,
    float4 *const __restrict__ dipole3sd) {
  int i, j, k, ch, inc, cnt;
  float xmin, ymin, zmin, xmax, ymax, zmax;
  float maxw;
  float4 midpoint, tangent, ctangent;
  float4 dipole1, dipole3, center, cdipole1, cdipole3;
  __shared__ int child[THREADS3 * 8];
  // __shared__ float4 s_midpt[THREADS3 * 8];
  // __shared__ float4 s_tangent[THREADS3 * 8];
  // this is the tangent vector's w component.
  //__shared__ float mass[THREADS3 * 8];

  inc = blockDim.x * gridDim.x;
  k = (bottom & (-WARPSIZE)) + threadIdx.x +
      blockIdx.x * blockDim.x; // align to warp size
  if (k < bottom)
    k += inc;

  int restart = k;
  for (j = 0; j < 3; j++) { // wait-free pre-passes
    // iterate over all cells assigned to thread
    while (k <= nnodesd) {
      if (midpointsd[k].w < 0.0f) {
        for (i = 0; i < 8; i++) {
          ch = childd[k * 8 + i];
          child[i * THREADS3 + threadIdx.x] = ch; // cache children
          if ((ch >= nbodiesd) && (midpointsd[ch].w < 0.0f)) {
            break;
          }
        }
        if (i == 8) {
          // all children are ready
          // cm = 0.0f;
          ctangent = {0.f, 0.f, 0.f, 0.f};
          cdipole1 = {0.f, 0.f, 0.f, 0.f};
          cdipole3 = {0.f, 0.f, 0.f, 0.f};
          xmin = CUDART_INF_F;
          ymin = CUDART_INF_F;
          zmin = CUDART_INF_F;
          xmax = -CUDART_INF_F;
          ymax = -CUDART_INF_F;
          zmax = -CUDART_INF_F;
          maxw = 0.0f;
          cnt = 0;
          for (i = 0; i < 8; i++) {
            ch = child[i * THREADS3 + threadIdx.x];
            if (ch >= 0) {
              // four reads due to missing copy constructor for "volatile
              // float4"
              midpoint = midpointsd[ch];
              tangent = tangentsd[ch];
              // const float chx = midpointsd[ch].x;
              // const float chy = midpointsd[ch].y;
              // const float chz = midpointsd[ch].z;
              // const float chw = midpointsd[ch].w;
              if (ch >= nbodiesd) { // count bodies (needed later)
                cnt += countd[ch];
              } else {
                cnt++;
              }
              // add child's contribution
              // cm += tangent.w;
              ctangent += tangent;
              xmin = fminf(xmin, midpoint.x);
              ymin = fminf(ymin, midpoint.y);
              zmin = fminf(zmin, midpoint.z);
              xmax = fmaxf(xmax, midpoint.x);
              ymax = fmaxf(ymax, midpoint.y);
              zmax = fmaxf(zmax, midpoint.z);
              //px += midpoint.x * tangent.w;
              //py += midpoint.y * tangent.w;
              //pz += midpoint.z * tangent.w;
              maxw = fmaxf(maxw, midpoint.w);
            }
          }
          //m = 1.0f / ctangent.w;
          center = make_float4(
            0.5f * (xmin + xmax),
            0.5f * (ymin + ymax),
            0.5f * (zmin + zmax),
            maxw);
          /*ctangent.w += 0.25f * ((xmax - xmin) * (xmax - xmin)
                                  + (ymax - ymin) * (ymax - ymin)
                                  + (zmax - zmin) * (zmax - zmin));*/
          

          for (i = 0; i < 8; i++) {
            ch = child[i * THREADS3 + threadIdx.x];
            if (ch >= 0) {
              // four reads due to missing copy constructor for "volatile
              // float4"
              midpoint = midpointsd[ch] - center;
              tangent = tangentsd[ch];
              dipole1 = dipole1sd[ch];
              dipole3 = dipole3sd[ch];
              // add D3 of the child to parent
              cdipole3 += dipole3;
              cdipole1 += dipole1;
              ctangent.w += tangent.w;
              // Add D1C r_c^T of the child to parent
              cdipole1.x += tangent.x * midpoint.x;
              cdipole1.y += tangent.x * midpoint.y;
              cdipole1.z += tangent.x * midpoint.z;
              cdipole3.x += tangent.z * midpoint.x;
              cdipole3.y += tangent.z * midpoint.y;
              cdipole3.z += tangent.z * midpoint.z;
              cdipole3.w += tangent.y * midpoint.x;
              ctangent.w += tangent.y * midpoint.y;
              cdipole1.w += tangent.y * midpoint.z;
            }
          }
          __threadfence();
          countd[k] = cnt;
          midpointsd[k] = center;
          tangentsd[k] = ctangent;
          dipole1sd[k] = cdipole1;
          dipole3sd[k] = cdipole3;
        }
      }
      k += inc; // move on to next cell
    }
    k = restart;
  }

  j = 0;
  // iterate over all cells assigned to thread
  while (k <= nnodesd) {
    if (midpointsd[k].w >= 0.0f) {
      k += inc;
    } else {
      if (j == 0) {
        j = 8;
        for (i = 0; i < 8; i++) {
          ch = childd[k * 8 + i];
          child[i * THREADS3 + threadIdx.x] = ch; // cache children
          if ((ch < nbodiesd) || ((midpointsd[ch].w) >= 0.0f)) {
            j--;
          }
        }
      } else {
        j = 8;
        for (i = 0; i < 8; i++) {
          ch = child[i * THREADS3 + threadIdx.x];
          if ((ch < nbodiesd) || ((midpointsd[ch].w) >= 0.0f)) {
            j--;
          }
        }
      }

      if (j == 0) {
        // all children are ready
        ctangent = {0.f, 0.f, 0.f, 0.f};
        cdipole1 = {0.f, 0.f, 0.f, 0.f};
        cdipole3 = {0.f, 0.f, 0.f, 0.f};
        xmin = CUDART_INF_F;
        ymin = CUDART_INF_F;
        zmin = CUDART_INF_F;
        xmax = -CUDART_INF_F;
        ymax = -CUDART_INF_F;
        zmax = -CUDART_INF_F;
        maxw = 0.0f;
        cnt = 0;
        for (i = 0; i < 8; i++) {
          ch = child[i * THREADS3 + threadIdx.x];
          if (ch >= 0) {
            // four reads due to missing copy constructor for "volatile float4"
            midpoint = midpointsd[ch];
            tangent = tangentsd[ch];
            // const float chx = midpointsd[ch].x;
            // const float chy = midpointsd[ch].y;
            // const float chz = midpointsd[ch].z;
            // const float chw = midpointsd[ch].w;
            if (ch >= nbodiesd) { // count bodies (needed later)
              cnt += countd[ch];
            } else {
              cnt++;
            }
            // add child's contribution
            ctangent += tangent;
            xmin = fminf(xmin, midpoint.x);
            ymin = fminf(ymin, midpoint.y);
            zmin = fminf(zmin, midpoint.z);
            xmax = fmaxf(xmax, midpoint.x);
            ymax = fmaxf(ymax, midpoint.y);
            zmax = fmaxf(zmax, midpoint.z);
            maxw = fmaxf(maxw, midpoint.w);
          }
        }
        //m = 1.0f / ctangent.w;
        center = make_float4(
          0.5f * (xmin + xmax),
          0.5f * (ymin + ymax),
          0.5f * (zmin + zmax),
          maxw);
        /*ctangent.w += 0.25f * ((xmax - xmin) * (xmax - xmin)
                                + (ymax - ymin) * (ymax - ymin)
                                + (zmax - zmin) * (zmax - zmin));*/
        for (i = 0; i < 8; i++) {
          ch = child[i * THREADS3 + threadIdx.x];
          if (ch >= 0) {
            // four reads due to missing copy constructor for "volatile
            // float4"
            midpoint = midpointsd[ch] - center;
            tangent = tangentsd[ch];
            dipole1 = dipole1sd[ch];
            dipole3 = dipole3sd[ch];
            // add D3 of the child to parent
            cdipole3 += dipole3;
            cdipole1 += dipole1;
            ctangent.w += tangent.w;
            // Add D1C r_c^T of the child to parent
            cdipole1.x += tangent.x * midpoint.x;
            cdipole1.y += tangent.x * midpoint.y;
            cdipole1.z += tangent.x * midpoint.z;
            cdipole3.x += tangent.z * midpoint.x;
            cdipole3.y += tangent.z * midpoint.y;
            cdipole3.z += tangent.z * midpoint.z;
            cdipole3.w += tangent.y * midpoint.x;
            ctangent.w += tangent.y * midpoint.y;
            cdipole1.w += tangent.y * midpoint.z;
          }
        }

        __threadfence();
        countd[k] = cnt;
        midpointsd[k] = center;
        tangentsd[k] = ctangent;
        dipole1sd[k] = cdipole1;
        dipole3sd[k] = cdipole3;
        k += inc;
      }
    }
  }
}

/******************************************************************************/
/*** sort bodies **************************************************************/
/******************************************************************************/

__global__ __launch_bounds__(THREADS4, FACTOR4) void SortKernelDipole(
    const int nnodesd, const int nbodiesd, const int bottom,
    int *const __restrict__ sortd, const int *const __restrict__ countd,
    volatile int *const __restrict__ startd, int *const __restrict__ childd) {
  int i, j, k, ch, dec, start;

  dec = blockDim.x * gridDim.x;
  k = nnodesd + 1 - dec + threadIdx.x + blockIdx.x * blockDim.x;

  // iterate over all cells assigned to thread
  while (k >= bottom) {
    start = startd[k];
    if (start >= 0) {
      j = 0;
      for (i = 0; i < 8; i++) {
        ch = childd[k * 8 + i];
        if (ch >= 0) {
          if (i != j) {
            // move children to front (needed later for speed)
            childd[k * 8 + i] = -1;
            childd[k * 8 + j] = ch;
          }
          j++;
          if (ch >= nbodiesd) {
            // child is a cell
            startd[ch] = start;  // set start ID of child
            start += countd[ch]; // add #bodies in subtree
          } else {
            // child is a body
            sortd[start] = ch; // record body in 'sorted' array
            start++;
          }
        }
      }
      k -= dec; // move on to next cell
    }
    __syncthreads(); // optional barrier for performance
  }
}


/******************************************************************************/
/*** compute linkages
 * ************************************************************/
/******************************************************************************/

__global__ __launch_bounds__(THREADS5, FACTOR5) void LinkageCalculationKernelDipole(
    const int nnodesds, const int nbodiesds, const int nbodiesdt, const float radiusd,
    const float itolsqd, const float epssqd, const float distratiosqd,
    const int *const __restrict__ childds,
    const float4 *const __restrict__ midpointsds,
    const float4 *const __restrict__ midpointsdt,
    const float4 *const __restrict__ tangentsds,
    const float4 *const __restrict__ tangentsdt,
    const float4 *const __restrict__ dipole1sds,
    const float4 *const __restrict__ dipole3sds,
    float *const __restrict__ linkagest) {
  int i, j, k, n, depth, base, sbase, diff, pd, nd;
  float lk, dx, dy, dz, tmp, distth;
  float3 alpha, beta, gamma, delta, dr, halftn, halfti;
  float la, lb, lc, ld, div1, div2, ac;
  __shared__ volatile int pos[MAXDEPTH * THREADS5 / WARPSIZE],
      node[MAXDEPTH * THREADS5 / WARPSIZE];
  __shared__ float dq[MAXDEPTH * THREADS5 / WARPSIZE];

  if (0 == threadIdx.x) {
    tmp = radiusd;
    // precompute values that depend only on tree level
    dq[0] = 0.25f * tmp * tmp * itolsqd; // 0.25 since we're using the children of each node.
    for (i = 1; i < MAXDEPTH; i++) {
      dq[i] = dq[i - 1] * 0.25f;
      dq[i - 1] += epssqd;
    }
    dq[i - 1] += epssqd;
  }
  __syncthreads();

  // figure out first thread in each warp (lane 0)
  base = threadIdx.x / WARPSIZE;
  sbase = base * WARPSIZE;
  j = base * MAXDEPTH;

  diff = threadIdx.x - sbase;
  // make multiple copies to avoid index calculations later
  if (diff < MAXDEPTH) {
    dq[diff + j] = dq[diff];
  }
  __syncthreads();

  // iterate over all bodies assigned to thread
  for (k = threadIdx.x + blockIdx.x * blockDim.x; k < nbodiesdt;
       k += blockDim.x * gridDim.x) {
    // i = sortdt[k]; // get permuted/sorted index
    i = k; // for linkages this is faster.
    // cache position info
    // target point
    const float4 pi = midpointsdt[i];
    // this is C1
    const float4 ti = tangentsdt[i];

    lk = 0.0f;

    // initialize iteration stack, i.e., push root node onto stack
    depth = j;
    if (sbase == threadIdx.x) {
      pos[j] = 0;
      node[j] = nnodesds * 8;
    }

    do {
      // stack is not empty
      pd = pos[depth];
      nd = node[depth];
      while (pd < 8) {
        // node on top of stack has more children to process
        n = childds[nd + pd]; // load child pointer
        pd++;

        if (n >= 0) {
          // source point
          const float4 pn = midpointsds[n];
          // this is D1
          const float4 tn = tangentsds[n];
          dx = pn.x - pi.x;
          dy = pn.y - pi.y;
          dz = pn.z - pi.z;
          tmp =
              dx * dx +
              (dy * dy +
               (dz * dz + epssqd)); // compute distance squared (plus softening)
          // Compute distance threshold
          distth = distratiosqd * (pn.w + pi.w + 2.f * sqrtf(pn.w * pi.w));
          distth = dq[depth] + distth + 2.f * sqrtf(distth * dq[depth]);
          // distth = fmaxf(8.f * tn.w * itolsqd, distratiosqd * (pn.w + pi.w));
          // n < nbodiesd checks if it's a body. the latter checks if cell is
          // far enough away.
          // direct (body) eval
          if (n < nbodiesds || __all_sync(0xffffffff, (tmp >= distth))) {
            // Store determinant into distth.
            // this is (D1 X C1) dot r0.
            distth = (tn.y * ti.z - tn.z * ti.y) * dx +
                     (tn.z * ti.x - tn.x * ti.z) * dy +
                     (tn.x * ti.y - tn.y * ti.x) * dz;
            if (n < nbodiesds) {
              // compute distance with arctan correction
              // denominator
              dr = make_float3(dx, dy, dz);
              // segment 2
              halftn = make_float3(0.5f * tn.x, 0.5f * tn.y, 0.5f * tn.z);
              // segment 1
              halfti = make_float3(0.5f * ti.x, 0.5f * ti.y, 0.5f * ti.z);
              alpha = dr - halftn + halfti; // start of 2 - start of 1
              beta = dr - halftn - halfti;  // start of 2 - end of 1
              gamma = dr + halftn - halfti; // end of 2 - end of 1
              la = length(alpha);
              lb = length(beta);
              lc = length(gamma);
              ac = dot(alpha, gamma);
              div1 = la * lb * lc + dot(alpha, beta) * lc + ac * lb +
                     dot(beta, gamma) * la;
              delta = dr + halftn + halfti; // end of 2 - start of 1
              ld = length(delta);
              div2 = la * ld * lc + dot(alpha, delta) * lc + ac * ld +
                     dot(delta, gamma) * la;
              // remember, determinant is in distth.
              // x'
              la = div1 * div2 - distth * distth;
              // ac is s'' = y' * determinant. We also want it negative when det = 0 and div1, div2 are both negative.
              ac = (distth == 0) ? copysignf(1.f, div1) + copysignf(1.f, div2) : (div1 + div2); // = y' * distth
              // y'
              lb = distth * ac;
              lc = (distth == 0) ? -copysignf(1.f, div1) : copysignf(1.f, distth);
              lk += (2.f * INV_FOUR_PI) * atan2f(lb, la) + (((ac < 0) || ((ac == 0) && (distth < 0))) ? lc : 0);

              //lk += (2.f * INV_FOUR_PI) *
              //      (atan2f(distth, div1) + atan2f(distth, div2));
            } else {
              const float4 dipole1n = dipole1sds[n];
              const float4 dipole3n = dipole3sds[n];
              div2 = 1./tmp;
              tmp = rsqrtf(tmp); // compute distance
              tmp = tmp * div2 * (INV_FOUR_PI);
              // div2 is |r_0|^-2, tmp is (4 pi)^(-1) |r_0|^-3
              // remember, determinant is in distth.
              //lk += tmp * distth;
              // compute D2
              // D2 (vector) = <D1.w-D3.y, D3.x-D1.z, D1.y-D3.w>
              delta.x = dipole1n.w - dipole3n.y;
              delta.y = dipole3n.x - dipole1n.z;
              delta.z = dipole1n.y - dipole3n.w;
              // compute C1 dot D2 
              distth -= (ti.x * delta.x + ti.y * delta.y + ti.z * delta.z);
              // compute D3 r0
              delta.x = dipole1n.x * dx + dipole1n.y * dy + dipole1n.z * dz; 
              delta.y = dipole3n.w * dx + tn.w * dy + dipole1n.w * dz;
              delta.z = dipole3n.x * dx + dipole3n.y * dy + dipole3n.z * dz;
              // compute 3 * (D3r0 cross C1) dot r0 / |r|^2.
              distth -= 3.f * div2 * ((delta.y * ti.z - delta.z * ti.y) * dx +
                       (delta.z * ti.x - delta.x * ti.z) * dy +
                       (delta.x * ti.y - delta.y * ti.x) * dz);
              // tmp is (4 pi)^-1 |r_0|^-3
              lk += tmp * distth;
            }
          } else {
            // push cell onto stack
            if (sbase == threadIdx.x) {
              pos[depth] = pd;
              node[depth] = nd;
            }
            depth++;
            pd = 0;
            nd = n * 8;
          }
        } else {
          pd = 8; // early out because all remaining children are also zero
        }
      }
      depth--; // done with this level
    } while (depth >= j);

    /*float4 acc = accVeld[i];
    if (stepd > 0) {
      // update velocity
      float2 v = veld[i];
      v.x += (ax - acc.x) * dthfd;
      v.y += (ay - acc.y) * dthfd;
      acc.w += (az - acc.z) * dthfd;
      veld[i] = v;
    }

    // save computed acceleration
    acc.x = ax;
    acc.y = ay;
    acc.z = az;
    accVeld[i] = acc; */
    linkagest[i] = lk;
  }
}


/******************************************************************************/
/*** compute writhe
 * ************************************************************/
/******************************************************************************/

__global__ __launch_bounds__(THREADS5, FACTOR5) void WritheCalculationKernelDipole(
  const int nnodesds, const int nbodiesds, const float radiusd,
  const float itolsqd, const float epssqd, const float distratiosqd,
  const int *const __restrict__ sortdt, const int *const __restrict__ childds,
  const float4 *const __restrict__ midpointsds,
  const float4 *const __restrict__ tangentsds,
  const float4 *const __restrict__ dipole1sds,
  const float4 *const __restrict__ dipole3sds,
  float *const __restrict__ writhest, const bool use_sorted_order) {
int i, j, k, n, depth, base, sbase, diff, pd, nd;
float lk, dx, dy, dz, tmp, distth;
float3 alpha, beta, gamma, delta, dr, halftn, halfti;
float la, lb, lc, ld, div1, div2, ac;
__shared__ volatile int pos[MAXDEPTH * THREADS5 / WARPSIZE],
    node[MAXDEPTH * THREADS5 / WARPSIZE];
__shared__ float dq[MAXDEPTH * THREADS5 / WARPSIZE];

if (0 == threadIdx.x) {
  tmp = radiusd;
  // precompute values that depend only on tree level
  dq[0] = 0.25f * tmp * tmp * itolsqd; // 0.25 since we're using the children of each node.
  for (i = 1; i < MAXDEPTH; i++) {
    dq[i] = dq[i - 1] * 0.25f;
    dq[i - 1] += epssqd;
  }
  dq[i - 1] += epssqd;
}
__syncthreads();

// figure out first thread in each warp (lane 0)
base = threadIdx.x / WARPSIZE;
sbase = base * WARPSIZE;
j = base * MAXDEPTH;

diff = threadIdx.x - sbase;
// make multiple copies to avoid index calculations later
if (diff < MAXDEPTH) {
  dq[diff + j] = dq[diff];
}
__syncthreads();

// iterate over all bodies assigned to thread
for (k = threadIdx.x + blockIdx.x * blockDim.x; k < nbodiesds;
     k += blockDim.x * gridDim.x) {
  if (use_sorted_order) {
    i = sortdt[k]; // get permuted/sorted index
  } else {
    i = k; // for linkages this is faster.
  }
  // cache position info
  // target point
  const float4 pi = midpointsds[i];
  // this is C1
  const float4 ti = tangentsds[i];

  lk = 0.0f;

  // initialize iteration stack, i.e., push root node onto stack
  depth = j;
  if (sbase == threadIdx.x) {
    pos[j] = 0;
    node[j] = nnodesds * 8;
  }

  do {
    // stack is not empty
    pd = pos[depth];
    nd = node[depth];
    while (pd < 8) {
      // node on top of stack has more children to process
      n = childds[nd + pd]; // load child pointer
      pd++;

      if (n >= 0) {
        // source point
        const float4 pn = midpointsds[n];
        // this is D1
        const float4 tn = tangentsds[n];
        dx = pn.x - pi.x;
        dy = pn.y - pi.y;
        dz = pn.z - pi.z;
        tmp =
            dx * dx +
            (dy * dy +
             (dz * dz + epssqd)); // compute distance squared (plus softening)
        // Compute distance threshold
        distth = distratiosqd * (pn.w + pi.w + 2.f * sqrtf(pn.w * pi.w));
        distth = dq[depth] + distth + 2.f * sqrtf(distth * dq[depth]);
        // distth = fmaxf(8.f * tn.w * itolsqd, distratiosqd * (pn.w + pi.w));
        // n < nbodiesd checks if it's a body. the latter checks if cell is
        // far enough away.
        // direct (body) eval
        if (n < nbodiesds || __all_sync(0xffffffff, (tmp >= distth))) {
          // Store determinant into distth.
          // this is (D1 X C1) dot r0.
          distth = (tn.y * ti.z - tn.z * ti.y) * dx +
                   (tn.z * ti.x - tn.x * ti.z) * dy +
                   (tn.x * ti.y - tn.y * ti.x) * dz;
          if (n < nbodiesds) {
            // compute distance with arctan correction
            // denominator
            dr = make_float3(dx, dy, dz);
            // segment 2
            halftn = make_float3(0.5f * tn.x, 0.5f * tn.y, 0.5f * tn.z);
            // segment 1
            halfti = make_float3(0.5f * ti.x, 0.5f * ti.y, 0.5f * ti.z);
            alpha = dr - halftn + halfti; // start of 2 - start of 1
            beta = dr - halftn - halfti;  // start of 2 - end of 1
            gamma = dr + halftn - halfti; // end of 2 - end of 1
            delta = dr + halftn + halfti; // end of 2 - start of 1
            la = length(alpha);
            lb = length(beta);
            lc = length(gamma);
            ac = dot(alpha, gamma);
            div1 = la * lb * lc + dot(alpha, beta) * lc + ac * lb +
                   dot(beta, gamma) * la;
            ld = length(delta);
            div2 = la * ld * lc + dot(alpha, delta) * lc + ac * ld +
                   dot(delta, gamma) * la;
            // remember, determinant is in distth.
            // x'
            la = div1 * div2 - distth * distth;
            // ac is s'' = y' * determinant. We also want it negative when det = 0 and div1, div2 are both negative.
            ac = (distth == 0) ? copysignf(1.f, div1) + copysignf(1.f, div2) : (div1 + div2); // = y' * distth
            // y'
            lb = distth * ac;
            lc = (distth == 0) ? -copysignf(1.f, div1) : copysignf(1.f, distth);
            lk += (2.f * INV_FOUR_PI) * atan2f(lb, la) + (((ac < 0) || ((ac == 0) && (distth < 0))) ? lc : 0);
            //lk += (2.f * INV_FOUR_PI) *
            //      (atan2f(distth, div1) + atan2f(distth, div2));
          } else {
            const float4 dipole1n = dipole1sds[n];
            const float4 dipole3n = dipole3sds[n];
            div2 = 1./tmp;
            tmp = rsqrtf(tmp); // compute distance
            tmp = tmp * div2 * (INV_FOUR_PI);
            // div2 is |r_0|^-2, tmp is (4 pi)^(-1) |r_0|^-3
            // remember, determinant is in distth.
            //lk += tmp * distth;
            // compute D2
            // D2 (vector) = <D1.w-D3.y, D3.x-D1.z, D1.y-D3.w>
            delta.x = dipole1n.w - dipole3n.y;
            delta.y = dipole3n.x - dipole1n.z;
            delta.z = dipole1n.y - dipole3n.w;
            // compute C1 dot D2 
            distth -= (ti.x * delta.x + ti.y * delta.y + ti.z * delta.z);
            // compute D3 r0
            delta.x = dipole1n.x * dx + dipole1n.y * dy + dipole1n.z * dz; 
            delta.y = dipole3n.w * dx + tn.w * dy + dipole1n.w * dz;
            delta.z = dipole3n.x * dx + dipole3n.y * dy + dipole3n.z * dz;
            // compute 3 * (D3r0 cross C1) dot r0 / |r|^2.
            distth -= 3.f * div2 * ((delta.y * ti.z - delta.z * ti.y) * dx +
                     (delta.z * ti.x - delta.x * ti.z) * dy +
                     (delta.x * ti.y - delta.y * ti.x) * dz);
            // tmp is (4 pi)^-1 |r_0|^-3
            lk += tmp * distth;
          }
        } else {
          // push cell onto stack
          if (sbase == threadIdx.x) {
            pos[depth] = pd;
            node[depth] = nd;
          }
          depth++;
          pd = 0;
          nd = n * 8;
        }
      } else {
        pd = 8; // early out because all remaining children are also zero
      }
    }
    depth--; // done with this level
  } while (depth >= j);

  /*float4 acc = accVeld[i];
  if (stepd > 0) {
    // update velocity
    float2 v = veld[i];
    v.x += (ax - acc.x) * dthfd;
    v.y += (ay - acc.y) * dthfd;
    acc.w += (az - acc.z) * dthfd;
    veld[i] = v;
  }

  // save computed acceleration
  acc.x = ax;
  acc.y = ay;
  acc.z = az;
  accVeld[i] = acc; */
  writhest[i] = lk;
}
}


/******************************************************************************/

static void CudaTest(const char *const msg) {
  cudaError_t e;

  cudaDeviceSynchronize();
  if (cudaSuccess != (e = cudaGetLastError())) {
    fprintf(stderr, "%s: %d\n", msg, e);
    fprintf(stderr, "%s\n", cudaGetErrorString(e));
    throw std::runtime_error("CudaTest failure.");
  }
}

/******************************************************************************/


// computes the w coordinates for each yarn
void PrepareYarn(std::vector<float4> &midpts,
                 std::vector<float4> &tangents,
                 std::vector<float4> &dipole1s,
                 std::vector<float4> &dipole3s) {
  if (midpts.size() != tangents.size()) {
    std::cout << "Preparing yarn: midpoint size and tangent size are different."
              << std::endl;
    throw std::runtime_error("PreparingYarn size error.");
  }
  for (int i = 0; i < tangents.size(); ++i) {
    using namespace std;
    midpts[i].w = tangents[i].x * tangents[i].x +
                  tangents[i].y * tangents[i].y + tangents[i].z * tangents[i].z;
    tangents[i].w = 0.f;
    dipole1s[i] = {0.f,0.f,0.f,0.f};
    dipole3s[i] = {0.f,0.f,0.f,0.f};
  }
}
bool checked_dev_properties_d = false;
void CheckDevicePropertiesForDevDipole(const cudaDeviceProp& deviceProp) {
  if ((deviceProp.major == 9999) && (deviceProp.minor == 9999)) {
    fprintf(stderr, "There is no CUDA capable device\n");
    throw std::runtime_error("PrepareTree CUDA device error.");
  }
  if (deviceProp.major < 3) {
    fprintf(stderr, "Need at least compute capability 3.0\n");
    throw std::runtime_error("PrepareTree CUDA cc error.");
  }
  if (deviceProp.warpSize != WARPSIZE) {
    fprintf(stderr, "Warp size must be %d\n", deviceProp.warpSize);
    throw std::runtime_error("PrepareTree CUDA warp size error.");
  }

  int blocks = deviceProp.multiProcessorCount;
  const int mTSM = deviceProp.maxThreadsPerMultiProcessor;
  printf("gpu: %s with %d SMs and %d mTpSM (%.1f MHz and %.1f MHz)\n",
         deviceProp.name, blocks, mTSM, deviceProp.clockRate * 0.001,
         deviceProp.memoryClockRate * 0.001);

  if ((WARPSIZE <= 0) || (WARPSIZE & (WARPSIZE - 1) != 0)) {
    fprintf(stderr, "Warp size must be greater than zero and a power of two\n");
    throw std::runtime_error("PrepareTree warp size error.");
  }
  if (MAXDEPTH > WARPSIZE) {
    fprintf(stderr, "MAXDEPTH must be less than or equal to WARPSIZE\n");
    throw std::runtime_error("PrepareTree max depth error.");
  }
  if ((THREADS1 <= 0) || (THREADS1 & (THREADS1 - 1) != 0)) {
    fprintf(stderr, "THREADS1 must be greater than zero and a power of two\n");
    throw std::runtime_error("PrepareTree THREADS1 error.");
  }
  checked_dev_properties_d = true;
}

// Call this before preparing the dipole trees.
void AllocateYarnRadiiAndBottoms(float** radii, int** bottoms, int* blocks, const int nyarns) {
  cudaGetLastError(); // reset error value
  if (cudaSuccess != cudaMalloc((void **)radii, sizeof(float) * (nyarns)))
    fprintf(stderr, "could not allocate radii\n");
  CudaTest("couldn't allocate radii");
  if (cudaSuccess != cudaMalloc((void **)bottoms, sizeof(int) * (nyarns)))
    fprintf(stderr, "could not allocate bottoms\n");
  CudaTest("couldn't allocate bottoms");
  cudaFuncSetCacheConfig(BoundingBoxKernelDipole, cudaFuncCachePreferShared);
  cudaFuncSetCacheConfig(TreeBuildingKernelDipole, cudaFuncCachePreferL1);
  cudaFuncSetCacheConfig(ClearKernelDipole1, cudaFuncCachePreferL1);
  cudaFuncSetCacheConfig(ClearKernelDipole2, cudaFuncCachePreferL1);
  cudaFuncSetCacheConfig(SummarizationKernelDipole, cudaFuncCachePreferShared);
  cudaFuncSetCacheConfig(SortKernelDipole, cudaFuncCachePreferL1);
  cudaFuncSetCacheConfig(LinkageCalculationKernelDipole, cudaFuncCachePreferEqual);
  //TreeMeta cpu_tree_meta;
  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, 0 /* device number */);
  if (!checked_dev_properties_d) CheckDevicePropertiesForDevDipole(deviceProp);

  *blocks = deviceProp.multiProcessorCount;
}

// returns a <midpts_dev, tangents_dev, sortl, childl, radius, nnodes> for the
// tree
YarnDipoleTree PrepareDipoleTree(const std::vector<float4> &midpts,
                     const std::vector<float4> &tangents,
                     const std::vector<float4> &dipole1s,
                     const std::vector<float4> &dipole3s,
                     float* radius_dev, int* bottom_dev,
                     const int blocks, const int min_threshold) {
  int nnodes, nbodies;
  float4 *midpts_dev;
  float4 *tangents_dev;
  float4 *dipole1s_dev = nullptr;
  float4 *dipole3s_dev = nullptr;
  int *sortl = nullptr;
  int *childl = nullptr;
  int *countl = nullptr;
  int *startl = nullptr;
  float3 *minl, *maxl;
  
  cudaGetLastError(); // reset error value
  nbodies = midpts.size();
  nnodes = nbodies * 2.5; // experimentally, 2.1 is enough.
  if (nnodes < 1024 * blocks)
    nnodes = 1024 * blocks;
  while ((nnodes & (WARPSIZE - 1)) != 0)
    nnodes++;
  nnodes--;

  if (cudaSuccess !=
    cudaMalloc((void **)&midpts_dev, sizeof(float4) * (nnodes + 1)))
  fprintf(stderr, "could not allocate midptsd\n");
  CudaTest("couldn't allocate midptsd");
  if (cudaSuccess !=
      cudaMalloc((void **)&tangents_dev, sizeof(float4) * (nnodes + 1)))
    fprintf(stderr, "could not allocate tangentsd\n");
  CudaTest("couldn't allocate tangentsd");
  if (cudaSuccess != cudaMemcpy(midpts_dev, midpts.data(),
                                sizeof(float4) * nbodies,
                                cudaMemcpyHostToDevice))
    fprintf(stderr, "copying of midptsd to device failed\n");
  CudaTest("midptsd copy to device failed");
  if (cudaSuccess != cudaMemcpy(tangents_dev, tangents.data(),
                                sizeof(float4) * nbodies,
                                cudaMemcpyHostToDevice))
    fprintf(stderr, "copying of tangentd to device failed\n");
  CudaTest("tangentd copy to device failed");

  if (nbodies >= min_threshold) {

    if (cudaSuccess !=
      cudaMalloc((void **)&childl, sizeof(int) * (nnodes + 1) * 8))
      fprintf(stderr, "could not allocate childd\n");
    CudaTest("couldn't allocate childd");
    if (cudaSuccess != cudaMalloc((void **)&countl, sizeof(int) * (nnodes + 1)))
      fprintf(stderr, "could not allocate countd\n");
    CudaTest("couldn't allocate countd");
    if (cudaSuccess != cudaMalloc((void **)&startl, sizeof(int) * (nnodes + 1)))
      fprintf(stderr, "could not allocate startd\n");
    CudaTest("couldn't allocate startd");
    if (cudaSuccess != cudaMalloc((void **)&sortl, sizeof(int) * (nnodes + 1)))
      fprintf(stderr, "could not allocate sortd\n");
    CudaTest("couldn't allocate sortd");

    if (cudaSuccess !=
        cudaMalloc((void **)&dipole1s_dev, sizeof(float4) * (nnodes + 1)))
      fprintf(stderr, "could not allocate dipole1s\n");
    CudaTest("couldn't allocate dipole1s");
    if (cudaSuccess !=
        cudaMalloc((void **)&dipole3s_dev, sizeof(float4) * (nnodes + 1)))
      fprintf(stderr, "could not allocate dipole3s\n");
    CudaTest("couldn't allocate dipole3s");
    if (cudaSuccess !=
        cudaMalloc((void **)&maxl, sizeof(float3) * blocks * FACTOR1))
      fprintf(stderr, "could not allocate maxd\n");
    CudaTest("couldn't allocate maxd");
    if (cudaSuccess !=
        cudaMalloc((void **)&minl, sizeof(float3) * blocks * FACTOR1))
      fprintf(stderr, "could not allocate mind\n");
    CudaTest("couldn't allocate mind");
    if (cudaSuccess != cudaMemcpy(dipole1s_dev, dipole1s.data(),
                                  sizeof(float4) * nbodies,
                                  cudaMemcpyHostToDevice))
      fprintf(stderr, "copying of dipole1 to device failed\n");
    CudaTest("dipole1 copy to device failed");
    if (cudaSuccess != cudaMemcpy(dipole3s_dev, dipole3s.data(),
                                  sizeof(float4) * nbodies,
                                  cudaMemcpyHostToDevice))
      fprintf(stderr, "copying of dipole3 to device failed\n");
    CudaTest("dipole3 copy to device failed");
    // Launch kernels.
    InitializationKernelDipole<<<1, 1>>>();
    // Kernel 1
    BoundingBoxKernelDipole<<<blocks * FACTOR1, THREADS1>>>(
        nnodes, nbodies, startl, childl, midpts_dev, maxl, minl, bottom_dev, radius_dev);
    // cudaDeviceSynchronize();
    // Kernel 2
    ClearKernelDipole1<<<blocks * 1, 1024>>>(nnodes, nbodies, childl);
    TreeBuildingKernelDipole<<<blocks * FACTOR2, THREADS2>>>(nnodes, nbodies, childl,
                                                      midpts_dev, bottom_dev, radius_dev);
    ClearKernelDipole2<<<blocks * 1, 1024>>>(nnodes, startl, midpts_dev, bottom_dev);
    // cudaDeviceSynchronize();
    // download data from Kernel 2


    // cudaFree(countl);
    // cudaFree(startl);
    cudaFree(minl);
    cudaFree(maxl);
  }
  return {midpts_dev, tangents_dev, dipole1s_dev, dipole3s_dev,
          sortl,      childl,   countl, startl,    0, nnodes};
}

// call this after PrepareDipoleTree.
void GrabYarnRadiiAndBottoms(std::vector<YarnDipoleTree>& yarn_dipole_trees,
    std::vector<int>& bottoms, float* radii_dev, int* bottoms_dev) {
  if (yarn_dipole_trees.size() != bottoms.size()) {
    throw std::runtime_error("GrabYarnRadiiAndBottoms yarn dipole trees and bottoms sizes are different.");
  }
  std::vector<float> radii(bottoms.size());
  if (cudaSuccess != cudaMemcpy(radii.data(), radii_dev,
                                sizeof(float) * radii.size(),
                                cudaMemcpyDeviceToHost))
    fprintf(stderr, "copying of radii to host failed\n");
  CudaTest("radii copy to host failed");
  if (cudaSuccess != cudaMemcpy(bottoms.data(), bottoms_dev,
                                sizeof(int) * bottoms.size(),
                                cudaMemcpyDeviceToHost))
    fprintf(stderr, "copying of bottoms to host failed\n");
  CudaTest("bottoms copy to host failed");
  for (int i = 0; i < radii.size(); ++i) {
    yarn_dipole_trees[i].radius = radii[i];
  }
  cudaFree(radii_dev);
  cudaFree(bottoms_dev);
}

// call this after PrepareDipoleTree and after GrabYarnRadiiAndBottoms
void ComputeKernel3(YarnDipoleTree& yarn_dipole_tree, const int nbodies, const int bottom, const int blocks) {
  if (yarn_dipole_tree.childl == nullptr) return;
  // Kernel 3
  SummarizationKernelDipole<<<blocks * FACTOR3, THREADS3>>>(
    yarn_dipole_tree.nnodes, nbodies, bottom, yarn_dipole_tree.countl, yarn_dipole_tree.childl, yarn_dipole_tree.midpts,
    yarn_dipole_tree.tangents, yarn_dipole_tree.dipole1s, yarn_dipole_tree.dipole3s);
  // cudaDeviceSynchronize();
  SortKernelDipole<<<blocks * FACTOR4, THREADS4>>>(
    yarn_dipole_tree.nnodes, nbodies, bottom, yarn_dipole_tree.sortl, yarn_dipole_tree.countl, yarn_dipole_tree.startl, yarn_dipole_tree.childl);
  // cudaDeviceSynchronize();
  cudaFree(yarn_dipole_tree.countl);
  cudaFree(yarn_dipole_tree.startl);
  yarn_dipole_tree.countl = nullptr;
  yarn_dipole_tree.startl = nullptr;
}


float ComputeLinkingNumberBHDipole(const YarnDipoleTree& source, int nbodiess,
    const YarnDipoleTree& target, int nbodiest, const int blocks,
    float epssq, float itolsq, float distratiosq) {
  thrust::device_vector<float> linkagest(nbodiest);
  // float epssq = 0.05f * 0.05f;
  // float itolsq = 1.0f / (0.5 * 0.5);
  LinkageCalculationKernelDipole<<<OVERSATURATE_FACTOR5 * blocks * FACTOR5,
                               THREADS5>>>(
      source.nnodes, nbodiess, nbodiest, source.radius, itolsq, epssq, distratiosq,
      source.childl, source.midpts, target.midpts, source.tangents, target.tangents,
      source.dipole1s, source.dipole3s, linkagest.data().get());
  // cudaDeviceSynchronize();
  return thrust::reduce(linkagest.begin(), linkagest.end());
}



float ComputeWritheBHDipole(const YarnDipoleTree& source, int nbodiess, const int blocks,
  float epssq, float itolsq, float distratiosq, bool use_sorted_order) {
thrust::device_vector<float> writhes(nbodiess);
// float epssq = 0.05f * 0.05f;
// float itolsq = 1.0f / (0.5 * 0.5);
WritheCalculationKernelDipole<<<OVERSATURATE_FACTOR5 * blocks * FACTOR5,
                             THREADS5>>>(
    source.nnodes, nbodiess, source.radius, itolsq, epssq, distratiosq, source.sortl,
    source.childl, source.midpts, source.tangents, 
    source.dipole1s, source.dipole3s, writhes.data().get(), use_sorted_order);
// cudaDeviceSynchronize();
return thrust::reduce(writhes.begin(), writhes.end());
}
