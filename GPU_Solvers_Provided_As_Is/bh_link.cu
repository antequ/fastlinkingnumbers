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
#define DISTRATIO2 (1000.f)
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
__device__ volatile TreeMeta tree_meta;

/******************************************************************************/
/*** initialize memory ********************************************************/
/******************************************************************************/

__global__ void InitializationKernel() { blkcntd = 0; }

/******************************************************************************/
/*** compute center and radius ************************************************/
/******************************************************************************/

__global__ __launch_bounds__(THREADS1, FACTOR1) void BoundingBoxKernel(
    const int nnodesd, const int nbodiesd, int *const __restrict__ startd,
    int *const __restrict__ childd, float4 *const __restrict__ midpointsd,
    float3 *const __restrict__ maxd, float3 *const __restrict__ mind) {
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
      tree_meta.radius = fmaxf(val, max.z - min.z) * 0.5f;

      // ADDED TO FIX tree radius to be constant across yarns. 
      // tree_meta.radius = 54.f; // see Alien sweater dims.

      // create root node
      k = nnodesd;
      tree_meta.bottom = k;

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

__global__ __launch_bounds__(1024, 1) void ClearKernel1(
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

__global__ __launch_bounds__(THREADS2, FACTOR2) void TreeBuildingKernel(
    const int nnodesd, const int nbodiesd,
    volatile int *const __restrict__ childd,
    const float4 *const __restrict__ midpointsd) {
  int i, j, depth, skip, inc;
  float x, y, z, r;
  float dx, dy, dz;
  int ch, n, cell, locked, patch;
  float radius;

  // cache root data
  radius = tree_meta.radius * 0.5f;
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

            cell = atomicSub((int *)&(tree_meta.bottom), 1) - 1;
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

__global__ __launch_bounds__(1024, 1) void ClearKernel2(
    const int nnodesd, int *const __restrict__ startd,
    float4 *const __restrict__ midpointsd) {
  int k, inc, bottom;

  bottom = tree_meta.bottom;
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
__global__ __launch_bounds__(THREADS3, FACTOR3) void SummarizationKernel(
    const int nnodesd, const int nbodiesd, const int bottom,
    volatile int *const __restrict__ countd,
    const int *const __restrict__ childd, float4 *const __restrict__ midpointsd,
    float4 *const __restrict__ tangentsd) {
  int i, j, k, ch, inc, cnt;
  float xmin, ymin, zmin, xmax, ymax, zmax;
  float maxw;
  float4 midpoint, tangent, ctangent;
  __shared__ int child[THREADS3 * 8];
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
          countd[k] = cnt;
          //m = 1.0f / ctangent.w;
          midpointsd[k] = make_float4(
              0.5f * (xmin + xmax),
              0.5f * (ymin + ymax),
              0.5f * (zmin + zmax),
              maxw);
          ctangent.w += 0.25f * ((xmax - xmin) * (xmax - xmin)
                                  + (ymax - ymin) * (ymax - ymin)
                                  + (zmax - zmin) * (zmax - zmin));
          tangentsd[k] = ctangent;
          __threadfence();
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
        countd[k] = cnt;
        // m = 1.0f / ctangent.w;
        midpointsd[k] = make_float4(
          0.5f * (xmin + xmax),
          0.5f * (ymin + ymax),
          0.5f * (zmin + zmax),
          maxw);
        ctangent.w += 0.25f * ((xmax - xmin) * (xmax - xmin)
                                + (ymax - ymin) * (ymax - ymin)
                                + (zmax - zmin) * (zmax - zmin));
        tangentsd[k] = ctangent;
        __threadfence();
        k += inc;
      }
    }
  }
}

/******************************************************************************/
/*** sort bodies **************************************************************/
/******************************************************************************/

__global__ __launch_bounds__(THREADS4, FACTOR4) void SortKernel(
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
/*** compute force ************************************************************/
/******************************************************************************/

__global__ __launch_bounds__(THREADS5, FACTOR5) void ForceCalculationKernel(
    const int nnodesd, const int nbodiesd, const float radiusd,
    const float dthfd, const float itolsqd, const float epssqd,
    const int *const __restrict__ sortd, const int *const __restrict__ childd,
    const float4 *const __restrict__ midpointsds,
    const float4 *const __restrict__ midpointsdt,
    const float4 *const __restrict__ tangentsds,
    const float4 *const __restrict__ tangentsdt,
    float4 *const __restrict__ accVeld) {
  int i, j, k, n, depth, base, sbase, diff, pd, nd;
  float ax, ay, az, dx, dy, dz, tmp, distth;
  __shared__ volatile int pos[MAXDEPTH * THREADS5 / WARPSIZE],
      node[MAXDEPTH * THREADS5 / WARPSIZE];
  __shared__ float dq[MAXDEPTH * THREADS5 / WARPSIZE];

  if (0 == threadIdx.x) {
    tmp = radiusd * 2;
    // precompute values that depend only on tree level
    dq[0] = tmp * tmp * itolsqd;
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
  for (k = threadIdx.x + blockIdx.x * blockDim.x; k < nbodiesd;
       k += blockDim.x * gridDim.x) {
    i = sortd[k]; // get permuted/sorted index
    // i = k;
    // cache position info
    const float4 pi = midpointsdt[i];
    // const float4 ti = tangentsdt[i];

    ax = 0.0f;
    ay = 0.0f;
    az = 0.0f;

    // initialize iteration stack, i.e., push root node onto stack
    depth = j;
    if (sbase == threadIdx.x) {
      pos[j] = 0;
      node[j] = nnodesd * 8;
    }

    do {
      // stack is not empty
      pd = pos[depth];
      nd = node[depth];
      while (pd < 8) {
        // node on top of stack has more children to process
        n = childd[nd + pd]; // load child pointer
        pd++;

        if (n >= 0) {
          const float4 pn = midpointsds[n];
          const float4 tn = tangentsds[n];
          dx = pn.x - pi.x;
          dy = pn.y - pi.y;
          dz = pn.z - pi.z;
          tmp =
              dx * dx +
              (dy * dy +
               (dz * dz + epssqd)); // compute distance squared (plus softening)
          // Compute distance threshold
          distth = fmaxf(dq[depth], DISTRATIO2 * (pn.w + pi.w));
          // n < nbodiesd checks if it's a body. the latter checks if cell is
          // far enough away.
          // direct (body) eval
          if (n < nbodiesd) {
            tmp = rsqrtf(tmp); // compute distance with arctan correction
            tmp = tn.w * tmp * tmp * tmp;
            ax += dx * tmp;
            ay += dy * tmp;
            az += dz * tmp;
          } else if (__all_sync(
                         0xffffffff,
                         (tmp >=
                          distth))) { // check if all threads agree that cell
                                      // is far enough away (or is a body)
            tmp = rsqrtf(tmp);        // compute distance
            tmp = tn.w * tmp * tmp * tmp;
            ax += dx * tmp;
            ay += dy * tmp;
            az += dz * tmp;
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
    accVeld[i].x = ax;
    accVeld[i].y = ay;
    accVeld[i].z = az;
  }
}

/******************************************************************************/
/*** compute linkages
 * ************************************************************/
/******************************************************************************/

__global__ __launch_bounds__(THREADS5, FACTOR5) void LinkageCalculationKernel(
    const int nnodesds, const int nbodiesds, const int nbodiesdt, const float radiusd,
    const float itolsqd, const float epssqd,
    const int *const __restrict__ sortdt, const int *const __restrict__ childds,
    const float4 *const __restrict__ midpointsds,
    const float4 *const __restrict__ midpointsdt,
    const float4 *const __restrict__ tangentsds,
    const float4 *const __restrict__ tangentsdt,
    float *const __restrict__ linkagest) {
  int i, j, k, n, depth, base, sbase, diff, pd, nd;
  float lk, dx, dy, dz, tmp, distth;
  float3 alpha, beta, gamma, delta, dr, halftn, halfti;
  float la, lb, lc, ld, div1, div2, ac;
  __shared__ volatile int pos[MAXDEPTH * THREADS5 / WARPSIZE],
      node[MAXDEPTH * THREADS5 / WARPSIZE];
  __shared__ float dq[MAXDEPTH * THREADS5 / WARPSIZE];

  if (0 == threadIdx.x) {
    tmp = radiusd * 2;
    // precompute values that depend only on tree level
    dq[0] = tmp * tmp * itolsqd;
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
    const float4 pi = midpointsdt[i];
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
          const float4 pn = midpointsds[n];
          const float4 tn = tangentsds[n];
          dx = pn.x - pi.x;
          dy = pn.y - pi.y;
          dz = pn.z - pi.z;
          tmp =
              dx * dx +
              (dy * dy +
               (dz * dz + epssqd)); // compute distance squared (plus softening)
          // Compute distance threshold
          distth = fmaxf(dq[depth], DISTRATIO2 * (pn.w + pi.w));
          // distth = fmaxf(8.f * tn.w * itolsqd, DISTRATIO2 * (pn.w + pi.w));
          // n < nbodiesd checks if it's a body. the latter checks if cell is
          // far enough away.
          // direct (body) eval
          if (n < nbodiesds || __all_sync(0xffffffff, (tmp >= distth))) {
            // Store determinant into distth.
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
              tmp = rsqrtf(tmp); // compute distance
              tmp = tmp * tmp * tmp * (INV_FOUR_PI);
              // remember, determinant is in distth.
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
/*** advance bodies ***********************************************************/
/******************************************************************************/

__global__ __launch_bounds__(THREADS6, FACTOR6) void IntegrationKernel(
    const int nbodiesd, const float dtimed, const float dthfd,
    float4 *const __restrict__ midpoints, float2 *const __restrict__ veld,
    float4 *const __restrict__ accVeld) {
  int i, inc;
  float dvelx, dvely, dvelz;
  float velhx, velhy, velhz;

  // iterate over all bodies assigned to thread
  inc = blockDim.x * gridDim.x;
  for (i = threadIdx.x + blockIdx.x * blockDim.x; i < nbodiesd; i += inc) {
    // integrate
    float4 acc = accVeld[i];
    dvelx = acc.x * dthfd;
    dvely = acc.y * dthfd;
    dvelz = acc.z * dthfd;

    float2 v = veld[i];
    velhx = v.x + dvelx;
    velhy = v.y + dvely;
    velhz = acc.w + dvelz;

    float4 p = midpoints[i];
    p.x += velhx * dtimed;
    p.y += velhy * dtimed;
    p.z += velhz * dtimed;
    midpoints[i] = p;

    v.x = velhx + dvelx;
    v.y = velhy + dvely;
    acc.w = velhz + dvelz;
    veld[i] = v;
    accVeld[i] = acc;
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

// random number generator (based on SPLASH-2 code at
// https://github.com/staceyson/splash2/blob/master/codes/apps/barnes/util.C)

static int randx = 7;

static double drnd() {
  const int lastrand = randx;
  randx = (1103515245 * randx + 12345) & 0x7FFFFFFF;
  return (double)lastrand / 2147483648.0;
}

/******************************************************************************/
// performs direct integration on the CPU
static void direct_integrate(int nbodies, float dtime, float dthf,
                             float4 *posMass, float2 *vel, float4 *accVel,
                             float epssqd) {
  int i;

#pragma omp parallel for
  for (i = 0; i < nbodies; ++i) {
    float dvelx, dvely, dvelz;
    float velhx, velhy, velhz;
    float dx, dy, dz;
    float ax, ay, az;
    int j;
    // compute forces
    ax = 0.f;
    ay = 0.f;
    az = 0.f;
    for (j = 0; j < nbodies; ++j) {
      if (j == i)
        continue;
      dx = posMass[j].x - posMass[i].x;
      dy = posMass[j].y - posMass[i].y;
      dz = posMass[j].z - posMass[i].z;
      float tmp = 1.f / sqrtf(dx * dx + dy * dy + dz * dz + epssqd);
      tmp = posMass[j].w * tmp * tmp * tmp;
      ax += dx * tmp;
      ay += dy * tmp;
      az += dz * tmp;
    }
    // integrate
    dvelx = ax * dthf;
    dvely = ay * dthf;
    dvelz = az * dthf;

    float2 v = vel[i];
    velhx = v.x + dvelx;
    velhy = v.y + dvely;
    velhz = accVel[i].w + dvelz;

    float4 p = posMass[i];
    p.x += velhx * dtime;
    p.y += velhy * dtime;
    p.z += velhz * dtime;
    posMass[i] = p;

    v.x = velhx + dvelx;
    v.y = velhy + dvely;
    accVel[i].w = velhz + dvelz;
    vel[i] = v;
    accVel[i].x = ax;
    accVel[i].y = ay;
    accVel[i].z = az;
  }
}

// computes the w coordinates for each yarn
void PrepareYarn(std::vector<float4> &midpts,
                 std::vector<float4> &tangents) {
  if (midpts.size() != tangents.size()) {
    std::cout << "Preparing yarn: midpoint size and tangent size are different."
              << std::endl;
    throw std::runtime_error("PreparingYarn size error.");
  }
  for (int i = 0; i < tangents.size(); ++i) {
    using namespace std;
    midpts[i].w = tangents[i].x * tangents[i].x +
                  tangents[i].y * tangents[i].y + tangents[i].z * tangents[i].z;
    tangents[i].w = 0.25 * midpts[i].w;
  }
}
bool checked_dev_properties = false;
void CheckDevicePropertiesForDev(const cudaDeviceProp& deviceProp) {
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
  checked_dev_properties = true;
}
// returns a <midpts_dev, tangents_dev, sortl, childl, radius, nnodes> for the
// tree
YarnTree PrepareTree(const std::vector<float4> &midpts,
                     const std::vector<float4> &tangents) {
  int nnodes, nbodies;
  int blocks;
  float4 *midpts_dev;
  float4 *tangents_dev;
  int *sortl;
  int *childl;
  int *countl;
  int *startl;
  float3 *minl, *maxl;
  TreeMeta cpu_tree_meta;
  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, 0 /* device number */);
  if (!checked_dev_properties) CheckDevicePropertiesForDev(deviceProp);

  blocks = deviceProp.multiProcessorCount;
  
  cudaFuncSetCacheConfig(BoundingBoxKernel, cudaFuncCachePreferShared);
  cudaFuncSetCacheConfig(TreeBuildingKernel, cudaFuncCachePreferL1);
  cudaFuncSetCacheConfig(ClearKernel1, cudaFuncCachePreferL1);
  cudaFuncSetCacheConfig(ClearKernel2, cudaFuncCachePreferL1);
  cudaFuncSetCacheConfig(SummarizationKernel, cudaFuncCachePreferShared);
  cudaFuncSetCacheConfig(SortKernel, cudaFuncCachePreferL1);
  cudaGetLastError(); // reset error value
  nbodies = midpts.size();
  nnodes = nbodies * 2.5; // experimentally, 2.1 is enough.
  if (nnodes < 1024 * blocks)
    nnodes = 1024 * blocks;
  while ((nnodes & (WARPSIZE - 1)) != 0)
    nnodes++;
  nnodes--;
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
      cudaMalloc((void **)&midpts_dev, sizeof(float4) * (nnodes + 1)))
    fprintf(stderr, "could not allocate midptsd\n");
  CudaTest("couldn't allocate midptsd");
  if (cudaSuccess !=
      cudaMalloc((void **)&tangents_dev, sizeof(float4) * (nnodes + 1)))
    fprintf(stderr, "could not allocate tangentsd\n");
  CudaTest("couldn't allocate tangentsd");
  if (cudaSuccess !=
      cudaMalloc((void **)&maxl, sizeof(float3) * blocks * FACTOR1))
    fprintf(stderr, "could not allocate maxd\n");
  CudaTest("couldn't allocate maxd");
  if (cudaSuccess !=
      cudaMalloc((void **)&minl, sizeof(float3) * blocks * FACTOR1))
    fprintf(stderr, "could not allocate mind\n");
  CudaTest("couldn't allocate mind");
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
  // Launch kernels.
  InitializationKernel<<<1, 1>>>();
  // Kernel 1
  BoundingBoxKernel<<<blocks * FACTOR1, THREADS1>>>(
      nnodes, nbodies, startl, childl, midpts_dev, maxl, minl);
  cudaDeviceSynchronize();
  // Kernel 2
  ClearKernel1<<<blocks * 1, 1024>>>(nnodes, nbodies, childl);
  TreeBuildingKernel<<<blocks * FACTOR2, THREADS2>>>(nnodes, nbodies, childl,
                                                     midpts_dev);
  ClearKernel2<<<blocks * 1, 1024>>>(nnodes, startl, midpts_dev);
  cudaDeviceSynchronize();
  // download data from Kernel 2
  if (cudaSuccess != cudaMemcpyFromSymbol(&cpu_tree_meta, tree_meta,
                                          sizeof(TreeMeta), 0,
                                          cudaMemcpyDeviceToHost))
    fprintf(stderr, "copying of tree_meta to device failed\n");
  CudaTest("tree_meta copy to device failed");
  // Kernel 3
  SummarizationKernel<<<blocks * FACTOR3, THREADS3>>>(
      nnodes, nbodies, cpu_tree_meta.bottom, countl, childl, midpts_dev,
      tangents_dev);
  cudaDeviceSynchronize();
  SortKernel<<<blocks * FACTOR4, THREADS4>>>(
      nnodes, nbodies, cpu_tree_meta.bottom, sortl, countl, startl, childl);
  cudaDeviceSynchronize();

  cudaFree(countl);
  cudaFree(startl);
  cudaFree(minl);
  cudaFree(maxl);
  return {midpts_dev, tangents_dev,         sortl,
          childl,     cpu_tree_meta.radius, nnodes};
}

float ComputeLinkingNumberBH(const YarnTree& source, int nbodiess, const YarnTree& target, int nbodiest, float epssq, float itolsq) {
  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, 0 /* device number */);
  if (!checked_dev_properties) CheckDevicePropertiesForDev(deviceProp);
  int blocks = deviceProp.multiProcessorCount;
  thrust::device_vector<float> linkagest(nbodiest);
  // float epssq = 0.05f * 0.05f;
  // float itolsq = 1.0f / (0.5 * 0.5);
  cudaFuncSetCacheConfig(LinkageCalculationKernel, cudaFuncCachePreferEqual);
  LinkageCalculationKernel<<<OVERSATURATE_FACTOR5 * blocks * FACTOR5,
                               THREADS5>>>(
      source.nnodes, nbodiess, nbodiest, source.radius, itolsq, epssq, target.sortl,
      source.childl, source.midpts, target.midpts, source.tangents, target.tangents,
      linkagest.data().get());
  cudaDeviceSynchronize();
  return thrust::reduce(linkagest.begin(), linkagest.end());
}

int main_old(int argc, char *argv[]) {
  int i, run, blocks;
  int nnodes, nbodies, step, timesteps;
  double runtime;
  float dtime, dthf, epssq, itolsq;
  float time, timing[8];
  cudaEvent_t start, stop;

  // accVel's x,y,z are acceleration,
  // vel is x and y velocity; accVel's w is z velocity.
  float4 *accVel;
  float2 *vel;
  int *sortl, *childl, *countl, *startl;
  float4 *accVell;
  float2 *vell;
  float3 *maxl, *minl;
  // x,y,z are position, w is mass. w is the fourth value.
  float4 *posMassl;
  float4 *tangentsl;
  float4 *posMass;
  double rsc, vsc, r, v, x, y, z, sq, scale;

  float4 *accVelc;
  float2 *velc;
  float4 *posMassc;

  // perform some checks

  printf("Based on ECL-BH v4.5\n");
  printf("Copyright (c) 2010-2020 Texas State University\n");
  fflush(stdout);

  if (argc != 5) {
    fprintf(stderr, "\n");
    fprintf(
        stderr,
        "arguments: number_of_bodies number_of_timesteps device skip_cpu\n");
    exit(-1);
  }

  int deviceCount;
  cudaGetDeviceCount(&deviceCount);
  if (deviceCount == 0) {
    fprintf(stderr, "There is no device supporting CUDA\n");
    exit(-1);
  }

  const int dev = atoi(argv[3]);
  if ((dev < 0) || (deviceCount <= dev)) {
    fprintf(stderr, "There is no device %d\n", dev);
    exit(-1);
  }
  cudaSetDevice(dev);

  bool skip_cpu = (atoi(argv[4]) != 0);
  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, dev);
  if ((deviceProp.major == 9999) && (deviceProp.minor == 9999)) {
    fprintf(stderr, "There is no CUDA capable device\n");
    exit(-1);
  }
  if (deviceProp.major < 3) {
    fprintf(stderr, "Need at least compute capability 3.0\n");
    exit(-1);
  }
  if (deviceProp.warpSize != WARPSIZE) {
    fprintf(stderr, "Warp size must be %d\n", deviceProp.warpSize);
    exit(-1);
  }

  blocks = deviceProp.multiProcessorCount;
  const int mTSM = deviceProp.maxThreadsPerMultiProcessor;
  printf("gpu: %s with %d SMs and %d mTpSM (%.1f MHz and %.1f MHz)\n",
         deviceProp.name, blocks, mTSM, deviceProp.clockRate * 0.001,
         deviceProp.memoryClockRate * 0.001);

  if ((WARPSIZE <= 0) || (WARPSIZE & (WARPSIZE - 1) != 0)) {
    fprintf(stderr, "Warp size must be greater than zero and a power of two\n");
    exit(-1);
  }
  if (MAXDEPTH > WARPSIZE) {
    fprintf(stderr, "MAXDEPTH must be less than or equal to WARPSIZE\n");
    exit(-1);
  }
  if ((THREADS1 <= 0) || (THREADS1 & (THREADS1 - 1) != 0)) {
    fprintf(stderr, "THREADS1 must be greater than zero and a power of two\n");
    exit(-1);
  }

  // set L1/shared memory configuration
  cudaFuncSetCacheConfig(BoundingBoxKernel, cudaFuncCachePreferShared);
  cudaFuncSetCacheConfig(TreeBuildingKernel, cudaFuncCachePreferL1);
  cudaFuncSetCacheConfig(ClearKernel1, cudaFuncCachePreferL1);
  cudaFuncSetCacheConfig(ClearKernel2, cudaFuncCachePreferL1);
  cudaFuncSetCacheConfig(SummarizationKernel, cudaFuncCachePreferShared);
  cudaFuncSetCacheConfig(SortKernel, cudaFuncCachePreferL1);
  cudaFuncSetCacheConfig(ForceCalculationKernel, cudaFuncCachePreferEqual);
  cudaFuncSetCacheConfig(IntegrationKernel, cudaFuncCachePreferL1);

  cudaGetLastError(); // reset error value
  for (run = 0; run < 1;
       run++) { // in case multiple runs are desired for timing
    for (i = 0; i < 8; i++)
      timing[i] = 0.0f;

    nbodies = atoi(argv[1]);
    if (nbodies < 1) {
      fprintf(stderr, "nbodies is too small: %d\n", nbodies);
      exit(-1);
    }
    if (nbodies > (1 << 30)) {
      fprintf(stderr, "nbodies is too large: %d\n", nbodies);
      exit(-1);
    }
    nnodes = nbodies * 2;
    if (nnodes < 1024 * blocks)
      nnodes = 1024 * blocks;
    while ((nnodes & (WARPSIZE - 1)) != 0)
      nnodes++;
    nnodes--;

    timesteps = atoi(argv[2]);
    dtime = 0.025;
    dthf = dtime * 0.5f;
    epssq = 0.05 * 0.05;
    itolsq = 1.0f / (0.5 * 0.5);

    // allocate memory

    if (run == 0) {
      printf("configuration: %d bodies, %d time steps\n", nbodies, timesteps);

      accVel = (float4 *)malloc(sizeof(float4) * nbodies);
      if (accVel == NULL) {
        fprintf(stderr, "cannot allocate accVel\n");
        exit(-1);
      }
      vel = (float2 *)malloc(sizeof(float2) * nbodies);
      if (vel == NULL) {
        fprintf(stderr, "cannot allocate vel\n");
        exit(-1);
      }
      posMass = (float4 *)malloc(sizeof(float4) * nbodies);
      if (posMass == NULL) {
        fprintf(stderr, "cannot allocate posMass\n");
        exit(-1);
      }

      accVelc = (float4 *)malloc(sizeof(float4) * nbodies);
      if (accVelc == NULL) {
        fprintf(stderr, "cannot allocate accVelc\n");
        exit(-1);
      }
      velc = (float2 *)malloc(sizeof(float2) * nbodies);
      if (velc == NULL) {
        fprintf(stderr, "cannot allocate velc\n");
        exit(-1);
      }
      posMassc = (float4 *)malloc(sizeof(float4) * nbodies);
      if (posMassc == NULL) {
        fprintf(stderr, "cannot allocate posMassc\n");
        exit(-1);
      }

      if (cudaSuccess !=
          cudaMalloc((void **)&childl, sizeof(int) * (nnodes + 1) * 8))
        fprintf(stderr, "could not allocate childd\n");
      CudaTest("couldn't allocate childd");
      if (cudaSuccess !=
          cudaMalloc((void **)&vell, sizeof(float2) * (nnodes + 1)))
        fprintf(stderr, "could not allocate veld\n");
      CudaTest("couldn't allocate veld");
      if (cudaSuccess !=
          cudaMalloc((void **)&accVell, sizeof(float4) * (nnodes + 1)))
        fprintf(stderr, "could not allocate accVeld\n");
      CudaTest("couldn't allocate accVeld");
      if (cudaSuccess !=
          cudaMalloc((void **)&countl, sizeof(int) * (nnodes + 1)))
        fprintf(stderr, "could not allocate countd\n");
      CudaTest("couldn't allocate countd");
      if (cudaSuccess !=
          cudaMalloc((void **)&startl, sizeof(int) * (nnodes + 1)))
        fprintf(stderr, "could not allocate startd\n");
      CudaTest("couldn't allocate startd");
      if (cudaSuccess !=
          cudaMalloc((void **)&sortl, sizeof(int) * (nnodes + 1)))
        fprintf(stderr, "could not allocate sortd\n");
      CudaTest("couldn't allocate sortd");

      if (cudaSuccess !=
          cudaMalloc((void **)&posMassl, sizeof(float4) * (nnodes + 1)))
        fprintf(stderr, "could not allocate posMassd\n");
      CudaTest("couldn't allocate posMassd");
      if (cudaSuccess !=
          cudaMalloc((void **)&tangentsl, sizeof(float4) * (nnodes + 1)))
        fprintf(stderr, "could not allocate tangentsd\n");
      CudaTest("couldn't allocate tangentsd");

      if (cudaSuccess !=
          cudaMalloc((void **)&maxl, sizeof(float3) * blocks * FACTOR1))
        fprintf(stderr, "could not allocate maxd\n");
      CudaTest("couldn't allocate maxd");
      if (cudaSuccess !=
          cudaMalloc((void **)&minl, sizeof(float3) * blocks * FACTOR1))
        fprintf(stderr, "could not allocate mind\n");
      CudaTest("couldn't allocate mind");
    }

    // generate input (based on SPLASH-2 code at
    // https://github.com/staceyson/splash2/blob/master/codes/apps/barnes/code.C)

    rsc = (3 * 3.1415926535897932384626433832795) / 16;
    vsc = sqrt(1.0 / rsc);
    for (i = 0; i < nbodies; i++) {
      float4 p;
      p.w = 1.0 / nbodies;
      r = 1.0 / sqrt(pow(drnd() * 0.999, -2.0 / 3.0) - 1);
      do {
        x = drnd() * 2.0 - 1.0;
        y = drnd() * 2.0 - 1.0;
        z = drnd() * 2.0 - 1.0;
        sq = x * x + y * y + z * z;
      } while (sq > 1.0);
      scale = rsc * r / sqrt(sq);
      p.x = x * scale;
      p.y = y * scale;
      p.z = z * scale;
      posMass[i] = p;

      do {
        x = drnd();
        y = drnd() * 0.1;
      } while (y > x * x * pow(1 - x * x, 3.5));
      v = x * sqrt(2.0 / sqrt(1 + r * r));
      do {
        x = drnd() * 2.0 - 1.0;
        y = drnd() * 2.0 - 1.0;
        z = drnd() * 2.0 - 1.0;
        sq = x * x + y * y + z * z;
      } while (sq > 1.0);
      scale = vsc * v / sqrt(sq);
      float2 v;
      v.x = x * scale;
      v.y = y * scale;
      accVel[i].w = z * scale;
      vel[i] = v;
    }

    if (cudaSuccess != cudaMemcpy(accVell, accVel, sizeof(float4) * nbodies,
                                  cudaMemcpyHostToDevice))
      fprintf(stderr, "copying of vel to device failed\n");
    CudaTest("vel copy to device failed");
    if (cudaSuccess !=
        cudaMemcpy(vell, vel, sizeof(float2) * nbodies, cudaMemcpyHostToDevice))
      fprintf(stderr, "copying of vel to device failed\n");
    CudaTest("vel copy to device failed");
    if (cudaSuccess != cudaMemcpy(posMassl, posMass, sizeof(float4) * nbodies,
                                  cudaMemcpyHostToDevice))
      fprintf(stderr, "copying of posMass to device failed\n");
    CudaTest("posMass copy to device failed");
    if (cudaSuccess != cudaMemcpy(tangentsl, posMass, sizeof(float4) * nbodies,
                                  cudaMemcpyHostToDevice))
      fprintf(stderr, "copying of tangent to device failed\n");
    CudaTest("tangent copy to device failed");

    memcpy(accVelc, accVel, sizeof(float4) * nbodies);
    memcpy(velc, vel, sizeof(float2) * nbodies);
    memcpy(posMassc, posMass, sizeof(float4) * nbodies);
    // run timesteps (launch GPU kernels)

    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    struct timeval starttime, endtime;
    gettimeofday(&starttime, NULL);

    cudaEventRecord(start, 0);
    InitializationKernel<<<1, 1>>>();
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&time, start, stop);
    timing[0] += time;
    // CudaTest("kernel 0 launch failed");

    for (step = 0; step < timesteps; step++) {
      cudaEventRecord(start, 0);
      BoundingBoxKernel<<<blocks * FACTOR1, THREADS1>>>(
          nnodes, nbodies, startl, childl, posMassl, maxl, minl);
      cudaEventRecord(stop, 0);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&time, start, stop);
      timing[1] += time;
      // CudaTest("kernel 1 launch failed");

      cudaEventRecord(start, 0);
      ClearKernel1<<<blocks * 1, 1024>>>(nnodes, nbodies, childl);
      TreeBuildingKernel<<<blocks * FACTOR2, THREADS2>>>(nnodes, nbodies,
                                                         childl, posMassl);
      ClearKernel2<<<blocks * 1, 1024>>>(nnodes, startl, posMassl);
      cudaEventRecord(stop, 0);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&time, start, stop);
      timing[2] += time;
      // CudaTest("kernel 2 launch failed");
      TreeMeta cpu_tree_meta;
      if (cudaSuccess != cudaMemcpyFromSymbol(&cpu_tree_meta, tree_meta,
                                              sizeof(TreeMeta), 0,
                                              cudaMemcpyDeviceToHost))
        fprintf(stderr, "copying of tree_meta to device failed\n");
      CudaTest("tree_meta copy to device failed");
      // if (cudaSuccess !=cudaMemcpyFromSymbol(&radius, radiusd, sizeof(float),
      // 0, cudaMemcpyDeviceToHost)) fprintf(stderr, "copying of radius to
      // device failed\n");  CudaTest("radius copy to device failed");

      cudaEventRecord(start, 0);
      SummarizationKernel<<<blocks * FACTOR3, THREADS3>>>(
          nnodes, nbodies, cpu_tree_meta.bottom, countl, childl, posMassl,
          tangentsl);
      cudaEventRecord(stop, 0);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&time, start, stop);
      timing[3] += time;
      // CudaTest("kernel 3 launch failed");

      cudaEventRecord(start, 0);
      SortKernel<<<blocks * FACTOR4, THREADS4>>>(
          nnodes, nbodies, cpu_tree_meta.bottom, sortl, countl, startl, childl);
      cudaEventRecord(stop, 0);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&time, start, stop);
      timing[4] += time;
      // CudaTest("kernel 4 launch failed");

      cudaEventRecord(start, 0);
      ForceCalculationKernel<<<OVERSATURATE_FACTOR5 * blocks * FACTOR5,
                               THREADS5>>>(
          nnodes, nbodies, cpu_tree_meta.radius, dthf, itolsq, epssq, sortl,
          childl, posMassl, posMassl, tangentsl, tangentsl, accVell);
      cudaEventRecord(stop, 0);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&time, start, stop);
      timing[5] += time;
      // CudaTest("kernel 5 launch failed");

      cudaEventRecord(start, 0);
      IntegrationKernel<<<blocks * FACTOR6, THREADS6>>>(
          nbodies, dtime, dthf, posMassl, vell, accVell);
      cudaEventRecord(stop, 0);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&time, start, stop);
      timing[6] += time;
      // CudaTest("kernel 6 launch failed");
    }
    CudaTest("kernel launch failed");
    printf("Finished GPU. Testing CPU now.\n");
    cudaEventRecord(start, 0);
    if (!skip_cpu) {
      for (step = 0; step < timesteps; ++step) {
        direct_integrate(nbodies, dtime, dthf, posMassc, velc, accVelc, epssq);
      }
    }
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&time, start, stop);
    timing[7] += time;
    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    gettimeofday(&endtime, NULL);
    runtime = (endtime.tv_sec + endtime.tv_usec / 1000000.0 - starttime.tv_sec -
               starttime.tv_usec / 1000000.0);

    printf("runtime: %.4lf s  (", runtime);
    time = 0;
    for (i = 1; i < 8; i++) {
      printf(" %.1f ", timing[i]);
      time += timing[i];
    }
    printf(") = %.1f ms\n", time);
  }

  // transfer final result back to CPU
  if (cudaSuccess != cudaMemcpy(accVel, accVell, sizeof(float4) * nbodies,
                                cudaMemcpyDeviceToHost))
    fprintf(stderr, "copying of accVel from device failed\n");
  CudaTest("vel copy from device failed");
  if (cudaSuccess !=
      cudaMemcpy(vel, vell, sizeof(float2) * nbodies, cudaMemcpyDeviceToHost))
    fprintf(stderr, "copying of vel from device failed\n");
  CudaTest("vel copy from device failed");
  if (cudaSuccess != cudaMemcpy(posMass, posMassl, sizeof(float4) * nbodies,
                                cudaMemcpyDeviceToHost))
    fprintf(stderr, "copying of posMass from device failed\n");
  CudaTest("posMass copy from device failed");

  // print output
  i = 0;
  int nprint = 10;
  // for (i = 0; i < nbodies; i++) {
  float diffsq = 0.f;
  float cpumagnsq = 0.f;
  for (i = 0; i < nprint; ++i) {
    printf("GPU body %d: %.2e %.2e %.2e\n", i, posMass[i].x, posMass[i].y,
           posMass[i].z);
    if (!skip_cpu) {
      printf("CPU body %d: %.2e %.2e %.2e\n", i, posMassc[i].x, posMassc[i].y,
             posMassc[i].z);
    }
  }
  if (!skip_cpu) {
#pragma omp parallel for reduction(+ : diffsq, cpumagnsq)
    for (i = 0; i < nbodies; ++i) {
      float dx = posMassc[i].x - posMass[i].x;
      float dy = posMassc[i].y - posMass[i].y;
      float dz = posMassc[i].z - posMass[i].z;
      cpumagnsq += posMassc[i].x * posMassc[i].x +
                   posMassc[i].y * posMassc[i].y +
                   posMassc[i].z * posMassc[i].z;
      diffsq += dx * dx + dy * dy + dz * dz;
    }
    printf("Relative L2 error: %.2e\n", sqrtf(diffsq / cpumagnsq));
  }

  free(vel);
  free(accVel);
  free(posMass);
  free(velc);
  free(accVelc);
  free(posMassc);

  cudaFree(childl);
  cudaFree(vell);
  cudaFree(accVell);
  cudaFree(countl);
  cudaFree(startl);
  cudaFree(sortl);
  cudaFree(posMassl);
  cudaFree(tangentsl);
  cudaFree(maxl);
  cudaFree(minl);

  return 0;
}
