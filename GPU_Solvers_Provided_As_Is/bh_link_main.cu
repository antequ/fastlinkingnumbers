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
// -O3" -use_fast_math bh_link_main.cu -o bhlink -I../common
#include "bh_link.cu"

int main(int argc, char *argv[]) {
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

  printf("ECL-BH v4.5\n");
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
