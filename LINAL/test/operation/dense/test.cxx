/*
  Copyright © 2011, Kyungjoo Kim
  All rights reserved.
  
  This file is part of LINAL.
  
  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:

  1. Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.

  2. Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.

  3. Neither the name of the owner nor the names of its contributors
    may be used to endorse or promote products derived from this software
    without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  POSSIBILITY OF SUCH DAMAGE.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cublas.h>

#define M 6
#define N 5

#define IDX2C(i,j,ld) (((j)*(ld))+(i))

void modify (double *m, int ldm, int n, int p, int q, double alpha,double beta)
{
  cublasDscal (n-p, alpha, &m[IDX2C(p,q,ldm)], ldm);
  cublasDscal (ldm-p, beta, &m[IDX2C(p,q,ldm)], 1);
}

int main(int argc, char *argv[])
{
  int i, j;
  cublasStatus stat;
  double* devPtrA;
  double* a = 0;
  a = (double *)malloc (M * N * sizeof (*a));
  if (!a) {
    printf ("host memory allocation failed");
    return 1;
  }
  for (j = 0; j < N; j++) {
    for (i = 0; i < M; i++) {
      a[IDX2C(i,j,M)] = i * M + j + 1;
    }
  }

  cublasInit();
  stat = cublasAlloc (M*N, sizeof(*a), (void**)&devPtrA);
  if (stat != CUBLAS_STATUS_SUCCESS) {
    printf ("device memory allocation failed");
    return 1;
  }
  cublasSetMatrix (M, N, sizeof(*a), a, M, devPtrA, M);
  modify (devPtrA, M, N, 1, 2, 16.0f, 12.0f);

  for (j = 0; j < N; j++) {
    for (i = 0; i < M; i++) {
      printf ("%8.2f", a[IDX2C(i,j,M)]);
    }
    printf ("\n");
  }
  printf ("\n");

  cublasGetMatrix (M, N, sizeof(*a), devPtrA, M, a, M);


  for (j = 0; j < N; j++) {
    for (i = 0; i < M; i++) {
      printf ("%8.2f", a[IDX2C(i,j,M)]);
    }
    printf ("\n");
  }
  printf ("\n");

  cublasFree (devPtrA);
  cublasShutdown();
  for (j = 0; j < N; j++) {
    for (i = 0; i < M; i++) {
      printf ("%8.2f", a[IDX2C(i,j,M)]);
    }
    printf ("\n");
  }

  return 0;
}
