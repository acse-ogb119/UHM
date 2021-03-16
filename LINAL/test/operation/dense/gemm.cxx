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
#include "dense_test.hxx"

int main(int argc, char **argv) {

  linal::Flat_    A,  B,  C, D, norm;
  linal::Hier_   hA, hB, hC;

  FLA_Obj *alpha, *beta;


  if (argc != 6) {
    printf("Try :: gemm [thread] [transa] [transb] [alpha] [beta]\n");
    return -1;
  }

  // ---------------------------------------
  // ** Initialization
  FLA_Init();

  int nthread = atoi( (argv[1]) );
  int transa, transb;

  int datatype = TEST_DATATYPE;

  if ( atoi( (argv[2]) ) ) transa = FLA_TRANSPOSE;
  else                     transa = FLA_NO_TRANSPOSE;
  if ( atoi( (argv[3]) ) ) transb = FLA_TRANSPOSE;
  else                     transb = FLA_NO_TRANSPOSE;
  if ( atoi( (argv[4]) ) ) alpha  = &FLA_ONE;
  else                     alpha  = &FLA_MINUS_ONE;
  if ( atoi( (argv[5]) ) ) beta   = &FLA_ONE;
  else                     beta   = &FLA_MINUS_ONE;

  // ---------------------------------------
  // ** MATRICES

  A.create   (datatype, N, N);
  B.create   (datatype, N, N);
  C.create   (datatype, N, N);
  D.create   (datatype, N, N);
  norm.create(datatype, 1,1);

  FLA_Random_matrix(~A);
  FLA_Random_matrix(~B);
  FLA_Random_matrix(~D);

  hA.create(A, BMN, BMN);
  hB.create(B, BMN, BMN);
  hC.create(D, BMN, BMN);

  FLA_Copy( ~D, ~C);

  // ---------------------------------------
  // ** FLAME
  FLA_Gemm( transa, transb, *alpha,
            ~A, ~B, *beta, ~C);

  // ---------------------------------------
  // ** LINAL
  omp_set_num_threads( nthread );
#pragma omp parallel
  {
#pragma omp single nowait
    {
#pragma omp task 
      linal::dense::gemm( transa, transb, *alpha,
                          hA, hB, *beta, hC);
    }
  }

  // ---------------------------------------
  // ** Check
  FLA_Axpy(FLA_MINUS_ONE, ~C, ~D);
  FLA_Norm1( ~D, ~norm);
  
  float  norm_f = (float)norm(0,0);
  double norm_d = (double)norm(0,0);

  int rval;

  printf("- TEST::");
  for (int i=0;i<argc;++i)
    printf(" %s ", argv[i] );
  printf("\n");

  switch (datatype) {
  case LINAL_SINGLE_REAL:
    if (norm_f < LINAL_ERROR_TOL) {
      printf("PASS::Norm :: %E \n", norm_f);   rval = 0;
    } else {
      printf("FAIL::Norm :: %E \n", norm_f);   rval = -1;
    }
    break;
  case LINAL_DOUBLE_REAL:
    if (norm_d < LINAL_ERROR_TOL) {
      printf("PASS::Norm :: %E \n", norm_d);   rval = 0;
    } else {
      printf("FAIL::Norm :: %E \n", norm_d);   rval = -1;
    }
    break;
  }

  // ---------------------------------------
  // ** Matrices
  A.free(); hA.free();
  B.free(); hB.free();
  C.free(); hC.free();
  D.free();

  norm.free();

  // ---------------------------------------
  // ** Finalization
  FLA_Finalize();
  return rval;
}
