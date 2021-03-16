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


  linal::Flat_  A,  B, norm;
  linal::Hier_ hA;

  if (argc != 3) {
    printf("Try :: trsm [thread] [uplo]\n");
    return -1;
  }

  int uplo;
  int nthread = atoi( (argv[1]) );

  int datatype = TEST_DATATYPE;

  // ---------------------------------------
  // ** Initialization   
  FLA_Init();

  if ( atoi( (argv[2]) ) ) uplo  = FLA_LOWER_TRIANGULAR;
  else                     uplo  = FLA_UPPER_TRIANGULAR;

  // ---------------------------------------
  // ** Matrices
  A.create   (datatype, N, N);
  B.create   (datatype, N, N);
  norm.create(datatype, 1, 1);

  FLA_Random_spd_matrix(uplo, ~A);

  hA.create(A, BMN, BMN);

  FLA_Copy( ~A, ~B );

  // ---------------------------------------
  // ** FLAME
  FLA_Chol( uplo, ~B );

  // ---------------------------------------
  // ** LINAL
  omp_set_num_threads( nthread );
#pragma omp parallel
  {
#pragma omp single nowait
    linal::dense::chol( uplo, hA );
  }

  // ---------------------------------------
  // ** Check
  FLA_Triangularize(uplo, FLA_NONUNIT_DIAG, ~A);
  FLA_Triangularize(uplo, FLA_NONUNIT_DIAG, ~B);

  FLA_Axpy(FLA_MINUS_ONE, ~A, ~B);
  FLA_Norm1( ~B, ~norm);

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
  B.free();
  
  norm.free();

  // ---------------------------------------
  // ** Finalization
  FLA_Finalize();

  return rval;
}
