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
#include "gpu_test.hxx"

int main(int argc, char **argv) {


  linal::Flat_  A,  B, C, norm;
  linal::Hier_ hA, hB;

  FLA_Obj *alpha;

  if (argc != 7) {
    printf("Try :: trsm [thread] [side] [uplo] [trans] [diag] [alpha]\n");
    return -1;
  }

  int side, uplo, trans, diag;
  int nthread = atoi( (argv[1]) ), ngpu = 1, n_flat_gpu = 3;

  int datatype = TEST_DATATYPE;

  // ---------------------------------------
  // ** Initialization  
  FLA_Init();
	
  if ( atoi( (argv[2]) ) ) side  = FLA_LEFT;
  else                     side  = FLA_RIGHT;
  if ( atoi( (argv[3]) ) ) uplo  = FLA_LOWER_TRIANGULAR;
  else                     uplo  = FLA_UPPER_TRIANGULAR;
  if ( atoi( (argv[4]) ) ) trans = FLA_TRANSPOSE;
  else                     trans = FLA_NO_TRANSPOSE;
  if ( atoi( (argv[5]) ) ) diag  = FLA_UNIT_DIAG;
  else                     diag  = FLA_NONUNIT_DIAG;
  if ( atoi( (argv[6]) ) ) alpha = &FLA_ONE;
  else                     alpha = &FLA_MINUS_ONE;

  // ---------------------------------------
  // ** Matrices
  A.create   (datatype, N, N);
  B.create   (datatype, N, N);
  C.create   (datatype, N, N);
  norm.create(datatype, 1, 1);

  FLA_Random_spd_matrix(uplo, ~A);
  FLA_Random_matrix(~B);

  hA.create(A, BMN, BMN);
  hB.create(B, BMN, BMN);

  FLA_Copy(~B, ~C);
	
  // ---------------------------------------
  // ** FLAME
  FLA_Trsm( side, uplo, trans, diag, *alpha, ~A, ~C);

  // ---------------------------------------
  // ** GPU
  omp_set_num_threads( nthread );
  linal::init_gpu(nthread, ngpu, n_flat_gpu);
  linal::set_computing_model( LINAL_GPU );
  linal::create_flat_gpu( datatype, BMN, BMN );

#pragma omp parallel
  {
#pragma omp single nowait
    linal::dense::trsm( side, uplo, trans, diag, *alpha, hA, hB );
  }

  linal::set_computing_model( LINAL_CPU );
  linal::free_flat_gpu();
  linal::finalize_gpu();
	
  // ---------------------------------------
  // ** Check
  FLA_Axpy(FLA_MINUS_ONE, ~B, ~C);
  FLA_Norm1( ~C, ~norm);

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
  C.free();
  
  // ---------------------------------------
  // ** Finalization
  norm.free();
  FLA_Finalize();

  return rval;
}
