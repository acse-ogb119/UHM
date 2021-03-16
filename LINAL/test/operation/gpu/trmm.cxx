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

  linal::Flat_    A1,  B1,  A2,  B2, norm;
  linal::Hier_   hA2, hB2;

  if (argc != 5) {
    printf("Try :: gemm [thread] [side] [uplo] [trans] \n");
    return -1;
  }

  int side, uplo, trans;
  int nthread = atoi( (argv[1]) ), ngpu = 1, n_flat_gpu = 3;

  int datatype = TEST_DATATYPE;
	
  // ---------------------------------------
  // ** Initialization   
  FLA_Init();
	
  if ( atoi( (argv[2]) ) ) side = FLA_LEFT;
  else                     side = FLA_RIGHT;
  if ( atoi( (argv[3]) ) ) uplo = FLA_LOWER_TRIANGULAR;
  else                     uplo = FLA_UPPER_TRIANGULAR;
  if ( atoi( (argv[4]) ) ) trans = FLA_NO_TRANSPOSE;
  else                     trans = FLA_TRANSPOSE;

  // ---------------------------------------
  // ** Matrices
  A1.create   (datatype, N, N);
  B1.create   (datatype, N, N);
  A2.create   (datatype, N, N);
  B2.create   (datatype, N, N);
  
  norm.create(datatype, 1,1);

  FLA_Random_matrix(~A1);
  FLA_Random_matrix(~B1);
  FLA_Copy( ~A1, ~A2 );
  FLA_Copy( ~B1, ~B2 );

  hA2.create(A2, BMN, BMN);
  hB2.create(B2, BMN, BMN);

  // ---------------------------------------
  // ** FLAME
  FLA_Trmm( side, uplo, trans, 
            FLA_NONUNIT_DIAG, FLA_ONE, 
            ~A1, ~B1 );

  // ---------------------------------------
  // ** GPU
  omp_set_num_threads( nthread );
  linal::init_gpu(nthread, ngpu, n_flat_gpu);
  linal::set_computing_model( LINAL_GPU );
  linal::create_flat_gpu( datatype, BMN, BMN );  

#pragma omp parallel
  {
#pragma omp single nowait
    linal::dense::trmm( side, uplo, trans,
                        FLA_NONUNIT_DIAG, FLA_ONE,
                        hA2, hB2 );
  }

  linal::set_computing_model( LINAL_CPU );
  linal::free_flat_gpu();
  linal::finalize_gpu();
	
  // ---------------------------------------
  // ** Check	
  FLA_Axpy(FLA_MINUS_ONE, ~B1, ~B2);
  FLA_Norm1( ~B2, ~norm);
  
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
  A1.free(); A2.free(); hA2.free();
  B1.free(); B2.free(); hB2.free();

  norm.free();
  // ---------------------------------------
  // ** Finalization
  FLA_Finalize();

  return rval;
}
