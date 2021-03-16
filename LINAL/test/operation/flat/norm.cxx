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
#include "linal.hxx"

#define LINAL_ERROR_TOL 1.0e-9

int main(int argc, char **argv) {


  int uplo, diag;
  double sample, norm, inv_norm_est, inv_norm, diff;
  linal::Flat_ A, B, C, p;
  
  uplo = FLA_LOWER_TRIANGULAR;
  diag = FLA_UNIT_DIAG;

  if (argc != 2) {
    printf("Try :: norm [N]\n");
    return -1;
  }

  int N = atoi( (argv[1]) );
  int datatype = LINAL_DOUBLE_REAL;

  // ---------------------------------------
  // ** Initialization  
  FLA_Init();


  // ---------------------------------------
  // ** Matrices

  // test sample
  A.create    (datatype, N, N);
  p.create    (LINAL_INT,  N, 1);

  // inverse matrix
  B.create    (datatype, N, N);
  C.create    (datatype, N, N);

  // random setup and factorization
  FLA_Random_matrix(~A);
  FLA_LU_piv(~A,~p);
  FLA_LU_nopiv(~A);

  // get inverse
  FLA_Copy(~A, ~B);
  FLA_Triangularize(uplo, diag, ~B);
  FLA_Obj_set_to_scalar(FLA_ZERO, ~C);
  for (int i=0;i<N;++i)
    C(i,i)= 1.0;

  FLA_Trsm( FLA_LEFT, uplo, FLA_NO_TRANSPOSE,
	    diag, FLA_ONE, ~B, ~C );

  // ---------------------------------------
  // ** Check

  // calculate norm
  norm         = linal::norm1(uplo, diag, A);
  inv_norm_est = linal::inv_norm1(uplo, diag, A);
  inv_norm     = linal::norm1(uplo, diag, C);
  diff         = inv_norm_est - inv_norm;
  
  int rval;
  printf("- TEST::");
  for (int i=0;i<argc;++i)
    printf(" %s ", argv[i] );
  printf("\n");

  if (diff < LINAL_ERROR_TOL) {
    printf("PASS::Norm %E, Inv_Norm_Est %E, Inv_norm %E, Diff %E ( %lf percents error)\n",
	   norm, inv_norm_est, inv_norm, diff, -1.0*diff/inv_norm*100);
    rval = 0;
  } else {
    printf("FAIL::Norm %E, Inv_Norm_Est %E, Inv_norm %E, Diff %E ( %lf percents error)\n",
	   norm, inv_norm_est, inv_norm, diff, -1.0*diff/inv_norm*100);
    rval = -1;
  }

  // ---------------------------------------
  // ** Matrices
  C.free();
  B.free();

  p.free();
  A.free();

  // ---------------------------------------
  // ** Finalization
  FLA_Finalize();

  return -1;
}
