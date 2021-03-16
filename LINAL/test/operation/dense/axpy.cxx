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

#define N 1000
#define LINAL_ERROR_TOL 1.0e-9

int main(int argc, char **argv) {
  FLA_Init();
  linal::Flat_    A,  B;
  linal::Hier_   hA, hB;

  double t, t_block, t_unblock;

  if (argc != 3) {
    printf("Try :: axpy [thread][BMN]\n");
    return 0;
  }

  int m, n;
  int nthread  = atoi( (argv[1]) );
  int BMN      = atoi( (argv[2]) );
  int datatype = LINAL_DOUBLE_REAL;

  A.create   (datatype, N, N);
  B.create   (datatype, N, N);

  FLA_Random_matrix(~A);
  FLA_Random_matrix(~B);

  hA.create(A, BMN, BMN);
  hB.create(B, BMN, BMN);

  t = FLA_Clock();
  FLA_Axpy( FLA_ONE, ~A, ~B );
  t_unblock =  FLA_Clock() - t;
  
  m = hA.get_m();
  n = hA.get_n();

  omp_set_num_threads( nthread );

  t =  FLA_Clock();
#pragma omp parallel
  {
#pragma omp single nowait
    {
      for (int j=0;j<n;++j) {
	for (int i=0;i<m;++i) {
#pragma omp task
	  FLA_Axpy( FLA_ONE, hA(i,j), hB(i,j) );

	}
      }
    }
  }
  t_block =  FLA_Clock() - t;
  
  printf("t_block %E, t_unblock %E, diff %lf\n",
	 t_block, t_unblock, t_block/t_unblock*100);

  A.free(); hA.free();
  B.free(); hB.free();

  FLA_Finalize();
  return 0;
}



