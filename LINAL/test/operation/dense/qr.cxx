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

#define BMN 192
#define N 500
#define LINAL_ERROR_TOL 1.0e-6

int main(int argc, char **argv) {
  FLA_Init();

  linal::Flat_ A1, A2, B, X1, X2, T1, T2, norm;
  linal::Hier_ hA2, hT2;;



  if (argc != 2) {
    printf("Try :: qr [thread]\n");
    return -1;
  }
  int datatype = LINAL_DOUBLE_REAL;
  int nthread = atoi( (argv[1]) );

  FLASH_Queue_set_num_threads(nthread);
  FLASH_Queue_set_sorting(TRUE);
  FLASH_Queue_set_caching(TRUE);

  omp_set_num_threads(nthread);
  printf("nthread from omp %d\n", omp_get_num_threads() );
  int ratio = 1;
  A1.create   (datatype, ratio*N, N);
  A2.create   (datatype, ratio*N, N);
  hA2.create  (A2, BMN, BMN);

  T1.create   (datatype,  N, N);
  T2.create   (datatype,  N, N);
  hT2.create  (T2, BMN, BMN);

  B.create    (datatype, ratio*N, 1);

  X1.create   (datatype, N, 1);
  X2.create   (datatype, N, 1);

  norm.create(datatype, 1, 1);

  FLA_Random_matrix(~A1);
  FLA_Copy(~A1, ~A2);
  FLA_Random_matrix(~B);
  
  double t_base, t_fla_qr, t_linal_qr;

  t_base = FLA_Clock();
  FLA_QR_UT(~A1,~T1);
  t_fla_qr = FLA_Clock() - t_base;


  //A1.disp("A1");
  //T1.disp("T1");
  printf("FLA QR DONE\n");

  t_base = FLA_Clock();

#pragma omp parallel
  {
#pragma omp single nowait
    {
#pragma omp task
      linal::dense::qr( hA2, hT2 );
    }
  }
  t_linal_qr = FLA_Clock() - t_base;
  printf("LINAL QR DONE\n");

  //A2.disp("A2");
  //T2.disp("T2");

  double normA, normT;
  FLA_Axpy(FLA_MINUS_ONE, ~A1, ~A2);
  FLA_Norm1( ~A2, ~norm);
  normA = norm(0,0);

  FLA_Axpy(FLA_MINUS_ONE, ~T1, ~T2);
  FLA_Norm1( ~T2, ~norm);
  normT = norm(0,0);

  int rval;
  printf("- TEST::");
  for (int i=0;i<2;++i)
    printf(" %s ", argv[i] );
  printf("\n");

  printf("fla %8.3lf linal %8.3lf\n", t_fla_qr, t_linal_qr);

  if ( (normA+normT) < LINAL_ERROR_TOL ) {
    printf("PASS::Norm :: A %E T %E \n", normA, normT);
    rval = 0;
  } else {
    printf("FAIL::Norm :: A %E T %E \n", normA, normT);
    rval = -1;
  }

  norm.free();

  X2.free();  X1.free();
  B.free();

  hT2.free(); T2.free(); T1.free();
  hA2.free(); A2.free(); A1.free();
      

  FLA_Finalize();

  return rval;
}
