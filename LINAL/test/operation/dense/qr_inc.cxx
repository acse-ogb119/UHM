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

#define BMN 96
#define N 5000
#define LINAL_ERROR_TOL 1.0e-6

int main(int argc, char **argv) {
  FLA_Init();

  linal::Flat_ A1, A2, B1, B2, X1, X2, T1, T2, W2, norm;
  linal::Hier_ hA1, hT1, hA2, hT2, hB1, hB2, hX1, hX2, hW2;
  FLA_Obj hA1_fla, hT1_fla, hX1_fla, hW1_fla, hB1_fla;

  
  if (argc != 2) {
    printf("Try :: qr_inc [thread]\n");
    return -1;
  }

  int datatype = LINAL_DOUBLE_REAL;
  int nthread = atoi( (argv[1]) );
  int blocksize = BMN;

  FLASH_Queue_set_num_threads(nthread);
  FLASH_Queue_set_sorting(TRUE);
  FLASH_Queue_set_caching(TRUE);

  omp_set_num_threads(nthread);

  int b_mn[2];
  b_mn[0] = b_mn[1] = blocksize;
  
  int ratio = 1;
  A1.create   (datatype, ratio*N, N);
  A2.create   (datatype, ratio*N, N);
  hA1.create  (A1, BMN, BMN);
  hA2.create  (A2, BMN, BMN);
  
  int T_size = (N/BMN+(N%BMN!=0))*BMN;
  T1.create   (datatype,  T_size, T_size);
  T2.create   (datatype,  T_size, T_size);
  hT1.create  (T1, BMN, BMN);
  hT2.create  (T2, BMN, BMN);

  int n_rhs = 1;
  B1.create    (datatype, ratio*N, n_rhs);
  hB1.create   (B1, BMN, BMN);

  B2.create    (datatype, ratio*N, n_rhs);
  hB2.create   (B2, BMN, BMN);

  X1.create   (datatype, N, n_rhs);
  hX1.create  (X1, BMN, BMN);

  X2.create   (datatype, N, n_rhs);
  hX2.create  (X2, BMN, BMN);


  W2.create   (datatype, BMN, BMN*hB2.get_n() );
  hW2.create  (W2, BMN, BMN);


  norm.create(datatype, 1, 1);

  FLA_Random_matrix(~A1);
  FLA_Copy(~A1, ~A2);

  FLA_Random_matrix(~B1);
  FLA_Copy(~B1, ~B2);

  //B1.disp("reference");
  
  FLASH_QR_UT_inc_create_hier_matrices( ~A1, 1, (dim_t*)b_mn, (dim_t)blocksize,
                                        &hA1_fla, &hT1_fla );

  FLASH_Obj_create_hier_copy_of_flat(~B1, 1, (dim_t*)b_mn, &hB1_fla);
  FLASH_Obj_create_hier_copy_of_flat(~X1, 1, (dim_t*)b_mn, &hX1_fla);
  FLASH_Apply_Q_UT_inc_create_workspace( hT1_fla, hB1_fla, &hW1_fla );

  /*
  linal::Hier_ A_tmp(hA1_fla), T_tmp(hT1_fla), B_tmp(hB1_fla), W_tmp(hW1_fla);

  //A_tmp.disp();
  //T_tmp.disp();
  //B_tmp.disp();
  W_tmp.disp();
  hW2.disp();
  */
  double t_base, t_fla_qr, t_linal_qr;

  t_base = FLA_Clock();
  //FLASH_Queue_begin();
  FLASH_QR_UT_inc( hA1_fla,hT1_fla );
  FLASH_QR_UT_inc_solve( hA1_fla, hT1_fla, hB1_fla, hX1_fla );
  FLASH_Apply_Q_UT_inc( FLA_LEFT, FLA_CONJ_TRANSPOSE,
                        FLA_FORWARD, FLA_COLUMNWISE,
                        hA1_fla, hT1_fla, hW1_fla, hB1_fla );

  FLASH_Apply_Q_UT_inc( FLA_LEFT, FLA_NO_TRANSPOSE,
                        FLA_FORWARD, FLA_COLUMNWISE,
                        hA1_fla, hT1_fla, hW1_fla, hB1_fla );

  //FLASH_Queue_end();
  t_fla_qr = FLA_Clock() - t_base;

  FLASH_Obj_flatten( hA1_fla, ~A1 );
  FLASH_Obj_flatten( hT1_fla, ~T1 );
  FLASH_Obj_flatten( hB1_fla, ~B1 );
  FLASH_Obj_flatten( hX1_fla, ~X1 );

  //B1.disp("FLASH");
  
  //A1.disp("A1");
  //T1.disp("T1");
  printf("FLA QR DONE\n");

  t_base = FLA_Clock();
#pragma omp parallel
  {
#pragma omp single nowait
    {
#pragma omp task
      linal::dense::qr_inc( hA2, hT2 );

    }
  }

  //B2.disp("LINAL begin");

#pragma omp parallel
  {
#pragma omp single nowait
    {
#pragma omp task
      linal::dense::apply_q_inc( FLA_LEFT, FLA_CONJ_TRANSPOSE,
                                 FLA_FORWARD, FLA_COLUMNWISE,
                                 hA2, hT2, hW2, hB2 );

    }
  }

  //B2.disp("LINAL int");

  t_linal_qr = FLA_Clock() - t_base;
  printf("LINAL QR DONE\n");

  FLA_Copy(~B2, ~X2);

#pragma omp parallel
  {
#pragma omp single nowait
    {
#pragma omp task
      linal::dense::apply_q_inc( FLA_LEFT, FLA_NO_TRANSPOSE,
                                 FLA_FORWARD, FLA_COLUMNWISE,
                                 hA2, hT2, hW2, hB2 );
    }
  }

  //B2.disp("LINAL");

#pragma omp parallel
  {
#pragma omp single nowait
    {
#pragma omp task
  linal::dense::trsm( FLA_LEFT, FLA_UPPER_TRIANGULAR, 
                      FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
                      FLA_ONE, hA2, hX2 );
    }
  }
  /*
  B1.disp("B1 after");  
  B2.disp("B2 after");

  X1.disp("X1 solved");
  X2.disp("X2 solved");

  //A2.disp("A2");
  //T2.disp("T2");
  */
  double normA, normT, normB, normX;
  FLA_Axpy(FLA_MINUS_ONE, ~A1, ~A2);
  FLA_Norm1( ~A2, ~norm);
  normA = norm(0,0);

  FLA_Axpy(FLA_MINUS_ONE, ~T1, ~T2);
  FLA_Norm1( ~T2, ~norm);
  normT = norm(0,0);

  FLA_Axpy(FLA_MINUS_ONE, ~B1, ~B2);
  FLA_Norm1( ~B2, ~norm);
  normB = norm(0,0);

  FLA_Axpy(FLA_MINUS_ONE, ~X1, ~X2);
  FLA_Norm1( ~X2, ~norm);
  normX = norm(0,0);

  int rval;
  printf("- TEST::");
  for (int i=0;i<2;++i)
    printf(" %s ", argv[i] );
  printf("\n");

  printf("fla %8.3lf linal %8.3lf\n", t_fla_qr, t_linal_qr);

  if ( (normA+normT+normB+normX) < LINAL_ERROR_TOL ) {
    printf("PASS::Norm :: A %E T %E B %E X %E\n", normA, normT, normB, normX);
    rval = 0;
  } else {
    printf("FAIL::Norm :: A %E T %E B %E X %E\n", normA, normT, normB, normX);
    rval = -1;
  }
  FLASH_Obj_free(&hW1_fla);
  FLASH_Obj_free(&hB1_fla);
  FLASH_Obj_free(&hT1_fla);
  FLASH_Obj_free(&hA1_fla);
  
  norm.free();
  hW2.free();  W2.free();
  hX2.free();  X2.free();
  hX1.free();  X1.free();
  hB2.free();  B2.free();
  hB1.free();  B1.free();

  hT2.free();  T2.free();
  hT1.free();  T1.free();
  hA2.free();  A2.free();
  hA1.free();  A1.free();


  FLA_Finalize();

  return rval;
}
