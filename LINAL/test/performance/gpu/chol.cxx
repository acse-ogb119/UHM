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
#include "gpu_perform.hxx"

int main(int argc, char **argv) {

  if (argc != 6) {
    printf("Try :: chol [thread] [ngpu] [blocksize] [ndof] [nitr]\n");
    return 0;
  }

  // ---------------------------------------
  // ** Initialization
  int ndof, blocksize, nthread, nitr, ngpu;
  nthread   = atoi( (argv[1]) );
  ngpu      = atoi( (argv[2]) );
  blocksize = atoi( (argv[3]) );
  ndof      = atoi( (argv[4]) );
  nitr      = atoi( (argv[5]) );

  printf("** TEST ENVIRONMENT **\n");
  printf("NDOF      = %d\n", ndof);
  printf("Blocksize = %d\n", blocksize);
  printf("N thread  = %d\n", nthread);
  printf("N gpu     = %d\n", ngpu);
  printf("Iteration = %d\n", nitr);

  int b_mn[2];
  b_mn[0] = b_mn[1] = blocksize;

  FLA_Init();

  double flop = linal::get_flop_chol( 0, ndof );

  // ---------------------------------------
  // ** Matrices
  double 
    t_base, t_flash_repack, t_linal_repack, t_temp, 
    t_flash_decompose, t_linal_decompose;

  linal::Flat_ A;
  A.create(TEST_DATATYPE, ndof, ndof);
  FLA_Random_matrix( ~A );

  FLA_Obj hA_fla;
  t_base = FLA_Clock();
  FLASH_Obj_create_hier_copy_of_flat(~A, 1, (dim_t*)b_mn, &hA_fla);
  t_flash_repack = FLA_Clock()-t_base;

  linal::Hier_ hA_linal;
  t_base = FLA_Clock();
  hA_linal.create(A, blocksize, blocksize);
  t_linal_repack = FLA_Clock()-t_base;
  
  // ---------------------------------------
  // ** FLASH
  FLASH_Queue_set_num_threads(nthread);
  FLASH_Queue_set_sorting(TRUE);
  FLASH_Queue_set_caching(TRUE);
  FLASH_Queue_enable_gpu();
  // FLASH_Queue_set_gpu_num_blocks(blocksize);

  t_flash_decompose = MAX_TIME;
  for (int q=0;q<nitr;++q) {
    printf("*** FLASH::LU CHOL BEGIN ***\n");
    t_base = FLA_Clock();
    FLASH_Queue_begin();
    FLASH_Chol(FLA_LOWER_TRIANGULAR, hA_fla);
    FLASH_Queue_end();
    t_temp = FLA_Clock()-t_base;
    printf("*** FLASH::LU CHOL END ***\n");

    t_flash_decompose = min(t_temp, t_flash_decompose);
  }

  FLASH_Queue_disable_gpu();

  // ---------------------------------------
  // ** LINAL
  switch (TEST_COMPUTING) {
  case LINAL_CPU_GPU:
    omp_set_num_threads(nthread+ngpu);
    break;
  default:
    omp_set_num_threads(nthread);
    break;
  }
  linal::init_gpu(nthread, ngpu);
  linal::set_computing_model( TEST_COMPUTING );
  linal::create_flat_gpu( TEST_DATATYPE, b_mn[0], b_mn[1] );

  t_linal_decompose = MAX_TIME;
  for (int q=0;q<nitr;++q) {
    printf("*** LINAL::LU CHOL BEGIN ***\n");
    t_base = FLA_Clock();
#pragma omp parallel 
    {
#pragma omp single nowait
      linal::dense::chol(FLA_LOWER_TRIANGULAR, hA_linal);
    }
    t_temp = FLA_Clock()-t_base;
    printf("*** LINAL::LU CHOL END ***\n");

    t_linal_decompose = min(t_temp, t_linal_decompose);
  }

  linal::set_computing_model( LINAL_CPU );
  linal::free_flat_gpu();
  linal::finalize_gpu();

  switch (TEST_COMPUTING) {
  case LINAL_CPU_GPU:  printf("CPU + GPU BLAS\n");    break;
  case LINAL_GPU:      printf("GPU BLAS\n");          break;
  default:	       printf("CPU ONLY\n");          break;
  }
	
  printf("----------------------------------------------\n");
  printf("*** Report Chol ***\n");
  printf("Ndof       = %d\n", ndof);
  printf("Blocksize  = %d\n", blocksize);
  printf("Nthread    = %d\n", nthread);
  printf("Ngpu       = %d\n", ngpu);
  printf("Niteration = %d\n", nitr);
  printf("----------------------------------------------\n");
  printf("Time Supermatrix   = %6.3lf [sec]\n", t_flash_decompose);
  printf("Time LINAL         = %6.3lf [sec]\n", t_linal_decompose);
  printf("----------------------------------------------\n");
  printf("FLOPS Supermatrix  = %6.3lf [Gflops]\n", 
	 flop/t_flash_decompose/1.0e9);
  printf("FLOPS LINAL        = %6.3lf [Gflops]\n", 
	 flop/t_linal_decompose/1.0e9);
  printf("----------------------------------------------\n");
  printf("With repacking cost \n");
  printf("FLOPS Supermatrix  = %6.3lf [Gflops]\n", 
	 flop/(t_flash_decompose+t_flash_repack)/1.0e9);
  printf("FLOPS LINAL        = %6.3lf [Gflops]\n", 
	 flop/(t_linal_decompose+t_linal_repack)/1.0e9);
  printf("----------------------------------------------\n");
  
  // ---------------------------------------
  // ** Matrix
  FLASH_Obj_free(&hA_fla);
  hA_linal.free();
  A.free();

  printf("*** TEST FINISHED ***\n");

  // ---------------------------------------
  // ** Finalization
  FLA_Finalize();
  return 0;
}
