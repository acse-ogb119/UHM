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
#include "linal/common.hxx"
#include "linal/const.hxx"
#include "linal/util.hxx"
#include "linal/matrix.hxx"
#include "linal/flat.hxx"
#include "linal/hier.hxx"

#ifdef LINAL_GPU_ENABLE
#include "cublas.h"
#include "linal/gpu.hxx"
#endif

#include "linal/operation.hxx"

namespace linal {

#ifdef LINAL_GPU_ENABLE  
  inline static void trmm_internal_gpu( int side, int uplo, int trans, int diag, 
                                        FLA_Obj alpha, FLA_Obj A, FLA_Obj B ) {
    int thread   = omp_get_thread_num();

    Flat_GPU_ A_gpu, B_gpu;

    get_flat_gpu( thread, 0 ).extract( A_gpu, A.m, A.n );
    get_flat_gpu( thread, 1 ).extract( B_gpu, B.m, B.n );
    
    A_gpu.set_matrix_from(A);
    B_gpu.set_matrix_from(B);
    
    FLA_Trmm_external_gpu( side, uplo, trans, diag, alpha,
                           ~A_gpu, A_gpu.get_buffer(),
                           ~B_gpu, B_gpu.get_buffer() );
    
    A_gpu.get_matrix_to(A);
    B_gpu.get_matrix_to(B);
  }
#endif
  
  int trmm_internal( int side, int uplo, int trans, int diag, 
                     FLA_Obj alpha, FLA_Obj A, FLA_Obj B ) {
    switch ( get_computing_model() ) {
#ifdef LINAL_GPU_ENABLE
    case LINAL_GPU:
      trmm_internal_gpu( side, uplo, trans, diag, alpha, A, B );
      break;

    case LINAL_CPU_GPU:
      {
        //if ( omp_get_thread_num() >= get_n_thread() )

        int thread   = omp_get_thread_num();
        int n_thread = get_n_thread();

        if ( thread >= n_thread )
          FLA_Trmm( side, uplo, trans, diag, alpha, A, B );
        else
          trmm_internal_gpu( side, uplo, trans, diag, alpha, A, B );
      }
      break;
#endif

    default:
      FLA_Trmm( side, uplo, trans, diag, alpha, A, B );
      break;
    }
    return true;
  }
}
