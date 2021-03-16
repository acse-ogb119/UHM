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
#include "linal/operation.hxx"


namespace linal {
  namespace dense {
    /*!
      General matrix matrix multiplication.
    */
    
    int gemm( int transa, int transb,
              FLA_Obj alpha, Hier_ A, Hier_ B,
              FLA_Obj beta,  Hier_ C ) {


      // for complex variable, 
      // only FLA_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE 
      // are allowed

      int val = transa*1000 + transb;
      // Use GEPP algorithm (-- use block rank one update)
      // scalar beta should be applied first otherwise sign is floping 
      scal( beta, C );

      switch (val) {
      case 400400:  // nt nt
        for (int p=0;p<A.get_n();++p) {
          for (int k2=0;k2<C.get_n();++k2) {
            for (int k1=0;k1<C.get_m();++k1) {

#pragma omp task firstprivate ( k1, k2, p )
              {
                gemm_internal( transa, transb,
                               alpha,   A(k1,p), B(p,k2),
                               FLA_ONE, C(k1,k2) );
              }
            }
          }
#pragma omp taskwait
        }
        break;
      case 401401:  // t t
      case 402402:  // conj_t conj_t
        for (int p=0;p<A.get_m();++p) {
          for (int k2=0;k2<C.get_n();++k2) {
            for (int k1=0;k1<C.get_m();++k1) {

#pragma omp task firstprivate ( k1, k2, p )
              {
                gemm_internal( transa, transb,
                               alpha,   A(p,k1), B(k2,p),
                               FLA_ONE, C(k1,k2) );
              }
            }
          }
#pragma omp taskwait
        }
        break;
      case 400401:  // nt t                                       
      case 400402:  // nt conj_t                                       
        for (int p=0;p<A.get_n();++p) {
          for (int k2=0;k2<C.get_n();++k2) {
            for (int k1=0;k1<C.get_m();++k1) {

#pragma omp task firstprivate ( k1, k2, p )
              {
                gemm_internal( transa, transb,
                               alpha,   A(k1,p), B(k2,p),
                               FLA_ONE, C(k1,k2) );
              }
            }
          }
#pragma omp taskwait
        }
        break;
      case 401400:  // t nt                               
      case 402400:  // conj_t nt                               
        for (int p=0;p<A.get_m();++p) {
          for (int k2=0;k2<C.get_n();++k2) {
            for (int k1=0;k1<C.get_m();++k1) {

#pragma omp task firstprivate ( k1, k2, p )
              {
                gemm_internal( transa, transb,
                               alpha,   A(p,k1), B(p,k2),
                               FLA_ONE, C(k1,k2) );
              }
            }
          }
#pragma omp taskwait
        }
        break;
      default: val = 0; break;
      }

      LINAL_ERROR(val,
                  ">> Input trans arguement is invalid");

      return true;
    }
  }
}
