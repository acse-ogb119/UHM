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

    // Size of T drive the QR_inc
    // The dimension of A and T should be same
    int qr_inc_var1(Hier_ A, Hier_ T) {
      if (!A.get_m() || !A.get_n()) return true;
      LINAL_ERROR( (A.get_m() == T.get_m()) &&
                   (A.get_n() == T.get_n()) ,
                  ">> Hier_ A and T should have same dimension");

      FLA_Obj ATL,   ATR,      A00, A01, A02, 
              ABL,   ABR,      A10, A11, A12,
                               A20, A21, A22;

      FLA_Obj TTL,   TTR,      T00, T01, T02, 
              TBL,   TBR,      T10, T11, W12,
                               T20, T21, T22;

      FLA_Part_2x2(~A,    &ATL, &ATR,
                          &ABL, &ABR,     0, 0, FLA_TL );

      FLA_Part_2x2(~T,    &TTL, &TTR,
                          &TBL, &TBR,     0, 0, FLA_TL );

      // factorize based on the dimension of TBR 
      while ( FLA_Obj_min_dim( TBR ) > 0 ){
        FLA_Repart_2x2_to_3x3( ATL,  ATR,       &A00,  &A01, &A02,
                                                &A10,  &A11, &A12,
                               ABL,  ABR,       &A20,  &A21, &A22,
                               1, 1, FLA_BR );

        FLA_Repart_2x2_to_3x3( TTL,  TTR,       &T00,  &T01, &T02,
                                                &T10,  &T11, &W12,
                               TBL,  TBR,       &T20,  &T21, &T22,
                               1, 1, FLA_BR );

        //------------------------------------------------------------
        Hier_ 
          A_11(A11), A_21(A21), A_12(A12), A_22(A22), 
          T_11(T11), T_21(T21), W_12(W12); 

        int b = FLA_Obj_min_dim( A_11(0,0) );

        // dimension of T is modified at this time
        T_11(0,0).m = b;
        T_11(0,0).n = b;
        
        FLA_QR_UT( A_11(0,0), T_11(0,0) );

	apply_q( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                 A_11, T_11, W_12, A_12);

        if (T_21.get_m() && T_21.get_n()) {
	  
	  for (int i=0;i<T_21.get_m();++i) {

            // make T_21 square
            T_21(i,0).m = b;
            T_21(i,0).n = b;

	    FLA_QR2_UT_internal( A_11(0,0),
				 A_21(i,0), T_21(i,0),
				 fla_qr2ut_cntl_leaf );


	    for (int j=0;j<W_12.get_n();++j) {

#pragma omp task firstprivate(A_21, T_21, W_12, A_12, A_22, i,j)
              {
                FLA_Apply_Q2_UT_internal( FLA_LEFT, FLA_CONJ_TRANSPOSE,
                                          FLA_FORWARD, FLA_COLUMNWISE,
                                          A_21(i,0), T_21(i,0),
                                          W_12(0,j), A_12(0,j),
                                          A_22(i,j),
                                          fla_apq2ut_cntl_leaf );
              }
            }

#pragma omp taskwait

	  }
	}
	
        //------------------------------------------------------------
        FLA_Cont_with_3x3_to_2x2( &ATL,  &ATR,       A00, A01,  A02,
                                                     A10, A11,  A12,
                                  &ABL,  &ABR,       A20, A21,  A22,
                                  FLA_TL );
        FLA_Cont_with_3x3_to_2x2( &TTL,  &TTR,       T00, T01,  T02,
                                                     T10, T11,  W12,
                                  &TBL,  &TBR,       T20, T21,  T22,
                                  FLA_TL );
      }

      return true;
    }
  }
}
  


	  /*
	  {
	    for (int i=0;i<T_21.get_m();++i) {
	      //#pragma omp task firstprivate(i)
	      FLA_QR2_UT_internal( A_11(0,0),
				   A_21(i,0), T_21(i,0),
				   fla_qr2ut_cntl_leaf );
              
	      for (int j=0;j<W_12.get_n();++j) {
		Flat_ T_21_i0(T_21(i,0));
		T_21_i0.set_n( T_21_i0.get_m() );
		
#pragma omp task firstprivate(i,j)
		FLA_Apply_Q2_UT_internal( FLA_LEFT, FLA_CONJ_TRANSPOSE,
					  FLA_FORWARD, FLA_COLUMNWISE,
					  A_21(i,0), ~T_21_i0,
					  W_12(0,j), A_12(0,j),
					  A_22(i,j),
					  fla_apq2ut_cntl_leaf );
	      }
#pragma omp taskwait
            }
          }
	  */
