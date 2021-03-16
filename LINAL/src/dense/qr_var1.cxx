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
    int qr_var1(Hier_ A, Hier_ T) {
      if (!A.get_m() || !A.get_n()) return true;

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

      while ( FLA_Obj_min_dim( ABR ) > 0 ){
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
          A_10(A10), A_11(A11), A_20(A20), A_21(A21), 
          T_01(T01), T_11(T11), W_12(W12); 

        FLA_Obj AB1;
        FLA_Merge_2x1( A11,
                       A21,   &AB1 );

        Hier_ AB_1(AB1);
        FLA_QR_UT( ~AB_1.flat(), ~T_11.flat() );

        if ( FLA_Obj_width( A12 ) > 0 )  {
          FLA_Obj AB2;
          FLA_Merge_2x1( A12,
                         A22,   &AB2 );

          Hier_ AB_2(AB2);

#pragma omp task
          {
            apply_q( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                     AB_1, T_11, W_12, AB_2);
          }
        }

        if ( A_10.get_m() && A_10.get_n() ) 
          FLA_Copyt( FLA_CONJ_TRANSPOSE, ~A_10.flat(), ~T_01.flat() );

#pragma omp task
        {
          trmm( FLA_RIGHT, FLA_LOWER_TRIANGULAR,
                FLA_NO_TRANSPOSE, FLA_UNIT_DIAG,
                FLA_ONE, A_11, T_01 );
          
          gemm( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, 
                FLA_ONE, A_20, A_21, FLA_ONE, T_01 );
        }
        
#pragma omp taskwait
	
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

            FLA_Apply_Q_UT( FLA_LEFT, FLA_CONJ_TRANSPOSE, 
			    FLA_FORWARD, FLA_COLUMNWISE,
			    ~AB_1.flat(), ~T_11.flat(), 
			    ~W_12.flat(), ~AB_2.flat());
	    



        if ( A_10.get_m() && A_10.get_n() ) {

          FLA_Copyt( FLA_CONJ_TRANSPOSE, ~A_10.flat(), ~T_01.flat() );

          FLA_Trmm( FLA_RIGHT, FLA_LOWER_TRIANGULAR,
                    FLA_NO_TRANSPOSE, FLA_UNIT_DIAG,
                    FLA_ONE, ~A_11.flat(), ~T_01.flat() );

          if ( A_20.get_m() && A_20.get_n() ) 
            FLA_Gemm( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, 
                      FLA_ONE, ~A_20.flat(), ~A_21.flat(), 
                      FLA_ONE, ~T_01.flat() );
        }
        */
