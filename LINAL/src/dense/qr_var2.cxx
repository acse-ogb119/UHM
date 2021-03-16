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
    int qr_var2(Hier_ A, Hier_ T) {
      if (!A.get_m() || !A.get_n()) return true;

      printf("qr_var2\n");
      //FLA_QR_UT_blk_var2( ~A.flat(), ~T.flat(), fla_qrut_cntl_leaf ); 
      //return true;


      FLA_Obj ATL,   ATR,      A00, A01, A02, 
              ABL,   ABR,      A10, A11, A12,
                               A20, A21, A22;
      FLA_Obj TL,    TR,       T0,  T1,  W12;

      FLA_Part_2x2(~A,    &ATL, &ATR,
                          &ABL, &ABR,     0, 0, FLA_TL );

      FLA_Part_1x2(~T,    &TL,  &TR,      0, FLA_LEFT );

      while ( FLA_Obj_min_dim( ABR ) > 0 ){
        FLA_Repart_2x2_to_3x3( ATL,  ATR,       &A00,  &A01, &A02,
                                                &A10,  &A11, &A12,
                               ABL,  ABR,       &A20,  &A21, &A22,
                               1, 1, FLA_BR );

        FLA_Repart_1x2_to_1x3( TL,   TR,        &T0,  &T1, &W12,
                               1, FLA_RIGHT );

        //------------------------------------------------------------
        FLA_Obj T1T, T2B;
        FLA_Part_2x1( T1,   &T1T, 
                            &T2B,   1, FLA_TOP );

        FLA_Obj AB1;
        FLA_Merge_2x1( A11,
                       A21,   &AB1 );

        Hier_ AB_1(AB1), T_1T(T1T);
        FLA_QR_UT( ~AB_1.flat(), ~T_1T.flat() );

        if ( FLA_Obj_width( A12 ) > 0 ) {
          FLA_Obj AB2;
          FLA_Merge_2x1( A12,
                         A22,   &AB2 );

          Hier_ AB_2(AB2), W_12(W12);
          apply_q( FLA_LEFT, FLA_CONJ_TRANSPOSE, 
                   FLA_FORWARD, FLA_COLUMNWISE,
                   AB_1, T_1T, W_12, AB_2 );

        }

#pragma omp taskwait
                
        //------------------------------------------------------------
        FLA_Cont_with_3x3_to_2x2( &ATL,  &ATR,       A00, A01,  A02,
                                                     A10, A11,  A12,
                                  &ABL,  &ABR,       A20, A21,  A22,
                                  FLA_TL );

        FLA_Cont_with_1x3_to_1x2( &TL,   &TR,        T0, T1,  W12,
                                  FLA_LEFT );
      }
      return true;
    }
  }
}
