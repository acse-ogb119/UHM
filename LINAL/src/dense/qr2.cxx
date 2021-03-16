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
    // FLA_QR2_UT_inc_blk_var1
    int qr2(Hier_ U, Hier_ D, Hier_ T) {
      FLA_Obj UTL,   UTR,      U00, U01, U02, 
              UBL,   UBR,      U10, U11, U12,
                               U20, U21, U22;

      FLA_Obj DL,    DR,       D0,  D1,  D2;
      FLA_Obj TL,    TR,       T0,  T1,  W12;

      FLA_Part_2x2( ~U,    &UTL, &UTR,
                           &UBL, &UBR,     0, 0, FLA_TL );
      FLA_Part_1x2( ~D,    &DL,  &DR,      0, FLA_LEFT );
      FLA_Part_1x2( ~T,    &TL,  &TR,      0, FLA_LEFT );

      while ( FLA_Obj_min_dim( UBR ) > 0 ){
        FLA_Repart_2x2_to_3x3( UTL,  UTR,       &U00,  &U01, &U02,
                                                &U10,  &U11, &U12,
                               UBL,  UBR,       &U20,  &U21, &U22,
                               1, 1, FLA_BR );
        FLA_Repart_1x2_to_1x3( DL,   DR,        &D0,  &D1, &D2,
                               1, FLA_RIGHT );
        FLA_Repart_1x2_to_1x3( TL,   TR,        &T0,  &T1, &W12,
                               1, FLA_RIGHT );
        //------------------------------------------------------------
        Hier_ U_11(U11), D_1(D1), T_1(T1);
        FLA_QR2_UT_internal( U11,
                             D1, T1, 
                             fla_qr2ut_cntl_leaf );

        if ( FLA_Obj_width( U12 ) > 0 ) {
          Hier_ U_12(U12), W_12(W12), D_2(D2);
          FLA_Copy( ~U_12.flat(), ~W_12.flat() );
          FLA_Gemm( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE,
                    FLA_ONE, D_1(0,0), ~D_2.flat(), FLA_ONE, ~W_12.flat() );

          FLA_Trsm( FLA_LEFT, FLA_UPPER_TRIANGULAR, 
                    FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG, 
                    FLA_ONE, T_1(0,0), ~W_12.flat() );
          
          FLA_Axpy( FLA_MINUS_ONE, ~W_12.flat(), ~U_12.flat() );

          FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, 
                    FLA_MINUS_ONE, D_1(0,0), ~W_12.flat(), FLA_ONE, ~D_2.flat() );
        }
        //------------------------------------------------------------
        FLA_Cont_with_3x3_to_2x2( &UTL,  &UTR,       U00, U01,  U02,
                                                     U10, U11,  U12,
                                  &UBL,  &UBR,       U20, U21,  U22,
                                  FLA_TL );
        FLA_Cont_with_1x3_to_1x2( &DL,   &DR,        D0, D1,  D2,
                                  FLA_LEFT );
        FLA_Cont_with_1x3_to_1x2( &TL,   &TR,        T0, T1,  W12,
                                  FLA_LEFT );
      }
      return true;
    }
  }
}
