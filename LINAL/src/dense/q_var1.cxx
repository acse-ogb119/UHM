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
    int apply_q_lhfc_var1(Hier_ A, Hier_ T, Hier_ W, Hier_ B) {
      if ( !B.get_m() || !B.get_n() ) return true;
      
      FLA_Obj ATL,   ATR,      A00, A01, A02, 
              ABL,   ABR,      A10, A11, A12,
                               A20, A21, A22;
      FLA_Obj TL,    TR,       T0,  T1,  T2;
      FLA_Obj BT,              B0,
              BB,              B1,
                               B2;

      FLA_Part_2x2( ~A,    &ATL, &ATR,
                           &ABL, &ABR,     0, 0, FLA_TL );
      FLA_Part_1x2( ~T,    &TL,  &TR,      0, FLA_LEFT );
      FLA_Part_2x1( ~B,    &BT, 
                           &BB,            0, FLA_TOP );

      while ( FLA_Obj_min_dim( ABR ) > 0 ){
        FLA_Repart_2x2_to_3x3( ATL,  ATR,   &A00,  &A01, &A02,
                                            &A10,  &A11, &A12,
                               ABL,  ABR,   &A20,  &A21, &A22,  1, 1, FLA_BR );
        FLA_Repart_1x2_to_1x3( TL,   TR,    &T0,  &T1, &T2,     1, FLA_RIGHT );
        FLA_Repart_2x1_to_3x1( BT,          &B0, 
                                            &B1, 
                               BB,          &B2,                1, FLA_BOTTOM );
        //------------------------------------------------------------
        Hier_  
          A_11(A11), A_21(A21), B_1(B1), B_2(B2), T_1(T1),
          T_1T, W_TL;

	T_1.extract(T_1T, 1, T_1.get_n(), 0, 0);
	W.extract(W_TL, 1, B_1.get_n(), 0, 0);

        Flat_ W_TL_flat = W_TL.flat(), B_1_flat = B_1.flat();

        W_TL_flat.set_m( B_1_flat.get_m() );
        W_TL_flat.set_n( B_1_flat.get_n() );

        FLA_Copyt( FLA_NO_TRANSPOSE, ~B_1_flat, ~W_TL_flat );
        
        trmm( FLA_LEFT, FLA_LOWER_TRIANGULAR,
              FLA_CONJ_TRANSPOSE, FLA_UNIT_DIAG,
              FLA_ONE, A_11, W_TL );

        gemm( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, 
              FLA_ONE, A_21, B_2, FLA_ONE, W_TL );

        trsm( FLA_LEFT, FLA_UPPER_TRIANGULAR,
              FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
              FLA_ONE, T_1T, W_TL );
        // -----
	gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
	      FLA_MINUS_ONE, A_21, W_TL, FLA_ONE, B_2 );

	trmm( FLA_LEFT, FLA_LOWER_TRIANGULAR,
	      FLA_NO_TRANSPOSE, FLA_UNIT_DIAG,
	      FLA_MINUS_ONE, A_11, W_TL );
	FLA_Axpyt( FLA_NO_TRANSPOSE, FLA_ONE, ~W_TL_flat, ~B_1_flat );

        //------------------------------------------------------------
        FLA_Cont_with_3x3_to_2x2( &ATL,  &ATR,   A00, A01,  A02,
                                                 A10, A11,  A12,
                                  &ABL,  &ABR,   A20, A21,  A22,      FLA_TL );
        FLA_Cont_with_1x3_to_1x2( &TL,   &TR,    T0, T1,  T2,         FLA_LEFT );
        FLA_Cont_with_3x1_to_2x1( &BT,           B0, 
                                                 B1, 
                                  &BB,           B2,                  FLA_TOP );
      }
      return true;
    }
    
    int apply_q2_lhfc_var1(Hier_ D, Hier_ T, Hier_ W, Hier_ C, Hier_ E) {
      FLA_Obj DL,    DR,       D0,  D1,  D2;
      FLA_Obj TL,    TR,       T0,  T1,  T2;
      FLA_Obj WL,    WR,       W0,  W1,  W2;
      FLA_Obj CT,              C0,
              CB,              C1,
                               C2;

      FLA_Part_1x2( ~D,    &DL,  &DR,      0, FLA_LEFT );
      FLA_Part_1x2( ~T,    &TL,  &TR,      0, FLA_LEFT );
      FLA_Part_1x2( ~W,    &WL,  &WR,      0, FLA_LEFT );
      FLA_Part_2x1( ~C,    &CT, 
                           &CB,            0, FLA_TOP );

      while ( FLA_Obj_width( DL ) < FLA_Obj_width( ~D ) ){
        FLA_Repart_1x2_to_1x3( DL,  DR,      &D0,  &D1, &D2, 1, FLA_RIGHT );
        FLA_Repart_1x2_to_1x3( TL,  TR,      &T0,  &T1, &T2, 1, FLA_RIGHT );
        FLA_Repart_1x2_to_1x3( WL,  WR,      &W0,  &W1, &W2, 1, FLA_RIGHT );
        FLA_Repart_2x1_to_3x1( CT,           &C0, 
                                             &C1, 
                               CB,           &C2,            1, FLA_BOTTOM );
        //------------------------------------------------------------
        //Hier_ C_1(C1), W_1(W1), D_1(D1);
        FLA_Copyt( FLA_NO_TRANSPOSE, C1, W1 );
        gemm_internal( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, 
                       FLA_ONE, D1, ~E, FLA_ONE, W1 );

        trsm_internal( FLA_LEFT, FLA_UPPER_TRIANGULAR,
                       FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
                       FLA_ONE, T1, W1 );

        FLA_Axpyt( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, W1, C1 );

        gemm_internal( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                       FLA_MINUS_ONE, D1, W1, FLA_ONE, ~E );

        //------------------------------------------------------------
        FLA_Cont_with_1x3_to_1x2( &DL,   &DR,   D0, D1,  D2, FLA_LEFT );
        FLA_Cont_with_1x3_to_1x2( &TL,   &TR,   T0, T1,  T2, FLA_LEFT );
        FLA_Cont_with_1x3_to_1x2( &WL,   &WR,   W0, W1,  W2, FLA_LEFT );
        FLA_Cont_with_3x1_to_2x1( &CT,          C0, 
                                                C1, 
                                  &CB,          C2,          FLA_TOP );

      }        
      return true;
    }
  }
}
