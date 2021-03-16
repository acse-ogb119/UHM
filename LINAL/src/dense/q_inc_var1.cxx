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
    // ------------------------------------------------------------
    // ** left, conj_transpose, forward, columnwise
    // ------------------------------------------------------------
    int apply_q_inc_lhfc_var1(Hier_ A, Hier_ T, Hier_ W_1, Hier_ B) {
      if ( !B.get_m() || !B.get_n() ) return true;
      LINAL_ERROR( (A.get_m() == T.get_m()) &&
                   (A.get_n() == T.get_n()),
                   ">> Hier_ A and T should have same dimension");
      LINAL_ERROR( (A.get_m() == B.get_m()),
                   ">> Hier_ A and B should have same length");
      LINAL_ERROR( (W_1.get_n() == B.get_n()),
                   ">> Hier_ W_1 and B should have same width");
      
      FLA_Obj ATL,   ATR,      A00, A01, A02, 
              ABL,   ABR,      A10, A11, A12,
                               A20, A21, A22;

      FLA_Obj TTL,   WTR,      T00, W01, W02, 
              TBL,   TBR,      T10, T11, W12,
                               T20, T21, T22;

      FLA_Obj BT,              B0,
              BB,              B1,
                               B2;

      FLA_Part_2x2(~A,    &ATL, &ATR,
                          &ABL, &ABR,     0, 0, FLA_TL );

      FLA_Part_2x2(~T,    &TTL, &WTR,
                          &TBL, &TBR,     0, 0, FLA_TL );

      FLA_Part_2x1(~B,    &BT,
                          &BB,            0, FLA_TOP );

      while ( FLA_Obj_min_dim( ABR ) > 0 ){

        FLA_Repart_2x2_to_3x3( ATL,  ATR,       &A00,  &A01, &A02,
                                                &A10,  &A11, &A12,
                               ABL,  ABR,       &A20,  &A21, &A22,
                               1, 1, FLA_BR );
        FLA_Repart_2x2_to_3x3( TTL,  WTR,       &T00,  &W01, &W02,
                                                &T10,  &T11, &W12,
                               TBL,  TBR,       &T20,  &T21, &T22,
                               1, 1, FLA_BR );
        FLA_Repart_2x1_to_3x1( BT,                &B0,
                                                  &B1,
                               BB,                &B2,        1, FLA_BOTTOM );
        //------------------------------------------------------------
        Hier_ 
          A_11(A11), A_21(A21), 
          T_11(T11), T_21(T21),
          B_1(B1), B_2(B2);

        apply_q( FLA_LEFT, FLA_CONJ_TRANSPOSE, 
                 FLA_FORWARD, FLA_COLUMNWISE,
                 A_11, T_11, W_1, B_1 );

        // do not change the order of loop
	for (int i=0;i<B_2.get_m();++i) {                 
	  for (int j=0;j<B_2.get_n();++j) {

#pragma omp task firstprivate(A_21, T_21, W_1, B_1, B_2, i, j) 
            FLA_Apply_Q2_UT_internal( FLA_LEFT, FLA_CONJ_TRANSPOSE, 
                                      FLA_FORWARD, FLA_COLUMNWISE,
                                      A_21(i,0), T_21(i,0),
                                      W_1(0,j),  B_1(0,j),
                                      B_2(i,j),
                                      fla_apq2ut_cntl_leaf );

          }

#pragma omp taskwait	  

        }

        //------------------------------------------------------------

        FLA_Cont_with_3x3_to_2x2( &ATL,  &ATR,       A00, A01,  A02,
                                                     A10, A11,  A12,
                                  &ABL,  &ABR,       A20, A21,  A22,
                                  FLA_TL );
        FLA_Cont_with_3x3_to_2x2( &TTL,  &WTR,       T00, W01,  W02,
                                                     T10, T11,  W12,
                                  &TBL,  &TBR,       T20, T21,  T22,
                                  FLA_TL );
        FLA_Cont_with_3x1_to_2x1( &BT,                B0,
                                                      B1,
                                  &BB,                B2,     FLA_TOP );
      }
      return true;
    }

    // ------------------------------------------------------------
    // ** left, no_transpose, forward, columnwise
    // ------------------------------------------------------------
    int apply_q_inc_lnfc_var1(Hier_ A, Hier_ T, Hier_ W_1, Hier_ B) {
      if ( !B.get_m() || !B.get_n() ) return true;

      LINAL_ERROR( (A.get_m() == T.get_m()) &&
                   (A.get_n() == T.get_n()),
                   ">> Hier_ A and T should have same dimension");
      LINAL_ERROR( (A.get_m() == B.get_m()),
                   ">> Hier_ A and B should have same length");
      LINAL_ERROR( (W_1.get_n() == B.get_n()),
                   ">> Hier_ W_1 and B should have same width");
      
      FLA_Obj ATL,   ATR,      A00, A01, A02, 
              ABL,   ABR,      A10, A11, A12,
                               A20, A21, A22;

      FLA_Obj TTL,   WTR,      T00, W01, W02, 
              TBL,   TBR,      T10, T11, W12,
                               T20, T21, T22;

      FLA_Obj BT,              B0,
              BB,              B1,
                               B2;

      FLA_Part_2x2(~A,    &ATL, &ATR,
                          &ABL, &ABR,     0, 0, FLA_BR );

      FLA_Part_2x2(~T,    &TTL, &WTR,
                          &TBL, &TBR,     0, 0, FLA_BR );

      FLA_Part_2x1(~B,    &BT,
                          &BB,            0, FLA_BOTTOM );

      while ( FLA_Obj_min_dim( ATL ) > 0 ){

        FLA_Repart_2x2_to_3x3( ATL,  ATR,       &A00,  &A01, &A02,
                                                &A10,  &A11, &A12,
                               ABL,  ABR,       &A20,  &A21, &A22,
                               1, 1, FLA_TL );
        FLA_Repart_2x2_to_3x3( TTL,  WTR,       &T00,  &W01, &W02,
                                                &T10,  &T11, &W12,
                               TBL,  TBR,       &T20,  &T21, &T22,
                               1, 1, FLA_TL );
        FLA_Repart_2x1_to_3x1( BT,                &B0,
                                                  &B1,
                               BB,                &B2,        
                               1, FLA_TOP );
        //------------------------------------------------------------
        Hier_ 
          A_11(A11), A_21(A21), 
          T_11(T11), T_21(T21),
          B_1(B1), B_2(B2);

        // do not change the order of loop
	for (int i=(B_2.get_m() - 1);i>=0;--i) {
	  for (int j=(B_2.get_n() - 1);j>=0;--j) {
            
#pragma omp task firstprivate(A_21, T_21, W_1, B_1, B_2, i, j) 
            FLA_Apply_Q2_UT_internal( FLA_LEFT, FLA_NO_TRANSPOSE, 
                                      FLA_FORWARD, FLA_COLUMNWISE,
                                      A_21(i,0), T_21(i,0),
                                      W_1(0,j),  B_1(0,j),
                                      B_2(i,j),
                                      fla_apq2ut_cntl_leaf );

          }

#pragma omp taskwait	  

        }

        apply_q( FLA_LEFT, FLA_NO_TRANSPOSE, 
                 FLA_FORWARD, FLA_COLUMNWISE,
                 A_11, T_11, W_1, B_1 );

        //------------------------------------------------------------

        FLA_Cont_with_3x3_to_2x2( &ATL,  &ATR,       A00, A01,  A02,
                                                     A10, A11,  A12,
                                  &ABL,  &ABR,       A20, A21,  A22,
                                  FLA_BR );
        FLA_Cont_with_3x3_to_2x2( &TTL,  &WTR,       T00, W01,  W02,
                                                     T10, T11,  W12,
                                  &TBL,  &TBR,       T20, T21,  T22,
                                  FLA_BR );
        FLA_Cont_with_3x1_to_2x1( &BT,                B0,
                                                      B1,
                                  &BB,                B2,     
                                  FLA_BOTTOM );
      }
      return true;
    }
  }
}
