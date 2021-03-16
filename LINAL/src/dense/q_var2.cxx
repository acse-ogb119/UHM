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
    int apply_q_l_var2( int trans, int direct, 
                        Hier_ A, Hier_ T, Hier_ W, Hier_ B ) {
      Flat_ A_flat, T_flat;
      A_flat = A.flat();      T_flat = T.flat();

      FLA_Obj BL,    BR,       B0,  B1,  B2;
      FLA_Obj WL,    WR,       W0,  W1,  W2;

      FLA_Part_1x2(~B,    &BL,  &BR,      0, FLA_LEFT );
      FLA_Part_1x2(~W,    &WL,  &WR,      0, FLA_LEFT );

      while ( FLA_Obj_width( BL ) < FLA_Obj_width(~B ) ){
        FLA_Repart_1x2_to_1x3( BL,   BR,        &B0,  &B1, &B2,
                               1, FLA_RIGHT );

        FLA_Repart_1x2_to_1x3( WL,   WR,        &W0,  &W1, &W2,
                               1, FLA_RIGHT );

        //------------------------------------------------------------
	Hier_ W_1(W1), B_1(B1);

#pragma omp task firstprivate( W_1, B_1 ) 
	FLA_Apply_Q_UT( FLA_LEFT, trans,
			direct, FLA_COLUMNWISE,
			~A_flat, ~T_flat, ~W_1.flat(), ~B_1.flat() );
                        
        //------------------------------------------------------------
        FLA_Cont_with_1x3_to_1x2( &BL,   &BR,        B0, B1,  B2,
                                  FLA_LEFT );

        FLA_Cont_with_1x3_to_1x2( &WL,   &WR,        W0, W1,  W2,
                                  FLA_LEFT );
      }

#pragma omp taskwait

      return true;
    }
    
    int apply_q2_l_var2( int trans, int direct, 
                            Hier_ D, Hier_ T, Hier_ W, Hier_ C, Hier_ E ) {
      if (T.get_m() && T.get_n()) {

        for (int i=0;i<T.get_m();++i) {

          for (int j=0;j<W.get_n();++j) {

#pragma omp task firstprivate(i,j)
            {
              FLA_Apply_Q2_UT_internal( FLA_LEFT, FLA_CONJ_TRANSPOSE,
                                        FLA_FORWARD, FLA_COLUMNWISE,
                                        D(i,0), T(i,0),
                                        W(0,j), C(0,j),
                                        E(i,j),
                                        fla_apq2ut_cntl_leaf );
            }
          }

#pragma omp taskwait

        }
      }
      return true;
    }
  }
}
