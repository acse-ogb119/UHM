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
      LU with incremental pivoting.
    */
    int lu_incpiv(Hier_ A, Hier_ p) {
      if (!A.get_m() || !A.get_n()) return true;

      LINAL_ERROR(A.get_m() == A.get_n(),
                  ">> Hier_ A is not square matrix");

      FLA_Obj ATL,   ATR,      A00,  A01,  A02,
              ABL,   ABR,      A10,  A11,  A12,
                              A20,  A21,  A22;
      FLA_Obj pt, pb, p0, p1, p2;

      FLA_Part_2x2( ~A, &ATL, &ATR,
                        &ABL, &ABR, 0, 0, FLA_TL );

      FLA_Part_2x1( ~p, &pt,
                        &pb, 0, FLA_TOP );

      while ( FLA_Obj_length( ATL ) < FLA_Obj_length( ~A ) ){
        FLA_Repart_2x2_to_3x3( ATL, ATR,     &A00, &A01, &A02,
                                             &A10, &A11, &A12,
                               ABL, ABR,     &A20, &A21, &A22,
                               1, 1, FLA_BR );
        FLA_Repart_2x1_to_3x1( pt, &p0,
                                   &p1,
                               pb, &p2, 1, FLA_BOTTOM );
        //------------------------------------------------------------
        Hier_ 
          A_10(A10), A_11(A11), A_12(A12), 
          A_21(A21), A_22(A22), p_1(p1);

        FLA_LU_piv( A_11(0,0) , p_1(0,0) );

#pragma omp task firstprivate( p_1, A_12 )
        apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE,
                      p_1, A_12 );

#pragma omp task firstprivate( p_1, A_10 )
        apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE,
                      p_1, A_10 );

#pragma omp taskwait

#pragma omp task firstprivate( A_11, A_12)
        trsm( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE,
              FLA_UNIT_DIAG, FLA_ONE, A_11, A_12 );

#pragma omp task firstprivate( A_11, A_21)
        trsm( FLA_RIGHT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE,
              FLA_NONUNIT_DIAG, FLA_ONE, A_11, A_21 );
      
#pragma omp taskwait
      
        gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_MINUS_ONE,
              A_21, A_12, FLA_ONE, A_22 );
      
        //------------------------------------------------------------
        FLA_Cont_with_3x3_to_2x2( &ATL,  &ATR,  A00, A01, A02,
                                                A10, A11, A12,
                                  &ABL, &ABR,   A20, A21, A22,
                                  FLA_TL );
        FLA_Cont_with_3x1_to_2x1( &pt, p0,
                                       p1,
                                  &pb, p2, FLA_TOP );
      }
      return true;
    }
  }
}
