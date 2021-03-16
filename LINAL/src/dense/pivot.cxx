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
      Apply permutation of p to matrix A. 
    */
    int apply_pivots(int side, int trans, Hier_ p, Hier_ A) {
      LINAL_ERROR( p.get_data_type() == LINAL_INT,
                   ">> Pivot matrix is not INT type");
    
      FLA_Obj AT, AB, A0, A1, A2;
      FLA_Obj pt, pb, p0, p1, p2;

      FLA_Part_2x1( ~A, &AT,
                        &AB, 0, FLA_TOP );
      FLA_Part_2x1( ~p, &pt,
                        &pb, 0, FLA_TOP );

      while ( FLA_Obj_length( pt ) < FLA_Obj_length( ~p ) ) {
        FLA_Repart_2x1_to_3x1( AT, &A0,
                                   &A1,
                               AB, &A2, 1, FLA_BOTTOM );
        FLA_Repart_2x1_to_3x1( pt, &p0,
                                   &p1,
                               pb, &p2, 1, FLA_BOTTOM );
        //------------------------------------------------------------
        FLA_Obj AA;
        FLA_Merge_2x1( A1,
                       A2, &AA );

        Hier_ HH(AA), p_1(p1);
        for (int i=0;i<HH.get_n();++i) {
          Hier_ T;
          HH.extract( T, HH.get_m(), 1, 0, i );

#pragma omp task firstprivate( p_1, T )
          FLA_Apply_pivots( side, trans, p_1(0,0), ~(T.flat()) );

        }

#pragma omp taskwait

        //------------------------------------------------------------
        FLA_Cont_with_3x1_to_2x1( &AT, A0,
                                       A1,
                                  &AB, A2, FLA_TOP );
        FLA_Cont_with_3x1_to_2x1( &pt, p0,
                                       p1,
                                  &pb, p2, FLA_TOP );
      }
      return true;
    }
  }
}

