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
    int apply_q(int side, int trans, int direct, int storev,
                Hier_ A, Hier_ T, Hier_ W, Hier_ B) {
      apply_q_l_var2( trans, direct, A, T, W, B );
      return true;
    }
    
    int apply_q2(int side, int trans, int direct, int storev,
                 Hier_ D, Hier_ T, Hier_ W, Hier_ C, Hier_ E) {
      
      apply_q2_l_var2( trans, direct, D, T, W, C, E );
      return true;
    }
    int apply_q_inc(int side, int trans, int direct, int storev,
                    Hier_ A, Hier_ T, Hier_ W_1, Hier_ B) {
      if (side == FLA_LEFT) {
        switch (trans) {
        case FLA_NO_TRANSPOSE: 
          apply_q_inc_lnfc_var1( A, T, W_1, B );
          break;
        case FLA_TRANSPOSE:
        case FLA_CONJ_TRANSPOSE:
          apply_q_inc_lhfc_var1( A, T, W_1, B );
          break;
        }
      } 
      return true;
    }
  }
}
