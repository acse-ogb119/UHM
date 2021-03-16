/*
  Copyright © 2011, Kyungjoo Kim
  All rights reserved.
  
  This file is part of UHM.
  
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
#include "uhm/common.hxx"
#include "uhm/const.hxx"

#include "uhm/matrix/uhm/matrix.hxx"
#include "uhm/matrix/uhm/fla.hxx"

namespace uhm {
  void Matrix_FLA_::check_chol_1() {
#ifdef UHM_HIER_MATRIX_ENABLE
    // ----------------------------------------------------------
    // ** Hier-Matrix 
    // ----------------------------------------------------------
    if (this->fs) {
      FLA_Copy( ~(this->flat.xt), ~(this->flat.rt) );
      linal::dense::trmm( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_CONJ_TRANSPOSE,
                          FLA_NONUNIT_DIAG, FLA_ONE, 
                          this->hier.ATL, this->hier.rt ); 
    }

    if (this->fs && this->ss) {
      linal::dense::gemm( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE,
                          this->hier.ABL, this->hier.xb, 
                          FLA_ONE, this->hier.rt );
    }
#else
    // ----------------------------------------------------------
    // ** Flat-Matrix 
    // ----------------------------------------------------------
    if (this->fs) {
      FLA_Copy( ~(this->flat.xt), ~(this->flat.rt) );
      FLA_Trmm( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_CONJ_TRANSPOSE,
		FLA_NONUNIT_DIAG, FLA_ONE, 
		~(this->flat.ATL), ~(this->flat.rt) ); 
    }

    if (this->fs && this->ss) {
      FLA_Gemm( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE,
		~(this->flat.ABL), ~(this->flat.xb), 
		FLA_ONE, ~(this->flat.rt) );
    }
#endif
  }
   
  // from leaf to root
  void Matrix_FLA_::check_chol_2() {
#ifdef UHM_HIER_MATRIX_ENABLE
    // ----------------------------------------------------------
    // ** Hier-Matrix 
    // ----------------------------------------------------------
    if (this->fs && this->ss) {
      linal::dense::gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE,
                          this->hier.ABL, this->hier.rt, 
                          FLA_ONE, this->hier.rb );
    }
    
    if (this->fs) {
      linal::dense::trmm( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE,
                          FLA_NONUNIT_DIAG, FLA_ONE, 
                          this->hier.ATL, this->hier.rt ); 
    }
    // rb should be merged for upper hierarchy
#else
    // ----------------------------------------------------------
    // ** Flat-Matrix 
    // ----------------------------------------------------------
    if (this->fs && this->ss) {
      FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE,
		~(this->flat.ABL), ~(this->flat.rt), 
		FLA_ONE, ~(this->flat.rb) );
    }
    
    if (this->fs) {
      FLA_Trmm( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE,
		FLA_NONUNIT_DIAG, FLA_ONE, 
		~(this->flat.ATL), ~(this->flat.rt) ); 
    }
    // rb should be merged for upper hierarchy
#endif
  }
}
