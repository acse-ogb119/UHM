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
  extern double flop;

  static inline int chol_flat( int fs, int ss, 
			       linal::Flat_ ATL, linal::Flat_ ATR,
			       linal::Flat_ ABL, linal::Flat_ ABR );

  static inline int chol_hier( int fs, int ss,
			       linal::Hier_ ATL, linal::Hier_ ATR,
			       linal::Hier_ ABL, linal::Hier_ ABR );
  
  void Matrix_FLA_::chol() {

#ifdef UHM_HIER_MATRIX_ENABLE
    // ----------------------------------------------------------
    // ** Hier-Matrix
    // ----------------------------------------------------------
    chol_hier( this->fs, this->ss, 
	       this->hier.ATL, this->hier.ATR,
	       this->hier.ABL, this->hier.ABR );
#else
    // ----------------------------------------------------------
    // ** Flat-Matrix
    // ----------------------------------------------------------
    chol_flat( this->fs, this->ss, 
	       this->flat.ATL, this->flat.ATR,
	       this->flat.ABL, this->flat.ABR );
#endif
  }

  static inline int chol_flat( int fs, int ss, 
			       linal::Flat_ ATL, linal::Flat_ ATR,
			       linal::Flat_ ABL, linal::Flat_ ABR ) {
    if (fs) { 
      FLA_Chol( FLA_LOWER_TRIANGULAR, ~ATL );
    }
  
    if (fs && ss) {
      FLA_Trsm( FLA_RIGHT, FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE,
		FLA_NONUNIT_DIAG, FLA_ONE, ~ATL, ~ABL );
      FLA_Syrk( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_MINUS_ONE,
		~ABL, FLA_ONE, ~ABR );
    }
    return true;
  }

  static inline int chol_hier( int fs, int ss, 
			       linal::Hier_ ATL, linal::Hier_ ATR,
			       linal::Hier_ ABL, linal::Hier_ ABR ) {
    if (fs) 
      linal::dense::chol( FLA_LOWER_TRIANGULAR, ATL);
  
    if (fs && ss) {
      linal::dense::trsm( FLA_RIGHT, FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE,
                          FLA_NONUNIT_DIAG, FLA_ONE, ATL, ABL );
      linal::dense::syrk( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_MINUS_ONE,
                          ABL, FLA_ONE, ABR );
    }
    return true;
  }
}
