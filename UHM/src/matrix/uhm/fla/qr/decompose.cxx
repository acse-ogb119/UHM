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
  static int qr_flat( int fs, int ss,
                      linal::Flat_ ATL, linal::Flat_ ATR,
                      linal::Flat_ ABL, linal::Flat_ ABR,
                      linal::Flat_ T );
  
  static int qr_hier( int fs, int ss, 
                      linal::Hier_ ATL, linal::Hier_ ATR,
                      linal::Hier_ ABL, linal::Hier_ ABR,
                      linal::Hier_ T );

  void Matrix_FLA_::_qr_create_T() {
    int b  = get_hier_block_size();

#ifdef UHM_HIER_MATRIX_ENABLE

    // ** QR_INC only works for hier matrix
#ifdef UHM_QR_INC_ENABLE
    int t_fs = ( this->fs/b > 0 ? ( ( this->fs/b + (this->fs%b>0) )*b ) : this->fs );
    int t_ss = ( this->ss/b > 0 ? ( ( this->ss/b + (this->ss%b>0) )*b ) : this->ss );

    flat.T.create_without_buffer(this->datatype, t_fs, t_fs);
    hier.T.create(flat.T, b, b);
#else
    flat.T.create_without_buffer(this->datatype, this->fs, this->fs);
    hier.T.create(flat.T, b, b);
#endif


#else
    flat.T.create_without_buffer(this->datatype, this->fs, this->fs);
#endif

    this->create_buffer(UHM_T);
  }
  
  void Matrix_FLA_::qr() {

    // ** create matrix without buffer
    if (!is_created(UHM_T)) 
      this->_qr_create_T();

    // ** decomposition
#ifdef UHM_HIER_MATRIX_ENABLE
    // ----------------------------------------------------------
    // ** Hier-Matrix
    // ----------------------------------------------------------
    qr_hier( this->fs, this->ss, 
             this->hier.ATL, this->hier.ATR,
             this->hier.ABL, this->hier.ABR,
             this->hier.T );
#else
    // ----------------------------------------------------------
    // ** Flat-Matrix - Level Matrix
    // ----------------------------------------------------------
    qr_flat( this->fs, this->ss, 
             this->flat.ATL, this->flat.ATR,
             this->flat.ABL, this->flat.ABR,
             this->flat.T );
#endif
  }
  
  static inline int qr_flat( int fs, int ss, 
                             linal::Flat_ ATL, linal::Flat_ ATR,
                             linal::Flat_ ABL, linal::Flat_ ABR,
                             linal::Flat_ T ) {
    
    if (fs) 
      FLA_QR_UT( ~ATL, ~T );
    
    if (fs && ss) {

#pragma omp task firstprivate(ATL, ATR, T)
      {
        linal::Flat_ W;

        W.create    ( ATL.get_data_type(), T.get_m(), ATR.get_n() );

        FLA_Apply_Q_UT( FLA_LEFT, FLA_CONJ_TRANSPOSE,
                        FLA_FORWARD, FLA_COLUMNWISE,
                        ~ATL, ~T, ~W, ~ATR ); 
        W.free();
      }

#pragma omp task firstprivate(ATL, ABL)
      FLA_Trsm( FLA_RIGHT, FLA_UPPER_TRIANGULAR,
                FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
                FLA_ONE, ~ATL, ~ABL );

#pragma omp taskwait

      FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                FLA_MINUS_ONE, ~ABL, ~ATR, FLA_ONE, ~ABR );
    }

    return true;
  }

  static inline int qr_hier( int fs, int ss, 
                             linal::Hier_ ATL, linal::Hier_ ATR,
                             linal::Hier_ ABL, linal::Hier_ ABR,
                             linal::Hier_ T ) {
    if (fs) {
#ifdef UHM_QR_INC_ENABLE
      linal::dense::qr_inc( ATL, T );
#else
      linal::dense::qr( ATL, T );
#endif
    }

    if (fs && ss) {

#ifdef UHM_QR_INC_ENABLE

#pragma omp task firstprivate(ATL, ATR, T)
      {
	linal::Flat_ W;
	linal::Hier_ W_1;
	
	int mb   = get_hier_block_size();
        int t_fs = ( fs/mb > 0 ? ( ( fs/mb + (fs%mb>0) )*mb ) : fs );

	W.create    ( ATL.get_data_type(), fs, ATR.flat().get_n() );
	W_1.create  ( W, mb, mb );
	
	
	linal::dense::apply_q_inc( FLA_LEFT, FLA_CONJ_TRANSPOSE,
				   FLA_FORWARD, FLA_COLUMNWISE,
				   ATL, T, W_1, ATR ); 
	
        W_1.free();
	W.free();
      }

#else

#pragma omp task firstprivate(ATL, ATR, T)
      {
	linal::Flat_ W;
	linal::Hier_ W_1;
	
	int mb  = get_hier_block_size();
	
	W.create    ( ATL.get_data_type(), T.flat().get_m(), ATR.flat().get_n() );
	W_1.create  ( W, mb, mb );
	
	linal::dense::apply_q( FLA_LEFT, FLA_CONJ_TRANSPOSE,
                               FLA_FORWARD, FLA_COLUMNWISE,
                               ATL, T, W_1, ATR ); 
	
	W_1.free();
	W.free();
      }

#endif

#pragma omp task firstprivate(ATL, ABL)
      {
	linal::dense::trsm( FLA_RIGHT, FLA_UPPER_TRIANGULAR,
			    FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
			    FLA_ONE, ATL, ABL );
      }
      
#pragma omp taskwait
      
      linal::dense::gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
			  FLA_MINUS_ONE, ABL, ATR, FLA_ONE, ABR );
    }

#pragma omp taskwait 
    
    return true;
  }
}
