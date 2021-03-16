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

// BLK variant selection
#undef LINAL_TRMM_BLK_VAR1
#undef LINAL_TRMM_BLK_VAR2

#define LINAL_TRMM_BLK_VAR1


namespace linal {
  namespace dense {
    
    int trmm( int side, int uplo, int trans, int diag,
              FLA_Obj alpha, Hier_ A, Hier_ B ) {
      if (!A.get_m() || !A.get_n() || !B.get_m() || !B.get_n() ) return true;

      int val=0;
      if (uplo == FLA_UPPER_TRIANGULAR) {
        val = side*1000+trans;
        switch (val) {
        case 210400: // left + no trans + upper
          l_nt_upper_trmm( diag, trans, alpha, A, B );break;
        case 210401: // left + trans + upper
        case 210402:
          l_t_upper_trmm ( diag, trans, alpha, A, B );break;
        case 211400: // right + no trans + uppper
          r_nt_upper_trmm( diag, trans, alpha, A, B );break;
        case 211401: // right + trans + upper
        case 211402:
          r_t_upper_trmm ( diag, trans, alpha, A, B);break;
        default: val = 0; break;
        }
      } else {
        val = side*1000+trans;
        switch (val) {
        case 210400: // left + no trans + lower
          l_nt_lower_trmm(diag, trans, alpha, A, B);break;
        case 210401: // left + trans + lower
        case 210402:
          l_t_lower_trmm (diag, trans, alpha, A, B);break;
        case 211400: // right + no trans + lower
          r_nt_lower_trmm(diag, trans, alpha, A, B);break;
        case 211401: // right + trans + lower
        case 211402:
          r_t_lower_trmm (diag, trans, alpha, A, B);break;
        default: val = 0; break;
        }
      }
      assert(val);

      return true;
    }

    int l_nt_upper_trmm( int diag, int trans,
                         FLA_Obj alpha, Hier_ A, Hier_ B) {
      FLA_Obj ATL, ATR, ABL, ABR,
        A00, A01,  A02, A10, A11, A12, A20, A21, A22,
        BT, BB, BL, BR,
        B0, B1, B2;

      LINAL_TR_LEFT(FLA_TL, FLA_TOP);

      while ( FLA_Obj_length( ATL ) < FLA_Obj_length( ~A ) ){
        LINAL_TR_LEFT_LOOP_HEAD(FLA_BR, FLA_BOTTOM);
        //------------------------------------------------------------
#ifdef LINAL_TRMM_BLK_VAR1
        {
          trmm_update( FLA_LEFT,
                       FLA_UPPER_TRIANGULAR,
                       FLA_NO_TRANSPOSE,
                       diag, alpha,
                       Hier_(A11), Hier_(B1) );
          
          gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                alpha, Hier_(A12), Hier_(B2),
                FLA_ONE, Hier_(B1) );
        }
#endif

#ifdef LINAL_TRMM_BLK_VAR2
        {
          gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                alpha, Hier_(A01), Hier_(B1),
                FLA_ONE, Hier_(B0) );

          trmm_update( FLA_LEFT,
                       FLA_UPPER_TRIANGULAR,
                       FLA_NO_TRANSPOSE,
                       diag, alpha,
                       Hier_(A11), Hier_(B1) );
        }
#endif

        //------------------------------------------------------------
        LINAL_TR_LEFT_LOOP_TAIL(FLA_TL, FLA_TOP);
      }
      return true;
    }
  
    int l_nt_lower_trmm( int diag, int trans, 
                         FLA_Obj alpha, Hier_ A, Hier_ B) {
      FLA_Obj ATL, ATR, ABL, ABR,
        A00, A01,  A02, A10, A11, A12, A20, A21, A22,
        BT, BB, BL, BR,
        B0, B1, B2;

      LINAL_TR_LEFT(FLA_BR, FLA_BOTTOM);

      while ( FLA_Obj_length( ABR ) < FLA_Obj_length( ~A ) ){
        LINAL_TR_LEFT_LOOP_HEAD(FLA_TL, FLA_TOP);
        //------------------------------------------------------------

#ifdef LINAL_TRMM_BLK_VAR1	
 	{ // block variant 1
 	  trmm_update( FLA_LEFT,
 		       FLA_LOWER_TRIANGULAR,
 		       FLA_NO_TRANSPOSE,
 		       diag, alpha,
 		       Hier_(A11), Hier_(B1) );
	  
 	  gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
 		alpha, Hier_(A10), Hier_(B0),
 		FLA_ONE, Hier_(B1) );
 	}
#endif

#ifdef LINAL_TRMM_BLK_VAR2
	{ // block variant 2
	  gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
		alpha, Hier_(A21), Hier_(B1),
		FLA_ONE, Hier_(B2) );

	  trmm_update( FLA_LEFT,
		       FLA_LOWER_TRIANGULAR,
		       FLA_NO_TRANSPOSE,
		       diag, alpha,
		       Hier_(A11), Hier_(B1) );
	}
#endif
        //------------------------------------------------------------
        LINAL_TR_LEFT_LOOP_TAIL(FLA_BR, FLA_BOTTOM);
      }

      return true;
    }
    int l_t_upper_trmm( int diag, int trans,
                        FLA_Obj alpha, Hier_ A, Hier_ B ) {
      FLA_Obj ATL, ATR, ABL, ABR,
        A00, A01,  A02, A10, A11, A12, A20, A21, A22,
        BT, BB, BL, BR,
        B0, B1, B2;

      LINAL_TR_LEFT(FLA_BR, FLA_BOTTOM);

      while ( FLA_Obj_length( ABR ) < FLA_Obj_length( ~A ) ){
        LINAL_TR_LEFT_LOOP_HEAD(FLA_TL, FLA_TOP);
        //------------------------------------------------------------
#ifdef LINAL_TRMM_BLK_VAR1	
        {
          trmm_update( FLA_LEFT,
                       FLA_UPPER_TRIANGULAR,
                       trans,
                       diag, alpha,
                       Hier_(A11), Hier_(B1) );
          
          gemm( trans, FLA_NO_TRANSPOSE,
                alpha, Hier_(A01), Hier_(B0),
                FLA_ONE, Hier_(B1) );
        }
#endif

#ifdef LINAL_TRMM_BLK_VAR2
        {
          gemm( trans, FLA_NO_TRANSPOSE,
                alpha, Hier_(A12), Hier_(B1),
                FLA_ONE, Hier_(B2) );

          trmm_update( FLA_LEFT,
                       FLA_UPPER_TRIANGULAR,
                       trans,
                       diag, alpha,
                       Hier_(A11), Hier_(B1) );
        }
#endif

        //------------------------------------------------------------
        LINAL_TR_LEFT_LOOP_TAIL(FLA_BR, FLA_BOTTOM);
      }

      return true;
    }
    int l_t_lower_trmm( int diag, int trans, 
                        FLA_Obj alpha, Hier_ A, Hier_ B ) {
      FLA_Obj ATL, ATR, ABL, ABR,
        A00, A01,  A02, A10, A11, A12, A20, A21, A22,
        BT, BB, BL, BR,
        B0, B1, B2;

      LINAL_TR_LEFT(FLA_TL, FLA_TOP);

      while ( FLA_Obj_length( ATL ) < FLA_Obj_length( ~A ) ){
        LINAL_TR_LEFT_LOOP_HEAD(FLA_BR, FLA_BOTTOM);
        //------------------------------------------------------------
#ifdef LINAL_TRMM_BLK_VAR1	
        {
          trmm_update( FLA_LEFT,
                       FLA_LOWER_TRIANGULAR,
                       trans,
                       diag, alpha,
                       Hier_(A11), Hier_(B1) );

          gemm( trans, FLA_NO_TRANSPOSE,
                alpha, Hier_(A21), Hier_(B2),
                FLA_ONE, Hier_(B1) );
        }
#endif

#ifdef LINAL_TRMM_BLK_VAR2
        {
          gemm( trans, FLA_NO_TRANSPOSE,
                alpha, Hier_(A10), Hier_(B1),
                FLA_ONE, Hier_(B0) );

          trmm_update( FLA_LEFT,
                       FLA_LOWER_TRIANGULAR,
                       trans,
                       diag, alpha,
                       Hier_(A11), Hier_(B1) );
        }
#endif

        //------------------------------------------------------------
        LINAL_TR_LEFT_LOOP_TAIL(FLA_TL, FLA_TOP);
      }

      return true;
    }

    int r_t_upper_trmm( int diag, int trans,
                        FLA_Obj alpha, Hier_ A, Hier_ B ) {
      FLA_Obj ATL, ATR, ABL, ABR,
        A00, A01,  A02, A10, A11, A12, A20, A21, A22,
        BT, BB, BL, BR,
        B0, B1, B2;

      LINAL_TR_RIGHT(FLA_TL, FLA_LEFT);

      while ( FLA_Obj_length( ATL ) < FLA_Obj_length( ~A ) ){
        LINAL_TR_RIGHT_LOOP_HEAD(FLA_BR, FLA_RIGHT);
        //------------------------------------------------------------
#ifdef LINAL_TRMM_BLK_VAR1
        {
          trmm_update( FLA_RIGHT,
                       FLA_UPPER_TRIANGULAR,
                       trans,
                       diag, alpha,
                       Hier_(A11), Hier_(B1) );
          
          gemm( FLA_NO_TRANSPOSE, trans,
                alpha, Hier_(B2), Hier_(A12),
                FLA_ONE, Hier_(B1) );
        }
#endif

#ifdef LINAL_TRMM_BLK_VAR2
        {
          gemm( FLA_NO_TRANSPOSE, trans,
                alpha, Hier_(B1), Hier_(A01),
                FLA_ONE, Hier_(B0) );

          trmm_update( FLA_RIGHT,
                       FLA_UPPER_TRIANGULAR,
                       trans,
                       diag, alpha,
                       Hier_(A11), Hier_(B1) );
        }
#endif
        

        //------------------------------------------------------------
        LINAL_TR_RIGHT_LOOP_TAIL(FLA_TL, FLA_LEFT);
      }

      return true;
    }
    int r_t_lower_trmm( int diag, int trans,
                        FLA_Obj alpha, Hier_ A, Hier_ B) {
      FLA_Obj ATL, ATR, ABL, ABR,
        A00, A01,  A02, A10, A11, A12, A20, A21, A22,
        BT, BB, BL, BR,
        B0, B1, B2;

      LINAL_TR_RIGHT(FLA_BR, FLA_RIGHT);

      while ( FLA_Obj_length( ABR ) < FLA_Obj_length( ~A ) ){
        LINAL_TR_RIGHT_LOOP_HEAD(FLA_TL, FLA_LEFT);
        //------------------------------------------------------------
#ifdef LINAL_TRMM_BLK_VAR1
        {
          trmm_update( FLA_RIGHT,
                       FLA_LOWER_TRIANGULAR,
                       trans,
                       diag, alpha,
                       Hier_(A11), Hier_(B1) );
          
          gemm( FLA_NO_TRANSPOSE, trans,
                alpha, Hier_(B0), Hier_(A10),
                FLA_ONE, Hier_(B1) );
        }
#endif

#ifdef LINAL_TRMM_BLK_VAR2
        { 
          gemm( FLA_NO_TRANSPOSE, trans,
                alpha, Hier_(B1), Hier_(A21),
                FLA_ONE, Hier_(B2) );
       
          trmm_update( FLA_RIGHT,
                       FLA_LOWER_TRIANGULAR,
                       trans,
                       diag, alpha,
                       Hier_(A11), Hier_(B1) );
        }
#endif

        //------------------------------------------------------------
        LINAL_TR_RIGHT_LOOP_TAIL(FLA_BR, FLA_RIGHT);
      }

      return true;
    }
    int r_nt_upper_trmm( int diag, int trans, 
                         FLA_Obj alpha, Hier_ A, Hier_ B ) {
      FLA_Obj ATL, ATR, ABL, ABR,
        A00, A01,  A02, A10, A11, A12, A20, A21, A22,
        BT, BB, BL, BR,
        B0, B1, B2;

      LINAL_TR_RIGHT(FLA_BR, FLA_RIGHT);

      while ( FLA_Obj_length( ABR ) < FLA_Obj_length( ~A ) ){
        LINAL_TR_RIGHT_LOOP_HEAD(FLA_TL, FLA_LEFT);
        //------------------------------------------------------------
#ifdef LINAL_TRMM_BLK_VAR1
        {
          trmm_update( FLA_RIGHT,
                       FLA_UPPER_TRIANGULAR,
                       FLA_NO_TRANSPOSE,
                       diag, alpha,
                       Hier_(A11), Hier_(B1) );
          
          gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                alpha, Hier_(B0), Hier_(A01),
                FLA_ONE, Hier_(B1) );
        }
#endif

#ifdef LINAL_TRMM_BLK_VAR2
        {
          gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                alpha, Hier_(B1), Hier_(A12),
                FLA_ONE, Hier_(B2) );

          trmm_update( FLA_RIGHT,
                       FLA_UPPER_TRIANGULAR,
                       FLA_NO_TRANSPOSE,
                       diag, alpha,
                       Hier_(A11), Hier_(B1) );
        }
#endif

        //------------------------------------------------------------
        LINAL_TR_RIGHT_LOOP_TAIL(FLA_BR, FLA_RIGHT);
      }

      return true;
    }

    int r_nt_lower_trmm( int diag, int trans,
                         FLA_Obj alpha, Hier_ A, Hier_ B ) {
      FLA_Obj ATL, ATR, ABL, ABR,
        A00, A01,  A02, A10, A11, A12, A20, A21, A22,
        BT, BB, BL, BR,
        B0, B1, B2;

      LINAL_TR_RIGHT(FLA_TL, FLA_LEFT);

      while ( FLA_Obj_length( ATL ) < FLA_Obj_length( ~A ) ){
        LINAL_TR_RIGHT_LOOP_HEAD(FLA_BR, FLA_RIGHT);
        //------------------------------------------------------------
#ifdef LINAL_TRMM_BLK_VAR1
        {
          trmm_update( FLA_RIGHT,
                       FLA_LOWER_TRIANGULAR,
                       FLA_NO_TRANSPOSE,
                       diag, alpha,
                       Hier_(A11), Hier_(B1) );
          
          gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                alpha, Hier_(B2), Hier_(A21),
                FLA_ONE, Hier_(B1) );
        }
#endif

#ifdef LINAL_TRMM_BLK_VAR2
        {
          gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                alpha, Hier_(B1), Hier_(A10),
                FLA_ONE, Hier_(B0) );

          trmm_update( FLA_RIGHT,
                       FLA_LOWER_TRIANGULAR,
                       FLA_NO_TRANSPOSE,
                       diag, alpha,
                       Hier_(A11), Hier_(B1) );
        }
#endif


        //------------------------------------------------------------
        LINAL_TR_RIGHT_LOOP_TAIL(FLA_TL, FLA_LEFT);
      }

      return true;
    }

    int trmm_update( int side, int uplo, int trans, int diag,
                     FLA_Obj alpha,
                     Hier_ A, Hier_ B ) {
      // Given input A is one 1x1 block matrix
      // B is column or row block matrix 
      // it apply TRMM on rowise or columnwise for TRMM block

      for (int k2=0;k2<B.get_n();++k2) {
        for (int k1=0;k1<B.get_m();++k1) {

#pragma omp task firstprivate( k1, k2 )
          trmm_internal( side, uplo, trans, diag, alpha, A(0,0), B(k1,k2) );

        }
      }

#pragma omp taskwait

      return true;
    }
  }
}
