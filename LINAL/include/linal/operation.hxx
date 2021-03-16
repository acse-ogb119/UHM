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
#ifndef LINAL_OPERATION_HXX
#define LINAL_OPERATION_HXX

/*! 
  Multi threaded level 3 blas operations are listed.
*/ 

extern fla_qrutinc_t* flash_qrutinc_cntl;

extern fla_qrut_t   * fla_qrut_cntl_leaf;
extern fla_apqut_t  * fla_apqut_cntl_leaf;

extern fla_qr2ut_t  * fla_qr2ut_cntl_leaf;
extern fla_apq2ut_t * fla_apq2ut_cntl_leaf;

extern fla_axpy_t   * flash_axpy_cntl_blas;
extern fla_copy_t   * flash_copy_cntl_blas;
extern fla_scal_t   * flash_scal_cntl_blas;
extern fla_appiv_t  * flash_appiv_cntl;

extern fla_chol_t   * flash_chol_cntl_leaf;
extern fla_lu_t     * flash_lu_nopiv_cntl_leaf;
extern fla_lu_t     * flash_lu_piv_cntl_leaf;

extern fla_qrut_t   * flash_qrut_cntl_leaf;
extern fla_apqut_t  * flash_apqut_cntl_leaf;

extern fla_qr2ut_t  * flash_qr2ut_cntl_leaf;
extern fla_apq2ut_t * flash_apq2ut_cntl_leaf;

extern fla_trsm_t   * flash_trsm_cntl_blas;
extern fla_gemm_t   * flash_gemm_cntl_blas;
extern fla_trmm_t   * flash_trmm_cntl_blas;

extern fla_lu_t     * flash_lu_nopiv_cntl;
extern fla_trsm_t   * flash_trsm_cntl_mm;
extern fla_gemm_t   * flash_gemm_cntl_mm_op;
extern fla_trmm_t   * flash_trmm_cntl_mm;

/*!
  Macro definition
 */
#define LINAL_TR_LEFT(quad, side)		\
  {						\
    FLA_Part_2x2( ~A, &ATL, &ATR,		\
	              &ABL, &ABR, 0, 0, quad );	\
    FLA_Part_2x1( ~B, &BT,			\
		      &BB, 0, side);		\
  }

#define LINAL_TR_RIGHT(quad, side)		\
  {						\
    FLA_Part_2x2( ~A, &ATL, &ATR,		\
		      &ABL, &ABR, 0, 0, quad );	\
    FLA_Part_1x2( ~B, &BL, &BR, 0, side);	\
  }

#define LINAL_TR_LEFT_LOOP_HEAD(quad, side)			\
  {								\
    FLA_Repart_2x2_to_3x3( ATL, ATR,     &A00, &A01, &A02,	\
			   &A10, &A11, &A12,			\
                           ABL, ABR,     &A20, &A21, &A22,	\
                           1, 1, quad );			\
    FLA_Repart_2x1_to_3x1( BT,           &B0,			\
			   &B1,					\
                           BB,           &B2,			\
                           1, side);				\
  }

#define LINAL_TR_LEFT_LOOP_TAIL(quad, side)			\
  {                                                             \
    FLA_Cont_with_3x3_to_2x2( &ATL,  &ATR,  A00, A01, A02,      \
			                    A10, A11, A12,	\
                              &ABL, &ABR,   A20, A21, A22,      \
                              quad );                           \
    FLA_Cont_with_3x1_to_2x1( &BT,          B0,                 \
			                    B1,			\
                              &BB,          B2,                 \
                              side);                            \
  }

#define LINAL_TR_RIGHT_LOOP_HEAD(quad, side)			\
  {								\
    FLA_Repart_2x2_to_3x3( ATL, ATR,     &A00, &A01, &A02,	\
			                 &A10, &A11, &A12,	\
                           ABL, ABR,     &A20, &A21, &A22,	\
                           1, 1, quad );			\
    FLA_Repart_1x2_to_1x3( BL, BR,       &B0, &B1, &B2,		\
                           1, side);				\
  }

#define LINAL_TR_RIGHT_LOOP_TAIL(quad, side)			\
  {                                                             \
    FLA_Cont_with_3x3_to_2x2( &ATL,  &ATR,  A00, A01, A02,      \
			                    A10, A11, A12,	\
                              &ABL, &ABR,   A20, A21, A22,      \
                              quad );                           \
    FLA_Cont_with_1x3_to_1x2( &BL, &BR,     B0, B1, B2,         \
                              side );                           \
  }


namespace linal {
  // ----------------------------------------------------------------
  // ** Graphviz
  extern void head_graphviz ( FILE *fp, char *name );
  extern void tail_graphviz ( FILE *fp );
  extern int  push_graphviz ( Hier_ A, int begin,
                              std::vector< std::pair<int,double> > &viz );
  extern void write_graphviz( FILE *fp, char *task, char *style, char *color,
			      std::vector< std::pair<int,double> > &in,
                              std::vector< std::pair<int,double> > &out );
  extern int lu_nopiv ( FILE *fp, int &start, int &end, Hier_ A );

  // ----------------------------------------------------------------
  // ** Internal blas routines
  extern void gemm_cnt_disp();

  extern int gemm_internal( int transa, int transb,
                            FLA_Obj alpha, FLA_Obj A, FLA_Obj B,
                            FLA_Obj beta,  FLA_Obj C );
  extern int trsm_internal( int side, int uplo, int trans, int diag,
                            FLA_Obj alpha, FLA_Obj A, FLA_Obj B );
  
  extern int trmm_internal( int side, int uplo, int trans, int diag,
                            FLA_Obj alpha, FLA_Obj A, FLA_Obj B );

  // ----------------------------------------------------------------
  // ** Flat only
  extern double norm1     ( int uplo, int diag, Flat_ A );
  extern double inv_norm1 ( int uplo, int diag, Flat_ A );  
  extern double inv_norm1 ( double sample, int uplo, int diag, Flat_ A );  



  // ----------------------------------------------------------------
  // ** Hier only
  namespace dense {
    extern int apply_pivots( int side, int trans, Hier_ piv, Hier_ A );
    extern int apply_q     ( int side, int trans, int direct, int storev,
                             Hier_ A, Hier_ T, Hier_ W, Hier_ B ); 
    extern int apply_q2    ( int side, int trans, int direct, int storev,
                             Hier_ D, Hier_ T, Hier_ W, Hier_ C, Hier_ E );
    extern int apply_q_inc ( int side, int trans, int direct, int storev,
                             Hier_ A, Hier_ T, Hier_ W_1, Hier_ B );
    
    extern int chol        ( int uplo, Hier_ A );
    extern int lu_nopiv    ( Hier_ A );
    extern int lu_incpiv   ( Hier_ A, Hier_ p );
    extern int lu_piv      ( Hier_ A, Hier_ p );
    extern int qr          ( Hier_ A, Hier_ T ); 
    extern int qr_inc      ( Hier_ A, Hier_ T ); 
    extern int qr_inc_var1 ( Hier_ A, Hier_ T ); 
    extern int qr_var1     ( Hier_ A, Hier_ T ); 
    extern int qr_var2     ( Hier_ A, Hier_ T ); 
    extern int qr2         ( Hier_ U, Hier_ D, Hier_ T ); 

    extern int trmm        ( int side, int uplo, int trans, int diag,
                             FLA_Obj alpha, Hier_ A, Hier_ B );
    extern int trsm        ( int side, int uplo, int trans,
                             int diag, FLA_Obj alpha, Hier_ A, Hier_ B );
    extern int gemm        ( int trans, int transb, FLA_Obj alpha,
                             Hier_ A, Hier_ B, FLA_Obj beta, Hier_ C );
    extern int syrk        ( int uplo, int trans, FLA_Obj alpha, 
                             Hier_ A, FLA_Obj beta,  Hier_ C ); 
    extern int scal        ( FLA_Obj alpha, Hier_ A );
    
    // ----------------------------------------------------------------
    // ** Operation inside
    extern int apply_q_l_var2   ( int trans, int direct, 
                                  Hier_ A, Hier_ T, Hier_ W, Hier_ B );
    extern int apply_q2_l_var2  ( int trans, int direct, 
                                  Hier_ D, Hier_ T, Hier_ W, Hier_ C, Hier_ E );

    extern int apply_q_inc_lhfc_var1  ( Hier_ A, Hier_ T, Hier_ W_1, Hier_ B );
    extern int apply_q_inc_lnfc_var1  ( Hier_ A, Hier_ T, Hier_ W_1, Hier_ B );

    extern int l_t_upper_trsm   ( int diag, int trans, Hier_ A, Hier_ B );
    extern int l_t_lower_trsm   ( int diag, int trans, Hier_ A, Hier_ B );
    extern int l_nt_upper_trsm  ( int diag, int trans, Hier_ A, Hier_ B );
    extern int l_nt_lower_trsm  ( int diag, int trans, Hier_ A, Hier_ B );
    extern int r_t_upper_trsm   ( int diag, int trans, Hier_ A, Hier_ B );
    extern int r_t_lower_trsm   ( int diag, int trans, Hier_ A, Hier_ B );
    extern int r_nt_upper_trsm  ( int diag, int trans, Hier_ A, Hier_ B );
    extern int r_nt_lower_trsm  ( int diag, int trans, Hier_ A, Hier_ B );
    extern int trsm_update      ( int side, int uplo, int trans, int diag,
                                  Hier_ A, Hier_ B );

    extern int l_t_upper_trmm   ( int diag, int trans, FLA_Obj alpha, Hier_ A, Hier_ B );
    extern int l_t_lower_trmm   ( int diag, int trans, FLA_Obj alpha, Hier_ A, Hier_ B );
    extern int l_nt_upper_trmm  ( int diag, int trans, FLA_Obj alpha, Hier_ A, Hier_ B );
    extern int l_nt_lower_trmm  ( int diag, int trans, FLA_Obj alpha, Hier_ A, Hier_ B );
    extern int r_t_upper_trmm   ( int diag, int trans, FLA_Obj alpha, Hier_ A, Hier_ B );
    extern int r_t_lower_trmm   ( int diag, int trans, FLA_Obj alpha, Hier_ A, Hier_ B );
    extern int r_nt_upper_trmm  ( int diag, int trans, FLA_Obj alpha, Hier_ A, Hier_ B );
    extern int r_nt_lower_trmm  ( int diag, int trans, FLA_Obj alpha, Hier_ A, Hier_ B );
    extern int trmm_update      ( int side, int uplo, int trans, int diag,
                                  FLA_Obj alpha, Hier_ A, Hier_ B );
    
    extern int nt_syrk          ( int uplo, FLA_Obj alpha, Hier_ A, Hier_ C );
    extern int t_syrk           ( int uplo, FLA_Obj alpha, Hier_ A, Hier_ C );
  }
}
#endif
