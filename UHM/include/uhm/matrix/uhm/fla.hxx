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
#ifndef UHM_MATRIX_UHM_FLA_HXX
#define UHM_MATRIX_UHM_FLA_HXX

// flame control tree
extern fla_axpy_t * flash_axpy_cntl_blas;
extern fla_copy_t * flash_copy_cntl_blas;
extern fla_scal_t * flash_scal_cntl_blas;

extern fla_lu_t   * flash_lu_nopiv_cntl_leaf;
extern fla_trsm_t * flash_trsm_cntl_blas;
extern fla_gemm_t * flash_gemm_cntl_blas;
extern fla_trmm_t * flash_trmm_cntl_blas;

extern fla_lu_t   * flash_lu_nopiv_cntl;
extern fla_trsm_t * flash_trsm_cntl_mm;
extern fla_gemm_t * flash_gemm_cntl_mm_op;
extern fla_trmm_t * flash_trmm_cntl_mm;

namespace uhm {

  typedef class Matrix_FLA_*    Matrix_FLA;

  int matrix_fla_valid(Matrix_FLA m);

  // --------------------------------------------------------------
  // ** Unassembled Matrix_FLA class
  template<class M>
  struct Mat_FLA_ {
    M ATL, ATR, ABL, ABR,  p, T;
    M  xt,  xb,  bt,  bb, rt, rb;
  };
  
  class Matrix_FLA_ : public Matrix_ {
  private:
    int cookie;
  protected:
    int fs, ss, n_rhs, datatype, cm; 
    
    Mat_FLA_<linal::Flat_> flat;
    linal::Flat_& _get_flat( int mat );

    Mat_FLA_<linal::Hier_> hier;
    linal::Hier_& _get_hier( int mat );

    linal::Flat_ back;

    void _init                  ( int datatype, int fs, int ss, int n_rhs );
    void _create_buffer         ( linal::Matrix_ &obj );
    void _free_buffer           ( linal::Matrix_ &obj );
    void _random                ( int uplo );

    void _qr_create_T();

    // for checking

  public:
    Matrix_FLA_();
    Matrix_FLA_( int datatype, int fs, int ss, int n_rhs );
    virtual ~Matrix_FLA_();

    virtual bool disp();
    virtual bool disp( FILE *stream );

    virtual bool write_to_ooc(FILE* stream, int mat);
    virtual bool read_from_ooc(FILE* stream, int mat);

    virtual bool export_matrix( FILE *stream, int mat );
    virtual bool export_matrix( int &m, int &n,
				std::vector< double > &val,
				int mat );
    virtual bool import_matrix( int m, int n, int lda,
			        std::vector< double > &val,
			        int mat);

    virtual int is_created ( int mat );
    virtual int is_buffer  ( int mat );
    virtual int is_complex_datatype ();

    virtual std::pair<int,int> get_dimension();

    virtual int  get_n_rhs();
    virtual int  get_compressed_mat_dim();


    virtual void create_without_buffer();
    virtual void free();

    virtual void copy_in ( int mat, void *buffer );
    virtual void copy_out( int mat, void *buffer );

    virtual void copy_in ( int mat, linal::Flat_ A );
    virtual void copy_out( int mat, linal::Flat_ B );

    virtual void create_buffer();
    virtual void free_buffer();

    virtual void create_buffer ( int mat );
    virtual void free_buffer   ( int mat );

    virtual void merge( Matrix src, 
		        int mat_s, int offm_s, int offn_s,
		        int mat_t, int offm_t, int offn_t, 
		        int m, int n, int is_erase );

    virtual void copy ( Matrix src, 
			int mat_s, int offm_s, int offn_s,
			int mat_t, int offm_t, int offn_t, 
			int m, int n, int is_erase );


    virtual void   random();
    virtual void   random_spd( int uplo );
    virtual void   triangularize( int uplo );
    virtual void   set_rhs( int is_leaf );
    virtual double get_residual();
    virtual double get_lower_triangular_norm();
    virtual void   check_solution();
    virtual void   improve_solution();

    // for checking lu_piv
    virtual void backup( int mat );
    virtual void restore( int mat, int is_merge );
    virtual void apply_pivots( int mat );
    virtual void set_zero( int mat );
    // --------------------------------------------------------------
    virtual void chol();
    virtual void solve_chol_1_x();
    virtual void solve_chol_2_x();
    virtual void check_chol_1();
    virtual void check_chol_2();
    virtual void solve_chol_1_r();
    virtual void solve_chol_2_r();
    // --------------------------------------------------------------
    virtual void lu_nopiv();
    virtual void solve_lu_nopiv_1_x();
    virtual void solve_lu_nopiv_2_x();
    virtual void check_lu_nopiv_1();
    virtual void check_lu_nopiv_2();
    virtual void solve_lu_nopiv_1_r();
    virtual void solve_lu_nopiv_2_r();
    // --------------------------------------------------------------
    virtual void lu_incpiv();
    virtual void lu_piv();
    virtual void solve_lu_piv_1_x();
    virtual void solve_lu_piv_2_x();
    virtual void check_lu_piv_1();
    virtual void check_lu_piv_2();
    virtual void solve_lu_piv_1_r();
    virtual void solve_lu_piv_2_r();
    // --------------------------------------------------------------

    virtual void qr();
    virtual void solve_qr_1_x();
    virtual void solve_qr_2_x();
    virtual void check_qr_1();
    virtual void check_qr_2();
    virtual void solve_qr_1_r();
    virtual void solve_qr_2_r();
    // --------------------------------------------------------------
   
    // friend interface for new object
    friend inline int matrix_fla_valid(Matrix_FLA m);
  };
  // --------------------------------------------------------------
  // ** Definition
  inline int matrix_fla_valid(Matrix_FLA m) {
    return (m && m->cookie == UHM_MATRIX_FLA_COOKIE);
  }
}

#endif
