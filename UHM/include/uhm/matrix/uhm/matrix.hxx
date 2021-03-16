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
#ifndef UHM_MATRIX_UHM_MATRIX_HXX
#define UHM_MATRIX_UHM_MATRIX_HXX

namespace uhm {

  typedef class Matrix_*    Matrix;

  extern void   set_hier_block_size(int size);
  extern int    get_hier_block_size();

  extern double matrix_buffer_used();
  extern double matrix_max_buffer_used();
  extern double matrix_flop();
  extern void   matrix_reset_flop();
  extern void   matrix_reset_buffer();

  // --------------------------------------------------------------
  // ** Abstract class for the interface 

  // connected to element object
  class Matrix_ {
  public:
    Matrix_() { }
    Matrix_( int datatype, int fs, int ss, int n_rhs ) { }
    //Matrix_( int datatype, int fs, int ss, int n_rhs,
    //         MPI_Comm comm ) { }
    virtual ~Matrix_() { }
    virtual bool disp()=0;
    virtual bool disp( FILE *stream )=0;

    virtual bool write_to_ooc(FILE* stream, int mat)=0;
    virtual bool read_from_ooc(FILE* stream, int mat)=0;

    virtual bool export_matrix( FILE *stream, int mat )=0;
    virtual bool export_matrix( int &m, int &n,
			        std::vector< double > &val,
			        int mat)=0;
    virtual bool import_matrix( int m, int n, int lda,
			        std::vector< double > &val,
			        int mat)=0;

    virtual int  is_created( int mat )=0;
    virtual int  is_buffer ( int mat )=0;
    virtual int  is_complex_datatype()=0;

    virtual std::pair<int,int> get_dimension()=0;

    virtual int  get_n_rhs()=0;
    virtual int  get_compressed_mat_dim()=0;

    virtual void create_without_buffer()=0;
    virtual void free()=0;

    virtual void copy_in ( int mat, void *buffer )=0;
    virtual void copy_out( int mat, void *buffer )=0;

    virtual void copy_in ( int mat, linal::Flat_ A )=0;
    virtual void copy_out( int mat, linal::Flat_ B )=0;

    virtual void create_buffer()=0;
    virtual void free_buffer()=0;

    virtual void create_buffer( int mat )=0;
    virtual void free_buffer  ( int mat )=0;

    virtual void merge( Matrix src, 
			int mat_s, int offm_s, int offn_s,
			int mat_t, int offm_t, int offn_t,
			int m, int n, int is_erase )=0;
    
    virtual void copy ( Matrix src, 
			int mat_s, int offm_s, int offn_s,
			int mat_t, int offm_t, int offn_t,
			int m, int n, int is_erase )=0;

    virtual void   random()=0;
    virtual void   random_spd( int uplo )=0;
    virtual void   triangularize( int uplo )=0;
    virtual void   set_rhs( int is_leaf )=0;
    virtual double get_residual()=0;
    virtual double get_lower_triangular_norm()=0;
    virtual void   check_solution()=0;
    virtual void   improve_solution()=0;

    // for checkin routine internally need
    virtual void backup( int mat )=0;
    virtual void restore( int mat, int is_merge )=0;
    virtual void apply_pivots( int mat )=0;
    virtual void set_zero( int mat )=0;
    // --------------------------------------------------------------
    virtual void chol()=0;
    virtual void solve_chol_1_x()=0;
    virtual void solve_chol_2_x()=0;
    virtual void check_chol_1()=0;
    virtual void check_chol_2()=0;
    virtual void solve_chol_1_r()=0;
    virtual void solve_chol_2_r()=0;
    // --------------------------------------------------------------
    virtual void lu_nopiv()=0;
    virtual void solve_lu_nopiv_1_x()=0;
    virtual void solve_lu_nopiv_2_x()=0;
    virtual void check_lu_nopiv_1()=0;
    virtual void check_lu_nopiv_2()=0;
    virtual void solve_lu_nopiv_1_r()=0;
    virtual void solve_lu_nopiv_2_r()=0;
    // --------------------------------------------------------------
    virtual void lu_incpiv()=0;
    virtual void lu_piv()=0;
    virtual void solve_lu_piv_1_x()=0;
    virtual void solve_lu_piv_2_x()=0;
    virtual void check_lu_piv_1()=0;
    virtual void check_lu_piv_2()=0;
    virtual void solve_lu_piv_1_r()=0;
    virtual void solve_lu_piv_2_r()=0;
    // --------------------------------------------------------------
    virtual void qr()=0;
    virtual void solve_qr_1_x()=0;
    virtual void solve_qr_2_x()=0;
    virtual void check_qr_1()=0;
    virtual void check_qr_2()=0;
    virtual void solve_qr_1_r()=0;
    virtual void solve_qr_2_r()=0;
    // --------------------------------------------------------------
  };
}

#endif
