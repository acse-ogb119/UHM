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
#include "uhm/util.hxx"

#include "uhm/matrix/uhm/matrix.hxx"
#include "uhm/matrix/uhm/fla.hxx"

namespace uhm {
  // --------------------------------------------------------------
  // ** Matrix
  extern linal::Flat_ nil_flat;
  extern linal::Hier_ nil_hier;

  extern double  buffer_used;
  extern double  max_buffer_used;
  extern double  flop;

  Matrix_FLA_::Matrix_FLA_() { /* shouldn't be called */  }
  Matrix_FLA_::Matrix_FLA_(int datatype, int fs, int ss, int n_rhs) {
    _init(datatype, fs, ss, n_rhs);
  }
  Matrix_FLA_::~Matrix_FLA_() { this->free(); }

  int Matrix_FLA_::is_created(int mat) { 
    return _get_flat(mat).is_created(); 
  }

  int Matrix_FLA_::is_buffer(int mat) { 
    return !_get_flat(mat).is_buffer_null(); 
  }

  int Matrix_FLA_::is_complex_datatype() {
    return ( this->datatype == UHM_COMPLEX );
  }

  std::pair<int,int> Matrix_FLA_::get_dimension() { 
    return (std::make_pair(this->fs, this->ss)); 
  }

  int Matrix_FLA_::get_n_rhs() { return this->n_rhs; }
  int Matrix_FLA_::get_compressed_mat_dim() { return this->cm; }

  void Matrix_FLA_::create_without_buffer() {
    int type = this->datatype;
    flat.ATL.create_without_buffer(type, this->fs, this->fs);
    flat.ATR.create_without_buffer(type, this->fs, this->ss);
    flat.ABL.create_without_buffer(type, this->ss, this->fs);
    flat.ABR.create_without_buffer(type, this->ss, this->ss);

    flat.p.create_without_buffer(FLA_INT, this->fs, 1);

    int b  = get_hier_block_size();

    flat.xt.create_without_buffer(type, this->fs, this->n_rhs);
    flat.xb.create_without_buffer(type, this->ss, this->n_rhs);
    flat.bt.create_without_buffer(type, this->fs, this->n_rhs);
    flat.bb.create_without_buffer(type, this->ss, this->n_rhs);
    flat.rt.create_without_buffer(type, this->fs, this->n_rhs);
    flat.rb.create_without_buffer(type, this->ss, this->n_rhs);
    
#ifdef UHM_HIER_MATRIX_ENABLE
    for (int i=UHM_ATL;i<UHM_T;++i) 
      if (is_created(i))
        _get_hier(i).create(_get_flat(i), b, b);

    for (int i=UHM_XT;i<UHM_END;++i) 
      if (is_created(i))
        _get_hier(i).create(_get_flat(i), b, b);
#endif
  }

  void Matrix_FLA_::free() {
#ifdef UHM_HIER_MATRIX_ENABLE
    for (int i=UHM_ATL;i<UHM_END;++i)
      _get_hier(i).free();
#endif

    for (int i=UHM_ATL;i<UHM_END;++i)
      _get_flat(i).free();
  }
  
  void Matrix_FLA_::copy_in(int mat, void *buffer) {
    linal::Flat_& obj = _get_flat(mat);
    FLA_Obj A;
    FLA_Obj_create_without_buffer( obj.get_data_type(),
				   obj.get_m(),
				   obj.get_n(),
				   &A );
    FLA_Obj_attach_buffer( buffer, 1, FLA_Obj_length(A), &A );
    FLA_Copy( A, ~obj );
    FLA_Obj_free_without_buffer( &A );
  }

  void Matrix_FLA_::copy_in(int mat, linal::Flat_ A) {
    FLA_Copy( ~A, ~(_get_flat(mat)) );
  }

  void Matrix_FLA_::copy_out(int mat, void *buffer) {
    linal::Flat_& obj = _get_flat(mat);
    FLA_Obj B;
    FLA_Obj_create_without_buffer( obj.get_data_type(),
				   obj.get_m(),
				   obj.get_n(),
				   &B );
    FLA_Obj_attach_buffer( buffer, 1, FLA_Obj_length(B), &B );
    FLA_Copy( ~obj, B );
    FLA_Obj_free_without_buffer( &B );
  }

  void Matrix_FLA_::copy_out(int mat, linal::Flat_ B) {
    FLA_Copy( ~(_get_flat(mat)), ~B );
  }

  void Matrix_FLA_::create_buffer() {
    for (int i=UHM_ATL;i<UHM_END;++i) 
      create_buffer(i);
  }

  void Matrix_FLA_::free_buffer() {
    for (int i=UHM_ATL;i<UHM_END;++i) 
      free_buffer(i);
  }

  void Matrix_FLA_::create_buffer(int mat) {
    linal::Flat_& obj = _get_flat(mat);
    if (obj.is_created() && obj.is_buffer_null()) 
      _create_buffer(obj);
  }

  void Matrix_FLA_::free_buffer(int mat) {
    _free_buffer(_get_flat(mat));
  }

  void Matrix_FLA_::check_solution() {
    // calculate | Ax -b |
    // do not calculate for schur complement
    if (this->fs) {
      FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE,
			~(this->flat.p), ~(this->flat.bt) );
      FLA_Axpy( FLA_MINUS_ONE, ~(this->flat.bt), ~(this->flat.rt) ) ;
    }
    if (this->ss) 
      FLA_Scal( FLA_ZERO, ~(this->flat.rb) );
  }

  void Matrix_FLA_::backup(int mat) {
    assert( mat == UHM_ATL || mat == UHM_XT  ||
            mat == UHM_BT  || mat == UHM_RT );
    linal::Flat_ obj = _get_flat(mat); 
    if (this->back.is_created()) this->back.free();
    this->back.create( obj.get_data_type(), 
		       obj.get_m(), obj.get_n() );
    FLA_Copy( ~obj, ~(this->back) );
  }

  void Matrix_FLA_::restore(int mat, int is_merge) {
    assert( mat == UHM_ATL || mat == UHM_XT  ||
            mat == UHM_BT  || mat == UHM_RT );
    linal::Flat_ obj = _get_flat(mat);
    if (!this->back.is_created()) return;
    
    // use should check whether the back is same object 
    // as he want to restore
    switch (is_merge) {
    case 0: FLA_Copy( ~(this->back), ~obj );break;
    case 1: FLA_Axpy( FLA_ONE, ~(this->back), ~obj );break;
    }
  }

  void Matrix_FLA_::apply_pivots( int mat ) {
    assert( mat == UHM_ATL || mat == UHM_XT  ||
	    mat == UHM_BT  || mat == UHM_RT );
    FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, 
		      ~(this->flat.p), ~(this->_get_flat(mat)) );
  }

  void Matrix_FLA_::set_zero( int mat ) {
    FLA_Obj_set_to_scalar( FLA_ZERO, ~(this->_get_flat(mat)) );
  }

  void Matrix_FLA_::improve_solution() {
    // Todo :: is not working yet
    // after solving (Ae = r), x += e
    if (this->fs)
      FLA_Axpy(FLA_MINUS_ONE, ~(this->flat.rt), ~(this->flat.xt) );
    if (this->ss)
      FLA_Axpy(FLA_MINUS_ONE, ~(this->flat.rb), ~(this->flat.xb) );
  }

  double Matrix_FLA_::get_residual() {
    if (!this->fs) return 0.0;

    double rval;
    linal::Flat_ norm;
    norm.create(LINAL_REAL, 1, 1);

    FLA_Norm1( ~(this->flat.rt), ~norm );
    rval = norm(0,0);

    norm.free();
    return rval;
  }

  double Matrix_FLA_::get_lower_triangular_norm() {
    if (!this->fs) return 0.0;

    double rval=0.0;
    linal::Flat_ sum_a, sum_b;
    sum_a.create(LINAL_REAL, 1, 1);
    sum_b.create(LINAL_REAL, 1, 1);

    for (int i=0;i<this->fs;++i) {
      linal::Flat_ tmp;
	
      // ATL :: diagonal members are included
      if (this->fs-i) {
	this->flat.ATL.extract( tmp, this->fs-i, 1, i, i );
	FLA_Norm1( ~tmp, ~sum_a );
      } else {
	sum_a(0,0) = 0.0;
      }
      
      // ABL
      if (this->ss) {
	this->flat.ABL.extract( tmp, this->ss, 1, 0, i );
	FLA_Norm1( ~tmp, ~sum_b );
      } else {
	sum_b(0,0) = 0.0;
      }

      rval = max(rval, sum_a(0,0)+sum_b(0,0));
    }

    sum_a.free();
    sum_b.free();

    return rval;
  }

  void Matrix_FLA_::random()             { _random(0); }
  void Matrix_FLA_::random_spd(int uplo) { _random(uplo); }
  void Matrix_FLA_::triangularize(int uplo) {
    // ATL and ABR are triangularized
    //
    if ( this->flat.ATL.is_created() && !this->flat.ATL.is_buffer_null() ) 
      FLA_Triangularize( uplo, FLA_NONUNIT_DIAG, ~(this->flat.ATL) );

    if ( this->flat.ABR.is_created() && !this->flat.ABR.is_buffer_null() )
      FLA_Triangularize( uplo, FLA_NONUNIT_DIAG, ~(this->flat.ABR) );
  }
  
  void Matrix_FLA_::merge(Matrix s, 
			  int mat_s, int offm_s, int offn_s,
			  int mat_t, int offm_t, int offn_t,
			  int m, int n, int is_erase) {
 
    if (!m || !n) return;

    // pointer casting into fla type
    Matrix_FLA src = (Matrix_FLA)s;
    assert(matrix_fla_valid(src));

    linal::Flat_ part_s, part_t;

    src->_get_flat(mat_s).extract(part_s, m, n, offm_s, offn_s);
    this->_get_flat(mat_t).extract(part_t, m, n, offm_t, offn_t);

    // merge
    FLA_Axpy(FLA_ONE, ~part_s, ~part_t);
    if (is_erase) FLA_Obj_set_to_scalar(FLA_ZERO, ~part_s);
  }

  void Matrix_FLA_::copy(Matrix s, 
			 int mat_s, int offm_s, int offn_s,
			 int mat_t, int offm_t, int offn_t,
			 int m, int n, int is_erase) {

    if (!m || !n) return;

    // pointer casting for fla 
    Matrix_FLA src = (Matrix_FLA)s;
    assert(matrix_fla_valid(src));
    
    linal::Flat_ part_s, part_t;
    
    src->_get_flat(mat_s).extract(part_s, m, n, offm_s, offn_s);
    this->_get_flat(mat_t).extract(part_t, m, n, offm_t, offn_t);

    FLA_Copy(~part_s, ~part_t);
    if (is_erase) FLA_Scal(FLA_ZERO, ~part_s);
  }

  void Matrix_FLA_::set_rhs(int is_leaf) {
    // prepare the rhs before solving sequence
    // when user want to solve system with different rhs
    // user can change rhs only then, set_rhs 
    // re-initialize the rhs living in upper hierarchy which 
    // is not leaf level : there could be gabage in upper hierarhcy level
    if (is_leaf) {
      // if leaf copy b to x
      if ( this->fs ) {
	assert(this->flat.xt.is_created() && this->flat.bt.is_created());
	FLA_Copy(~(this->flat.bt), ~(this->flat.xt));
      }
      if ( this->ss) {
	assert(this->flat.xb.is_created() && this->flat.bb.is_created());
	FLA_Copy(~(this->flat.bb), ~(this->flat.xb));
      }
    } else {
      // for non leaf, set to zero
      if ( this->fs && 
	   !this->flat.xt.is_buffer_null() &&
	   !this->flat.bt.is_buffer_null() ) {
	FLA_Obj_set_to_scalar( FLA_ZERO, ~(this->flat.xt) );
	FLA_Obj_set_to_scalar( FLA_ZERO, ~(this->flat.bt) );
      }
      if ( this->ss && 
	   !this->flat.xb.is_buffer_null() &&
	   !this->flat.bb.is_buffer_null() ) {
	FLA_Obj_set_to_scalar( FLA_ZERO, ~(this->flat.xb) );
	FLA_Obj_set_to_scalar( FLA_ZERO, ~(this->flat.bb) );
      }
    }
  }

  bool Matrix_FLA_::disp() { return this->disp(stdout); }
  bool Matrix_FLA_::disp(FILE* stream) {
    fprintf(stream, " - Matrix - \n");
    fprintf(stream, "   dimension [ %d , %d ], n_rhs [ %d ]\n", 
	   this->get_dimension().first, this->get_dimension().second,
	   this->get_n_rhs());
    //    this->flat.ATL.disp(stream, "- ATL -");
    this->flat.p.disp(stream, "- p -");
    this->flat.xt.disp(stream, "- xt -");
    //this->flat.xb.disp(stream, "- xb -");

    this->flat.bt.disp(stream, "- bt -");
    //this->flat.bb.disp(stream, "- bb -");

    this->flat.rt.disp(stream, "- rt -");
    //this->flat.rb.disp(stream, "- rb -");
    return true;
  }

  bool Matrix_FLA_::export_matrix(FILE* stream, int mat) {
    linal::Flat_& obj = this->_get_flat(mat);
    fprintf(stream, "### matrix\n");
    fprintf(stream, "%d\n", mat);
    fprintf(stream, "%d %d\n", obj.get_m(), obj.get_n());

    if (this->is_complex_datatype()) {
      for (int k2=0;k2<obj.get_n();++k2) {
	for (int k1=0;k1<obj.get_m();++k1) {
	  fprintf(stream, "% .14E % .14E\n", obj(k1,k2,0), obj(k1,k2,1));
	}
      }
    } else {
      for (int k2=0;k2<obj.get_n();++k2) {
	for (int k1=0;k1<obj.get_m();++k1) {
	  fprintf(stream, "% .14E\n", obj(k1,k2));
	}
      }
    }
    return true;
  }

  bool Matrix_FLA_::export_matrix(int &m, int &n,
				  std::vector< double > &val,
				  int mat) {
    linal::Flat_& obj = this->_get_flat(mat);
    //if (obj.is_buffer_null()) return true;

    m = obj.get_m(); 
    n = obj.get_n();

    if (this->is_complex_datatype()) {
      for (int k2=0;k2<obj.get_n();++k2) {
	for (int k1=0;k1<obj.get_m();++k1) {
	  val.push_back( obj(k1,k2,0) );
	  val.push_back( obj(k1,k2,1) );
	}
      }
    } else {
      for (int k2=0;k2<obj.get_n();++k2) {
	for (int k1=0;k1<obj.get_m();++k1) {
	  val.push_back( obj(k1,k2) );
	}
      }
    }
    return true;
  }

  bool Matrix_FLA_::import_matrix(int m, int n, int lda,
				  std::vector< double > &val,
				  int mat) {
    // skip if m or n is 0
    if (!m || !n) return true;

    linal::Flat_& obj = this->_get_flat(mat);
    assert( m == obj.get_m() &&
            n == obj.get_n() );

    if (this->is_complex_datatype()) 
      for (int k2=0;k2<obj.get_n();++k2) 
	for (int k1=0;k1<obj.get_m();++k1) {
	  obj(k1,k2,0) = val.at(lda*2*k2 + k1*2);
	  obj(k1,k2,1) = val.at(lda*2*k2 + k1*2+1);
	}
    else 
      for (int k2=0;k2<obj.get_n();++k2) 
	for (int k1=0;k1<obj.get_m();++k1) 
          obj(k1,k2) = val.at(lda*k2 + k1);
    
    //obj.disp("solution import_matrix");

    return true;
  }

  bool Matrix_FLA_::write_to_ooc(FILE* stream, int mat) {
    if (is_created(mat)) {
      linal::Flat_& obj = this->_get_flat(mat);
      long size = obj.get_buffer_size();
      if (size)
	write_buffer_to_file(stream, size, 
			     (char*)obj.get_buffer()); 
    }
    return true;
  }
  bool Matrix_FLA_::read_from_ooc(FILE* stream, int mat) {
    if (is_created(mat)) {
      linal::Flat_& obj = this->_get_flat(mat);
      long size = obj.get_buffer_size();
      if (size)
	read_buffer_from_file(stream, size,
			      (char*)obj.get_buffer()); 
    }
    return true;
  }

  // --------------------------------------------------------------
  // ** Protected
  void Matrix_FLA_::_init(int datatype, int fs, int ss, int n_rhs) {
    this->cookie     = UHM_MATRIX_FLA_COOKIE;
    this->fs         = fs;
    this->ss         = ss;
    this->n_rhs      = n_rhs;
    this->datatype   = datatype;
    this->cm         = fs;
  }

  void Matrix_FLA_::_create_buffer(linal::Matrix_ &obj) {
    if (obj.is_base_null()) return;
    if (!obj.is_hier()) {
      buffer_used     += obj.get_buffer_size();
      max_buffer_used  = max(max_buffer_used, buffer_used);
      obj.create_buffer();
    }
  }
  void Matrix_FLA_::_free_buffer(linal::Matrix_ &obj) {
    if (obj.is_base_null()) return;
    if (!obj.is_hier()) {
      buffer_used -= obj.get_buffer_size();
      obj.free_buffer();
    }
  }

  void Matrix_FLA_::_random(int uplo) {
    // Random matrix on A
    for (int i=UHM_ATL;i<UHM_P;++i) {
      linal::Flat_& obj = _get_flat(i);
      if ( obj.is_created() && !obj.is_buffer_null() ) {
	if (uplo && (i==UHM_ATL || i==UHM_ABR)) {
	  FLA_Random_spd_matrix( uplo, ~obj );
          FLA_Axpy(FLA_ONE, ~obj, ~obj);
        } else {
	  FLA_Random_matrix( ~obj );
        }
      }
    }

    // Random matrix on B
    for (int i=UHM_BT;i<UHM_RT;i++) {
      linal::Flat_& obj = _get_flat(i);
      if ( obj.is_created() && !obj.is_buffer_null() )
	FLA_Random_matrix( ~obj );
    }
    
    // set_rhs need for calculation
  }


  linal::Flat_& Matrix_FLA_::_get_flat(int mat) {
    switch (mat) {
    case UHM_ATL: return this->flat.ATL;break;
    case UHM_ATR: return this->flat.ATR;break;
    case UHM_ABL: return this->flat.ABL;break;
    case UHM_ABR: return this->flat.ABR;break;
    case UHM_P:   return this->flat.p;  break;
    case UHM_T:   return this->flat.T;  break;
    case UHM_XT:  return this->flat.xt; break;
    case UHM_XB:  return this->flat.xb; break;
    case UHM_BT:  return this->flat.bt; break;
    case UHM_BB:  return this->flat.bb; break;
    case UHM_RT:  return this->flat.rt; break;
    case UHM_RB:  return this->flat.rb; break;
    }
    return nil_flat;
  }

  linal::Hier_& Matrix_FLA_::_get_hier(int mat) {
    switch (mat) {
    case UHM_ATL: return this->hier.ATL;break;
    case UHM_ATR: return this->hier.ATR;break;
    case UHM_ABL: return this->hier.ABL;break;
    case UHM_ABR: return this->hier.ABR;break;
    case UHM_P:   return this->hier.p;  break;
    case UHM_T:   return this->hier.T;  break;
    case UHM_XT:  return this->hier.xt; break;
    case UHM_XB:  return this->hier.xb; break;
    case UHM_BT:  return this->hier.bt; break;
    case UHM_BB:  return this->hier.bb; break;
    case UHM_RT:  return this->hier.rt; break;
    case UHM_RB:  return this->hier.rb; break;
    }
    return nil_hier;
  }
}
