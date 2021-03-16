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
#ifndef UHM_MATRIX_UHM_HELPER_HXX
#define UHM_MATRIX_UHM_HELPER_HXX

namespace uhm {
  // --------------------------------------------------------------
  typedef class  Helper_* Helper;
  typedef class  Element_* Element;
  typedef struct Mapper_* Mapper;
  typedef struct Mapper_* Mapper;

  struct Mapper_ { int offs_c, offs_p, fs_p, n_dof; };

  // ** Helper matrix
  class Helper_ {
  protected:
    int     cookie;
    Element p, c;

    std::vector<Mapper_> mapper;

    void _merge_A    ();
    void _branch_ABR ();

    void _merge_rhs  (int kind, int is_pivot_applied);
    void _branch_rhs (int kind);

    int  _get_rhs    (int col);
    int  _get_A      (int col, int row);

  public:
    Helper_();
    Helper_(Element p, Element c);
    virtual ~Helper_();

    bool disp();
    bool disp(FILE *stream);

    void set_mapper();

    void merge_A();
    void branch_ABR();

    void merge_rhs_x();
    void merge_rhs_b();
    void merge_rhs_r();

    void branch_rhs_x();
    void branch_rhs_b();
    void branch_rhs_r();
  };

  // --------------------------------------------------------------
  // ** Definition
  inline Helper_::Helper_() {  }
  inline Helper_::Helper_(Element p, Element c) {
    assert( element_valid(p)       && element_valid(c) &&
	    p->is_matrix_created() && c->is_matrix_created() );

    this->cookie = UHM_HELPER_COOKIE;
    this->p      = p;    
    this->c      = c;
  }
  inline Helper_::~Helper_() { }

  // --------------------------------------------------------------
  inline void Helper_::set_mapper() {
    Mapper_ q;

    std::map< Node, int > factor, schur;

    // ** dump vector into map
    for (int i=0;i<this->p->factor.size();++i) 
      factor.insert( this->p->factor.at(i) );
    for (int i=0;i<this->p->schur.size();++i)
      schur.insert( this->p->schur.at(i) );

    std::vector<Mapper_> bijection;
    bijection.reserve(this->c->schur.size());

    // ** visit child schur nodes
    std::map< Node, int>::iterator pit;
    std::vector< std::pair<Node, int> >::iterator cit;
    for (cit=this->c->schur.begin();cit!=this->c->schur.end();++cit) {
      q.offs_c = cit->second;
      
      // find the node living in parent element
      pit = factor.find(cit->first);
      
      if (pit == factor.end()) {
	pit = schur.find(cit->first);
	assert(pit != schur.end());
	q.fs_p = 1; 
      } else {
	q.fs_p = 0; 
      }
      
      // set offset and dof
      q.offs_p = pit->second;
      q.n_dof  = cit->first->get_n_dof();
      
      // push into bijection
      bijection.push_back(q);
    }


    // condense
    this->mapper.clear();

    if (!bijection.size()) return;
    
    q = bijection.at(0);
    for (int i=0;i<bijection.size();++i) {
      if (i) {
	int flag = (bijection.at(i).fs_p   - bijection.at(i-1).fs_p);	  
	int diff = (bijection.at(i).offs_p - bijection.at(i-1).offs_p);
	
	if (!flag && (diff == bijection.at(i-1).n_dof)) {
	  this->mapper.back().n_dof += bijection.at(i).n_dof;
	} else {
	  this->mapper.push_back( bijection.at(i) );
	}
      } else {
	this->mapper.push_back( bijection.at(i) );
      }
    }
  }

  inline void Helper_::merge_A()     { _merge_A(); }
  inline void Helper_::branch_ABR()  { _branch_ABR(); }

  inline void Helper_::merge_rhs_x() { _merge_rhs(0,   false); }
  inline void Helper_::merge_rhs_b() { _merge_rhs(10,  false); }
  inline void Helper_::merge_rhs_r() { _merge_rhs(100, true); }

  inline void Helper_::branch_rhs_x() { _branch_rhs(0); }
  inline void Helper_::branch_rhs_b() { _branch_rhs(10); }
  inline void Helper_::branch_rhs_r() { _branch_rhs(100); }

  inline bool Helper_::disp() { 
    this->disp(stdout);
    return true;
  }
  inline bool Helper_::disp(FILE *stream) {
    fprintf(stream, "- Helper -\n");
    fprintf(stream, "  parent [ %d ], child [ %d ] \n",
	    this->p->get_id(), this->c->get_id());
    
    std::vector<Mapper_>::iterator it;
    for (it=this->mapper.begin();it<this->mapper.end();++it) 
      fprintf(stream, 
	      "  side :: [ %d ], offset :: ( %d , %d ), n_dof :: < %d >\n",
	      it->fs_p, it->offs_p, it->offs_c, it->n_dof);
    
    return true;
  }

  // --------------------------------------------------------------
  // ** Protected 
  inline void Helper_::_merge_A() {
    if (!this->mapper.size()) return;

    Matrix parent = this->p->get_matrix();
    Matrix child  = this->c->get_matrix();
    int is_erase  = false;

    int n         = this->mapper.size();
    //int n = 8;

#pragma unroll(UHM_UNROLL_N)
    for (int j = 0 ;j < n; ++j) {

#pragma unroll(UHM_UNROLL_N)
      for (int i = 0 ;i < n; ++i) {

        int mat     = _get_A(this->mapper.at(i).fs_p,
                             this->mapper.at(j).fs_p);
        int ioffs_c = this->mapper.at(i).offs_c; int joffs_c = this->mapper.at(j).offs_c;
        int ioffs_p = this->mapper.at(i).offs_p; int joffs_p = this->mapper.at(j).offs_p;
        int in_dof  = this->mapper.at(i).n_dof;  int jn_dof  = this->mapper.at(j).n_dof;

        parent->merge( child,
                       UHM_ABR,
                       ioffs_c, joffs_c,
                       mat,
                       ioffs_p, joffs_p,
                       in_dof, jn_dof,
                       is_erase);
      }
    }
  }
  inline void Helper_::_branch_ABR() {
    if (!this->mapper.size()) return;

    Matrix parent = this->p->get_matrix();
    Matrix child  = this->c->get_matrix();
    int is_erase  = false;

    int n         = this->mapper.size();

#pragma unroll(UHM_UNROLL_N)
    for (int j = 0 ;j < n; ++j) {

#pragma unroll(UHM_UNROLL_N)
      for (int i = 0 ;i < n; ++i) {

        int mat     = _get_A(this->mapper.at(i).fs_p, this->mapper.at(j).fs_p);
        int ioffs_c = this->mapper.at(i).offs_c; int joffs_c = this->mapper.at(j).offs_c;
        int ioffs_p = this->mapper.at(i).offs_p; int joffs_p = this->mapper.at(j).offs_p;
        int in_dof  = this->mapper.at(i).n_dof;  int jn_dof  = this->mapper.at(j).n_dof;

        child->copy( parent,
                     mat,
                     ioffs_p, joffs_p,
                     UHM_ABR,
                     ioffs_c, joffs_c,
                     in_dof, jn_dof,
                     is_erase);
      }
    }
  }
  inline void Helper_::_merge_rhs(int kind, int is_pivot_applied) {
    if (!this->mapper.size()) return;

    // for checking routine it is necessary to apply pivot
    // at this level -- tricky...
    if (is_pivot_applied &&
        this->p->get_matrix()->is_created( _get_rhs(kind) ) &&
        this->p->get_matrix()->is_buffer( _get_rhs(kind) ) ){
      this->p->get_matrix()->backup   ( _get_rhs(kind) );
      this->p->get_matrix()->set_zero ( _get_rhs(kind) );
    }

    Matrix parent   = this->p->get_matrix();
    Matrix child    = this->c->get_matrix();
    int    is_erase = true;
    int    n_rhs    = p->get_matrix()->get_n_rhs();

    int n           = this->mapper.size();

#pragma unroll(UHM_UNROLL_N)
    for (int j = 0; j < n ; ++j) {
      int mat_src = _get_rhs(1+kind);
      int mat_tgt = _get_rhs(this->mapper.at(j).fs_p+kind);

      int joffs_c = this->mapper.at(j).offs_c;
      int joffs_p = this->mapper.at(j).offs_p;

      int jn_dof  = this->mapper.at(j).n_dof;

      parent->merge( child, 
		     mat_src,
                     joffs_c, 0,
                     mat_tgt,
                     joffs_p, 0,
                     jn_dof, n_rhs, 
		     is_erase);

      // after merge set child bb zero
    }

    if (is_pivot_applied &&
        this->p->get_matrix()->is_created( _get_rhs(kind) ) &&
        this->p->get_matrix()->is_buffer( _get_rhs(kind) ) ){
      int is_merge = true;

      this->p->get_matrix()->apply_pivots( _get_rhs(kind) );
      this->p->get_matrix()->restore( _get_rhs(kind), is_merge );
    }
  }
  inline void Helper_::_branch_rhs(int kind) {
    if (!this->mapper.size()) return;

    Matrix parent   = this->p->get_matrix();
    Matrix child    = this->c->get_matrix();
    int    is_erase = false;
    int    n_rhs    = p->get_matrix()->get_n_rhs();

    int n           = this->mapper.size();

#pragma unroll(UHM_UNROLL_N)
    for (int j = 0; j < n ; ++j) {

      int mat_src = _get_rhs(this->mapper.at(j).fs_p+kind);
      int mat_tgt = _get_rhs(1+kind);

      int joffs_c = this->mapper.at(j).offs_c;
      int joffs_p = this->mapper.at(j).offs_p;

      int jn_dof  = this->mapper.at(j).n_dof;

      child->copy( parent,
                   mat_src,
                   joffs_p, 0,
                   mat_tgt,
                   joffs_c, 0,
                   jn_dof, n_rhs, 
		   is_erase );
    }
  }

  inline int Helper_::_get_rhs(int col) {
    switch (col) {
    case 0:   return UHM_XT;
    case 1:   return UHM_XB;
    case 10:  return UHM_BT;
    case 11:  return UHM_BB;
    case 100: return UHM_RT;
    case 101: return UHM_RB;
    }
    return 0;
  }
  inline int Helper_::_get_A(int col, int row) {
    switch (col*10+row) {
    case 0:  return UHM_ATL;
    case 1:  return UHM_ATR;
    case 10: return UHM_ABL;
    case 11: return UHM_ABR;
    }
    return 0;
  }
}


#endif
