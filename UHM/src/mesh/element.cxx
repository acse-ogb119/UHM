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

#include "uhm/object.hxx"
#include "uhm/mesh/node.hxx"
#include "uhm/mesh/element.hxx"

#include "uhm/matrix/uhm/matrix.hxx"

#include "uhm/util.hxx"


namespace uhm {
  // --------------------------------------------------------------
  // ** Element
  Element_::Element_()                { _init( 0,   0); }
  Element_::Element_(int id)          { _init(id,   0); }
  Element_::Element_(int id, int gen) { _init(id, gen); }
  Element_::~Element_() { 
    if (this->is_matrix_created()) 
      delete this->get_matrix();
  }

  // --------------------------------------------------------------
  void Element_::reset_parent()   { this->parent = nil_element; }
  void Element_::reset_children() { this->children.clear(); }
  void Element_::reset_nodes()    { this->nodes.clear(); }
  void Element_::reset_factor()   { this->factor.clear(); }
  void Element_::reset_schur()    { this->schur.clear(); }

  void Element_::update_generation() {
    if (this->is_leaf()) return;

    int gen = 0;
    std::vector< Element >::iterator it;
    for (it=this->children.begin();it!=this->children.end();++it)
      gen = min(gen, (*it)->get_generation());

    this->generation = gen - 1;
  }

  void Element_::set_matrix(Matrix hm) { this->hm     = hm; }
  void Element_::set_reuse(int flag)   { this->reuse  = flag; }
  void Element_::set_marker(int index, int marker){ 
    assert(index > -1 && index < 2);
    this->marker[index] = marker;
  }
  void Element_::set_parent(Element p) {
    assert(element_valid(p));
    this->parent = p;
  }

  int  Element_::get_generation()  { return this->generation; }
  int  Element_::get_height() {
    Element e = this;

    int cnt = 0;
    while (!e->is_orphan()) { e = e->get_parent(); ++cnt; }
    return cnt;
  }
  
  Element Element_::get_parent()  { return this->parent; }
  Element Element_::get_child(int c) { 
    assert(this->children.size() > c);
    return this->children.at(c);
  }


  Matrix Element_::get_matrix()          { return this->hm; }
  int    Element_::get_marker(int index) { 
    assert(index>-1 && index<2);
    return this->marker[index];
  }

  int  Element_::get_n_children()      { return this->children.size(); }
  int  Element_::get_n_nodes()         { return this->nodes.size(); }
  int  Element_::get_n_factor_nodes()  { return this->factor.size(); }
  int  Element_::get_n_schur_nodes()   { return this->schur.size(); }
  
  std::pair<int,int> Element_::get_n_dof() {
    std::pair<int,int> ret(0,0);
    std::map< Node, int >::iterator it;
    for (it=this->nodes.begin();it!=this->nodes.end();++it) {
      switch (it->second) {
      case UHM_SEPARATED_FACTOR : 
	ret.first  += it->first->get_n_dof();
	break;
      case UHM_SEPARATED_SCHUR  :
	ret.second += it->first->get_n_dof();
	break;
      }
    }
    return ret;
  }

  void Element_::estimate_cost( int method, int datatype, int n_rhs,
                                double &flop_decompose, 
                                double &flop_solve,
                                unsigned int &n_nonzero_factor,
                                double &buffer ) {

    double is_complex = ( datatype == UHM_COMPLEX );
    std::pair<int,int> n_dof = this->get_n_dof();
    int fs = n_dof.first, ss = n_dof.second;

    switch (method) {
    case UHM_CHOL:
      flop_decompose = ( linal::get_flop_chol      (is_complex, fs    ) +         //chol:ATL
                         linal::get_flop_trsm_lower(is_complex, fs, ss) +         //trsm:ABL<-ATL
                         linal::get_flop_syrk      (is_complex, fs, ss) );        //syrk:ABR<-ABL

      flop_solve = ( linal::get_flop_trsm_lower(is_complex,    fs, n_rhs    ) +   //trsm:BT<-ATL
                     linal::get_flop_gemm      (is_complex,    ss, n_rhs, fs) +   //gemm
                     linal::get_flop_trsm_upper(is_complex, n_rhs,    fs    ) +   //trsm:BT<-ATL
                     linal::get_flop_gemm      (is_complex,    fs, n_rhs, ss) );  //gemm
      
      buffer = ( linal::get_memory(is_complex, fs, fs) +                          //ATL
                 linal::get_memory(is_complex, ss, fs) );                         //ABL

      n_nonzero_factor = fs*fs + ss*fs;

      break;
    case UHM_LU_NOPIV:
    case UHM_LU_PIV:
    case UHM_LU_INCPIV:
      flop_decompose = ( linal::get_flop_lu        (is_complex, fs, fs) +         //lu:  ATL
                         linal::get_flop_trsm_lower(is_complex, fs, ss) +         //trsm:ABL<-ATL
                         linal::get_flop_trsm_upper(is_complex, ss, fs) +         //trsm:ATR<-ATL
                         linal::get_flop_gemm      (is_complex, ss, ss, fs) );    //gemm

      flop_solve = ( linal::get_flop_trsm_lower(is_complex,    fs, n_rhs    ) +   //trsm:BT<-ATL
                     linal::get_flop_gemm      (is_complex,    ss, n_rhs, fs) +   //gemm
                     linal::get_flop_trsm_upper(is_complex, n_rhs,    fs    ) +   //trsm:BT<-ATL
                     linal::get_flop_gemm      (is_complex,    fs, n_rhs, ss) );  //gemm

      buffer = ( linal::get_memory(is_complex, fs, fs) +                          //ATL
                 linal::get_memory(is_complex, ss, fs) +                          //ABL
                 linal::get_memory(is_complex, fs, ss) );                         //ATR

      n_nonzero_factor = fs*fs + 2.0*ss*fs;
      break;
    case UHM_QR:
      // I am not sure about this...
      flop_decompose = ( linal::get_flop_qr        (is_complex, fs, fs) +         //lu:  ATL
                         linal::get_flop_trsm_lower(is_complex, fs, ss) +         //trsm:ABL<-ATL
                         linal::get_flop_q         (is_complex, ss, fs, fs) +     //q:   ATR<-ATL
                         linal::get_flop_gemm      (is_complex, ss, ss, fs) );    //gemm

      flop_solve = ( linal::get_flop_q         (is_complex,    fs, n_rhs, fs) +   //q:BT<-ATL
                     linal::get_flop_gemm      (is_complex,    ss, n_rhs, fs) +   //gemm
                     linal::get_flop_trsm_upper(is_complex, n_rhs,    fs    ) +   //trsm:BT<-ATL
                     linal::get_flop_gemm      (is_complex,    fs, n_rhs, ss) );  //gemm

      buffer = ( linal::get_memory(is_complex, fs, fs) +                          //ATL
                 linal::get_memory(is_complex, fs, fs) +                          //Q
                 linal::get_memory(is_complex, ss, fs) +                          //ABL
                 linal::get_memory(is_complex, fs, ss) );                         //ATR

      n_nonzero_factor = 2.0*fs*fs + 2.0*ss*fs;

      break;
    }
  }
  
  bool Element_::is_orphan()         { return (this->parent == nil_element);}
  bool Element_::is_leaf()           { return this->children.empty(); }
  bool Element_::is_nodes_separated(){
    int n_factor=0, n_schur=0;
    std::map< Node,int >::iterator it;
    for (it=this->nodes.begin();it!=this->nodes.end();++it) {
      switch (it->second) {
      case UHM_SEPARATED_FACTOR : ++n_factor; break; 
      case UHM_SEPARATED_SCHUR  : ++n_schur;  break;
      }
    }
    return (this->nodes.size() == (n_factor+n_schur));
  }
  bool Element_::is_nodes_arranged() { 
    return ( this->nodes.size() && 
	     (this->factor.size() || this->schur.size()) );
  }

  bool Element_::is_matrix_created()  { return ( this->hm != nil_matrix ); }
  bool Element_::is_matrix_reusable() { return this->reuse; }

  void Element_::collect_leaf_children( int n_max, int &n_leaves, 
					Element *leaves ) {
    std::vector< Element > tmp;
    this->collect_leaf_children( tmp );

    n_leaves = tmp.size();
    assert(n_max >= n_leaves);

    for (int i=0;i<n_leaves;++i)
      leaves[i] = tmp.at(i);
  }

  void Element_::collect_leaf_children( std::vector< Element > &leaves ) {
    leaves.clear();
    if (this->is_leaf()) return;
    
    std::vector< Element > tmp;
    tmp.push_back( this );
    
    for (int i=0;i<tmp.size();++i) 
      for (int j=0;j<tmp.at(i)->get_n_children();++j) 
	tmp.push_back( tmp.at(i)->get_child(j) );
    
    for (int i=0;i<tmp.size();++i)
      if (tmp.at(i)->is_leaf()) leaves.push_back(tmp.at(i));
  }
  
  void Element_::add_child(Element c) { 
    assert(element_valid(c));
#pragma omp critical
    this->children.push_back(c);
  }

  void Element_::add_node(Node n) { this->add_node(n,0);  }
  void Element_::add_node(Node n, int separated) {
    assert(node_valid(n));

    std::pair< std::map< Node, int >::iterator, bool> ret;

    ret = this->nodes.insert( std::pair<Node,int>(n,separated) );
    ret.first->second = separated;

    // if newly added add the owner
    if (ret.second) n->add_owner(this);
  }

  // manual set up for factor and schur nodes
  void Element_::add_factor(Node n, int offs) {
    this->factor.push_back(std::pair<Node,int>(n, offs));
  }
  void Element_::add_schur(Node n, int offs) {
    this->schur.push_back(std::pair<Node,int>(n, offs));
  }

  void Element_::separate_nodes() {
    std::map< Node,int >::iterator it;
    for (it=this->nodes.begin();it!=this->nodes.end();++it) {

      // if node is owned by no element, error!!!!
      assert(it->first->get_n_owner());
      
      // if the node is owned by one element, it can be eliminated
      if (it->first->get_n_owner() == 1 &&
          it->first->get_kind() == UHM_NODE_KIND_DEFAULT) 
	it->second = UHM_SEPARATED_FACTOR;
      else 
	it->second = UHM_SEPARATED_SCHUR;
    }
  }

  void Element_::arrange_nodes() {
    // assumption :: nodes are separated and numbered
    // clear storage
    this->reset_factor();
    this->reset_schur();

    std::map < int , Node > factor_nodes, schur_nodes;
    std::pair< int , Node > in;

    // ** sort with global offsets
    std::map< Node, int  >::iterator nit;
    for (nit=this->nodes.begin();nit!=this->nodes.end();++nit) {

      // node has n_dof
      if (nit->first->get_n_dof()) {
	assert(nit->second==UHM_SEPARATED_FACTOR || 
	       nit->second==UHM_SEPARATED_SCHUR);
	
	// key is offset value so, tmp contains the nodes ascending order
	in.first  = nit->first->get_offset();
	in.second = nit->first;
	
	switch (nit->second) {
	case UHM_SEPARATED_FACTOR : factor_nodes.insert( in ); break;
	case UHM_SEPARATED_SCHUR  : schur_nodes.insert( in );  break;
	}
      }
    }

    int offs;
    std::map< int, Node >::iterator tit;

    // ** arrange the nodes according to the sorted order
    // add to factor
    offs = 0;
    for (tit=factor_nodes.begin();tit!=factor_nodes.end();++tit) {
      this->add_factor( tit->second, offs );
      offs += (tit->second->get_n_dof());
    }

    // add to schur 
    offs = 0;
    for (tit=schur_nodes.begin();tit!=schur_nodes.end();++tit) {
      this->add_schur( tit->second, offs );
      offs += (tit->second->get_n_dof());
    }
  }

  void Element_::merge_nodes_from_children() {
    // if this is leaf, nothing to do. Else, clear 'nodes' container
    if (this->is_leaf()) return;

    this->reset_nodes();
    
    // collect all children
    std::vector< Element >::iterator it;
    for (it=this->children.begin();it!=this->children.end();it++) 
      this->merge_nodes((*it));
  }

  void Element_::merge_nodes(Element c) {
    assert(element_valid(c));

    std::map< Node, int >::iterator it;
    for (it=c->nodes.begin();it!=c->nodes.end();++it) {
      if (it->second == UHM_SEPARATED_SCHUR) {

	Node n = it->first;
	this->add_node( n );
	n->remove_owner( c );
      }
    }
  }

  void Element_::restore_connectivity() {
    if (!this->is_leaf()) return;
    
    // for leaf elements, restore connectivity to initial mesh
    std::map< Node,int >::iterator it;
    for (it=this->nodes.begin();it!=this->nodes.end();++it) {

      // separated nodes are recovered
      it->second = UHM_NOT_SEPARATED;

      // restore the connectivity as initial mesh status
      Node n = it->first;
      n->clean_connectivity();
      n->add_owner(this);

      // reset factor and schur
      this->reset_factor();
      this->reset_schur();
    }
  }

  void Element_::numbering() {
    std::map< Node, int >::iterator it;
    for (it=this->nodes.begin();it!=this->nodes.end();++it) {
      if (it->second == UHM_SEPARATED_FACTOR) {
	Node n = it->first;
	n->set_offset(get_g_offset());
	add_g_offset(n->get_n_dof());
      }
    }
  }

  bool Element_::disp() { return this->disp(stdout); }
  //bool Element_::disp(FILE *stream) { return this->disp(stream, 0); }
  bool Element_::disp(FILE *stream) {
    // if mode exist, mode 0 - brief, mode 1 - detail
    
    fprintf(stream, "-Element-\n");

    int n_rank = 0;
//     if (this->comm != MPI_COMM_NULL) 
//       MPI_Comm_size(this->comm, &n_rank);
//     else 
//       n_rank = -1;
      
    // general data
    fprintf(stream, "  id [ %d ] , generation [ %d ], is_reusable [ %d ], Commsize [ %d ]\n",
	    this->get_id(), this->get_generation(), 
	    this->is_matrix_reusable(), n_rank);

    // tree info
    if (!this->is_orphan()) {
      assert(element_valid(this->parent));
      fprintf(stream, "  parent [ %d ]\n", this->parent->get_id());
    }

    if (this->get_n_children()) {
      std::vector< Element >::iterator it;
      fprintf(stream, "  %d children [ ", this->get_n_children());
      for (it=this->children.begin();it<this->children.end();++it) {
	if (element_valid(*it)) 
	  fprintf(stream, " %d ", (*it)->get_id());
      }
      fprintf(stream, " ]\n");
    }
    
    // nodal conncectivity
    if (this->get_n_nodes()) {
      printf("  n_nodes : %d, n_factor : %d, n_schur : %d\n",
	     this->get_n_nodes(), 
	     this->get_n_factor_nodes(), this->get_n_schur_nodes());
      if (this->is_nodes_arranged()) {
	std::vector< std::pair<Node, int> >::iterator it;
	fprintf(stream, "  %d factor (id):<offset> [ ", 
		this->get_n_factor_nodes());
	for (it=this->factor.begin();it!=this->factor.end();++it) 
	  fprintf(stream, " (%d, %d):<%d> ", 
		 (it->first)->get_id().first, 
		 (it->first)->get_id().second,
		 (it->second));
	fprintf(stream, " ]\n");
	fprintf(stream, "  %d schur (id):<offset> [ ", 
		this->get_n_schur_nodes());
	for (it=this->schur.begin();it!=this->schur.end();it++) 
	  fprintf(stream, " (%d, %d):<%d> ", 
		 (it->first)->get_id().first, 
		 (it->first)->get_id().second,
		 (it->second));
	fprintf(stream, " ]\n");
      } else {
	std::map< Node,int >::iterator it;
	fprintf(stream, "  %d nodes (id):<separation>[ ", this->get_n_nodes());
	for (it=this->nodes.begin();it!=this->nodes.end();it++) 
	  if (node_valid(it->first)) 
	    fprintf(stream, " (%d, %d):<%d> ", 
		    it->first->get_id().first, 
		    it->first->get_id().second,
		    it->second);
	fprintf(stream, " ]\n");
      }
    }
    return true;
  }

  bool Element_::write_graphviz_hier(FILE *fp, 
				     int is_leaf2root, double max_n_dof) {
    std::pair<int,int> n_dof = this->get_n_dof();
    double size = (n_dof.first + n_dof.second)/max_n_dof;

    fprintf(fp, "%d [label=\"Elt_%d\",shape=box,width=%lf,height=%lf,style=filled,color=peru]; %d -> {",
	    this->get_id(), this->get_id(), 
	    size,size, 
	    this->get_id());
    if (is_leaf2root) {
      if (!this->is_orphan())
	fprintf(fp,"%d;", this->get_parent()->get_id());
    } else {
      for (int i=0;i<this->get_n_children();++i)
	fprintf(fp,"%d;", this->get_child(i)->get_id());
    }
    fprintf(fp,"};\n");
    return true;
  }

  bool Element_::export_matrix(int n_rhs, 
			       int &m, int &n,
			       std::vector< std::pair<int,int> > &ij, 
			       int mat) {
    // export matrix unassembled form for the leaf
    // assert(this->is_leaf());

    std::vector< std::pair<Node, int> >::iterator k1, k2, 
      k1_begin, k2_begin, k1_end, k2_end;

    // export connectivity only use e->get_matrix()->export_matrix()
    // for values
    //linal::Flat_ matrix;
    std::pair< int, int > n_dof = this->get_n_dof();

    switch (mat) {
    case UHM_ATL:
      k2_begin = this->factor.begin();
      k2_end   = this->factor.end();
      k1_begin = this->factor.begin();
      k1_end   = this->factor.end();

      m = n_dof.first; 
      n = n_dof.first;

      if (!m) return true;
      break;
    case UHM_ATR:
      k2_begin = this->schur.begin();
      k2_end   = this->schur.end();
      k1_begin = this->factor.begin();
      k1_end   = this->factor.end();

      m = n_dof.first; 
      n = n_dof.second;

      if (!m || !n) return true;
      break;
    case UHM_ABL:
      k2_begin = this->factor.begin();
      k2_end   = this->factor.end();
      k1_begin = this->schur.begin();
      k1_end   = this->schur.end();

      m = n_dof.second; 
      n = n_dof.first;

      if (!m || !n) return true;
      break;
    case UHM_ABR:
      k2_begin = this->schur.begin();
      k2_end   = this->schur.end();
      k1_begin = this->schur.begin();
      k1_end   = this->schur.end();

      m = n_dof.second; 
      n = n_dof.second;

      if (!m) return true;
      break;
    case UHM_BT:
    case UHM_XT:
    case UHM_RT:

      k2_begin = this->factor.begin();
      k2_end   = this->factor.end();
      k1_begin = this->factor.begin();
      k1_end   = this->factor.end();

      m = n_dof.first;
      n = n_rhs;

      if (!m) return true;
      break;
    case UHM_BB:
    case UHM_XB:
    case UHM_RB:

      k2_begin = this->schur.begin();
      k2_end   = this->schur.end();
      k1_begin = this->schur.begin();
      k1_end   = this->schur.end();

      m = n_dof.second;
      n = n_rhs;

      if (!m) return true;
      break;
    default:
      fprintf(stderr,"export_matrix: not support\n");
      return false;
    }

    if (mat < UHM_P) {
      for (k2=k2_begin;k2!=k2_end;++k2) {
	int k2_dof      = (k2->first)->get_n_dof();
	int k2_offs_tgt = (k2->first)->get_offset();
	
	for (int l2=0;l2<k2_dof;++l2) {

	  for (k1=k1_begin;k1!=k1_end;++k1) {
	    int k1_dof      = (k1->first)->get_n_dof();
	    int k1_offs_tgt = (k1->first)->get_offset();
	  
	    for (int l1=0;l1<k1_dof;++l1) {
	      std::pair<int,int> in;
	      in.first  = k1_offs_tgt+l1;
	      in.second = k2_offs_tgt+l2;
	      
	      ij.push_back( in );
	    }
	  }
	}
      }  
    } else if (mat > UHM_T && mat < UHM_END) {
      for (int j=0;j<n_rhs;++j) {
	for (k1=k1_begin;k1!=k1_end;++k1) {
	  
	  int k1_dof      = (k1->first)->get_n_dof();
	  int k1_offs_tgt = (k1->first)->get_offset();
	  
	  for (int l1=0;l1<k1_dof;++l1) {
	    std::pair<int,int> in;
	    in.first  = k1_offs_tgt+l1;
	    in.second = j;

	    ij.push_back( in );
	  }
	}
      }
    }
    return true;
  }

  bool Element_::sparse_pattern(std::set< std::pair<int,int> > &s, int mat) {

    // get the sparse pattern from unassembled matrix
    int m=0, n=0;
    std::vector< std::pair<int,int> > ij;
    this->export_matrix( 0, m, n, ij, mat );

    std::vector< std::pair<int,int> >::iterator it;
    for (it=ij.begin();it!=ij.end();++it)
      s.insert(*it);
    
    return true;
  }

  bool Element_::export_matrix(FILE *stream, int n_rhs, int mat) {

    fprintf(stream, "###  element id : %d\n", this->get_id());
    fprintf(stream, "%d\n", mat);

    int m=0, n=0;
    std::vector< std::pair<int,int> > ij;
    this->export_matrix( n_rhs, m, n, ij, mat );

    fprintf(stream, "%d %d\n", m, n);

    std::vector< std::pair<int,int> >::iterator it;
    for (it=ij.begin();it!=ij.end();++it)
      fprintf(stream, "%10d %10d\n", it->first, it->second);
    
    return true;
  }

  
}
