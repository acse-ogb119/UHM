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

#include "uhm/object.hxx"

#include "uhm/operation/scheduler.hxx"
#include "uhm/operation/element.hxx"

#include "uhm/mesh/node.hxx"
#include "uhm/mesh/element.hxx"

#include "uhm/matrix/uhm/matrix.hxx"

#include "uhm/mesh/mesh.hxx"

#include "uhm/interf/sparse.hxx"


namespace uhm {
  // --------------------------------------------------------------
  // ** Callable from C
  //void mesh_new   ( Mesh &m ) { m = new Mesh_; }
  //void mesh_delete( Mesh &m ) { delete m; }
  // --------------------------------------------------------------
  // ** Mesh
  Mesh_::Mesh_()                       { _init( 0, 0); }
  Mesh_::Mesh_(int id)                 { _init(id, 0); }
  Mesh_::Mesh_(int id, int id_element) { _init(id, id_element); }
  Mesh_::~Mesh_()                      {  
    this->reset();
  }

  // --------------------------------------------------------------
  bool Mesh_::is_elements_separated() {
    int cnt = 0;
    std::map< int, Element_ >::iterator eit;
    for (eit=this->elements.begin();eit!=this->elements.end();++eit) {
      Element e = &(eit->second);
      if (!e->is_nodes_separated()) {
	e->disp();
	++cnt;
      }
    }
    if (cnt) {
      printf("elements in mesh : %d, elements not separated : %d\n",
	     (int)this->elements.size(), cnt);
      return false;
    }
    return true;
  }

  bool Mesh_::is_nodes_numbered() {
    int cnt = 0;
    std::set< int > numbers;
    std::map< int, Element_ >::iterator eit;
    for (eit=this->elements.begin();eit!=this->elements.end();++eit) {
      Element e = &(eit->second);

      std::pair< std::set<int>::iterator,bool > ret;
      std::map< Node,int >::iterator nit;
      for (nit=e->nodes.begin();nit!=e->nodes.end();++nit) {
	if (nit->second == UHM_SEPARATED_FACTOR && 
	    nit->first->get_n_dof()) {
	  ret = numbers.insert( nit->first->get_offset() );
	  if (ret.second == false) {
	    nit->first->disp();
	    ++cnt;
	  }
	}
      }
    }
    if (cnt) {
      printf("nodes in mesh: %d, nodes not numbered : %d\n",
	     (int)this->nodes.size(), cnt);
      return false;
    }
    return true;
  }

  void Mesh_::backup(Mesh backup) {
    assert(mesh_valid(backup));
    backup->reset();

    // ** backup only store the connectivity without matrices
    //    the backup mesh is used to compare the previous connectivity
    //    for 'updating factorization'

    std::map< int, Element_ >::iterator eit;
    for (eit=this->elements.begin();eit!=this->elements.end();eit++) {
      Element a = &(eit->second);
      Element b = backup->insert_element(a->get_id());

      /*
        { // add nodes as separated
	std::map< Node, int >::iterator nit;
	for (nit=a->nodes.begin();nit!=a->nodes.end();++nit) {
        Node n = nit->first;
        b->add_node( n, nit->second );
	}
        }
      */
      { // add factor and schur nodes with encoded value
	std::vector< std::pair<Node, int> >::iterator nit;
	
	// store the pointer and node info 
	for (nit=a->factor.begin();nit!=a->factor.end();++nit) { 
	  Node n = nit->first;
	  b->add_factor( n, ( n->get_n_dof() + 
			      n->get_p() + 
			      n->get_kind() ) );
	}
	
	for (nit=a->schur.begin();nit!=a->schur.end();++nit) {
	  Node n = nit->first;
	  b->add_schur( n, ( n->get_n_dof() + 
			     n->get_p() + 
			     n->get_kind() ) );
	}
      }
    }
  }

  void Mesh_::check_reuse(Mesh backup) {
    assert(mesh_valid(backup));

    std::map< int, Element_ >::iterator eit;

    // initially set reuse flag true
    for (eit=this->elements.begin();eit!=this->elements.end();++eit) {
      Element a = &(eit->second);
      a->set_reuse(true);
    }
    
    // check all elements
    for (eit=this->elements.begin();eit!=this->elements.end();++eit) {

      // a - current , b - backup
      Element a = &(eit->second);
      
      bool flag = true;
      
      // check existence of backup
      Element b = backup->find_element(a->get_id());
      
      if (b == nil_element) {
	flag = false;
      } else {
	std::vector< std::pair<Node, int> >::iterator ait;
	std::vector< std::pair<Node, int> >::iterator bit;
	
	if ( a->get_n_factor_nodes() == b->get_n_factor_nodes() &&
	     a->get_n_schur_nodes()  == b->get_n_schur_nodes() ) {
	  
	  // compare factor nodes with b
	  if (flag) {
	    for (int i=0;i<a->get_n_factor_nodes();++i) {
	      Node n  = a->factor.at(i).first;
	      int val = n->get_n_dof() + n->get_p() + n->get_kind();
	      
	      if ( n   != b->factor.at(i).first  ||
		   val != b->factor.at(i).second ) {
		flag = false; break;
	      }
	    }
	  }
	  // compare schur nodes with b
	  if (flag) {
	    for (int i=0;i<a->get_n_schur_nodes();++i) {
	      Node n  = a->schur.at(i).first;
	      int val = n->get_n_dof() + n->get_p() + n->get_kind();
	      
	      if ( n   != b->schur.at(i).first  ||
		   val != b->schur.at(i).second ) {
		flag = false; break;
	      }
	    }
	  }
	} else {
	  flag = false;
	}
      }
      a->set_reuse(flag);
      
      // cannot reuse upper hierarchy, too
      if (!a->is_matrix_reusable()) {
	while (a != nil_element) {
	  a->set_reuse(false);
	  a = a->get_parent();
	}
      }
    }
  }

  void Mesh_::reset() {
    this->unlock();

    this->remove_all_elements();
    this->remove_all_nodes();

    this->id_element = 0;

    //     if (this->comm != MPI_COMM_NULL)
    //       MPI_Comm_free(&this->comm);
    //     this->comm = MPI_COMM_NULL;
  }
  
  Scheduler Mesh_::get_scheduler() { return &this->scheduler; }

  int  Mesh_::get_n_nodes()    { return this->nodes.size(); }
  int  Mesh_::get_n_elements() { return this->elements.size(); }

  void Mesh_::remove_all_nodes()    { this->nodes.clear(); }
  void Mesh_::remove_orphan_nodes() { 
    std::map< std::pair<int,int>, Node_ >::iterator it;
    for (it=this->nodes.begin();it!=this->nodes.end();it++) {
      if (!(it->second.get_n_owner())) this->nodes.erase(it);
    }
  }
  void Mesh_::remove_all_elements() { this->elements.clear(); }
  void Mesh_::remove_lower_elements(int gen) {
    std::map< int, Element_ >::iterator it;

    // first iteration : cut the tree
    for (it=this->elements.begin();it!=this->elements.end();it++) {
      Element parent = it->second.get_parent();
      if (parent && parent->get_generation() < gen) 
	it->second.set_parent(nil_element);
    }
    
    // second iteration : delete neg-generation
    for (it=this->elements.begin();it!=this->elements.end();it++) 
      if (it->second.get_generation() < gen) this->elements.erase(it);
  }

  Node Mesh_::add_node(int id, int n_dof) {
    return (this->add_node(std::pair<int,int>(id, 0), n_dof, n_dof, 0));
  }
  Node Mesh_::add_node(int id, int n_dof, int p) {
    return (this->add_node(std::pair<int,int>(id, 0), n_dof, p, 0));
  }
  Node Mesh_::add_node(std::pair<int,int> id, int n_dof) {
    return (this->add_node(id, n_dof, n_dof, 0));
  }
  Node Mesh_::add_node(std::pair<int,int> id, int n_dof, int p) {
    return (this->add_node(id, n_dof, p, 0));
  }
  Node Mesh_::add_node(std::pair<int,int> id, int n_dof, int p, int kind) {
    std::pair< std::map< std::pair<int,int>, Node_ >::iterator, bool> ret;

#pragma omp critical 
    {
      ret = this->nodes.insert(std::make_pair(id, Node_(id, n_dof, p, kind)));
      
      // if node is already exist, check it is same node
      if (!ret.second) 
	assert((*ret.first).second.get_n_dof() == n_dof &&
	       (*ret.first).second.get_p()     == p &&
	       (*ret.first).second.get_kind()  == kind);
    }

    return &((*ret.first).second);
  }

  Node Mesh_::find_node(int id) {
    return (this->find_node(std::pair<int,int>(id, 0)));
  }
  Node Mesh_::find_node(std::pair<int,int> id) {
    std::map< std::pair<int,int>, Node_ >::iterator it;
    it = this->nodes.find(id);

    // if node is already exist, check it is same node
    if (it != this->nodes.end()) return &(it->second);
    else return nil_node;
  }
  
  bool Mesh_::remove_node(int id) {
    return (this->remove_node(std::pair<int,int>(id, 0)));
  }
  bool Mesh_::remove_node(std::pair<int,int> id) {
    std::map< std::pair<int,int>, Node_ >::iterator it;
    it = this->nodes.find(id);
    
    // if node exist, erase it
    if (it != this->nodes.end()) {
      this->nodes.erase(it);
      return true;
    } else return false;
  }

  Element Mesh_::get_root() {
    std::map< int, Element_ >::iterator it;

    // if empty
    if (!this->elements.size()) return nil_element;

    Element e = &(this->elements.begin()->second);
    while ( !e->is_orphan() ) { e = e->get_parent(); }
    return e;
  }

  Element Mesh_::add_element() {
    return (this->add_element(0));
  }
  Element Mesh_::add_element(int gen) {
    std::pair<std::map< int, Element_ >::iterator, bool> ret;

#pragma omp critical
    {
      ret = this->elements.insert(std::pair<int,Element_>(this->id_element,
							  Element_(this->id_element, gen)));
      assert(ret.second);
      this->id_element++;
    }

    return &((ret.first)->second);
  }

  Element Mesh_::insert_element(int id) {
    std::pair<std::map< int, Element_ >::iterator, bool> ret;
    ret = this->elements.insert(std::pair<int,Element_>(id,
                                                        Element_(id, 0)));
    assert(ret.second);
    return &((ret.first)->second);
  }

  void Mesh_::adjust_element_numbering() {
    int id = 0;
    std::map< int, Element_ >::iterator eit;
    for (eit=this->elements.begin();eit!=this->elements.end();eit++) {
      Element e = &(eit->second);
      id = max(id, e->get_id());
    }
    this->id_element = id + 1;
  }

  Element Mesh_::find_element(int id) {
    std::map< int, Element_ >::iterator it;
    it = this->elements.find(id);
    if (it != this->elements.end()) return &(it->second);
    else return nil_element;
  }

  bool Mesh_::remove_element(int id) {
    std::map< int, Element_ >::iterator it;
    it = this->elements.find(id);
    if (it != this->elements.end()) {
      this->elements.erase(it);
      return true;
    } else return false;
  }

  Element Mesh_::refine_element(int id, int n_children) {
    return this->refine_element(id, true, n_children);
  }
  Element Mesh_::refine_element(int id, int is_binary, int n_children) {
    Element e = this->find_element(id);
    assert(n_children > 1 && e != nil_element);
    
    if (is_binary) {
      std::vector< Element > tree;
      tree.push_back( e );

      for (int i=0;i<n_children;++i) {
	int gen = (tree.at(i)->get_generation() + 1);

	for (int j=0;j<2;++j) {
	  Element c;
	  c = this->add_element(gen);
	  c->set_parent(tree.at(i));
	  tree.at(i)->add_child(c);
	  tree.push_back( c );
	}
      }

    } else {

      // add new element and set up family relation
      for (int i=0;i<n_children;i++) {
	Element c;
	int gen = (e->get_generation() + 1);
	c = this->add_element(gen);
	c->set_parent(e);
	e->add_child(c);
      }
    }

    return e;
  }

  Element Mesh_::unrefine_element(int id) {
    Element e = this->find_element(id);
    assert(e != nil_element);

    for (int i=0;i<e->get_n_children();i++) {
      Element c;
      c = e->get_child(i);
      if (c->is_leaf()) this->unrefine_element(c->get_id());
      this->remove_element(c->get_id());
    }
    e->reset_children();

    return e;
  }

  void Mesh_::lock() {
    this->unlock();

    // scheduler ready
    Scheduler s = this->get_scheduler();
    s->load(this);

    s->execute_leaves_seq(&op_restore_connectivity);
    s->execute_elements_seq(&op_update_connectivity, true);

    reset_g_offset();
    s->execute_elements_seq(&op_numbering, true);

    // arrange is not working
    s->execute_elements_seq(&op_arrange_nodes, true);
    
    // locked
    this->locker = true;
  }

  void Mesh_::unlock() { 
    this->get_scheduler()->unload();
    this->locker = false;
  }
  
  int  Mesh_::is_locked() { return this->locker; }

  bool Mesh_::disp() { return this->disp(stdout); }
  bool Mesh_::disp(int mode) { return this->disp(stdout, mode); }
  bool Mesh_::disp(FILE *stream) { return this->disp(stream, UHM_DISP_ALL); }
  bool Mesh_::disp(FILE *stream, int mode) {
    int n, e, m;
    switch (mode) {
    case UHM_DISP_ALL:        n=1;e=1;m=1;break;
    case UHM_DISP_NODE:       n=1;e=0;m=0;break;
    case UHM_DISP_ELEMENT:    n=0;e=1;m=0;break;
    case UHM_DISP_MATRIX:     n=0;e=1;m=1;break;
    }

    fprintf(stream, "-Mesh-\n");

    // general data
    fprintf(stream, "  id [ %d ] , n_nodes [ %d ], n_elements [ %d ]\n",
	    this->get_id(), this->get_n_nodes(), 
	    this->get_n_elements());
    
    // node disp
    if (n && this->get_n_nodes()) {
      fprintf(stream, "- Node -\n");
      fprintf(stream, "-------------------------------------------\n");
      std::map< std::pair<int,int>, Node_ >::iterator it;
      for (it=this->nodes.begin();it!=this->nodes.end();it++) 
	(*it).second.disp(stream);
      fprintf(stream, "-------------------------------------------\n");
    }
    
    // element disp
    if (e && this->get_n_elements()) {
      fprintf(stream, "- Element -\n");
      fprintf(stream, "-------------------------------------------\n");
      std::map< int, Element_ >::iterator it;
      for (it=this->elements.begin();it!=this->elements.end();it++) {
	(*it).second.disp(stream);
        // matrix disp
	if (m && (*it).second.is_matrix_created()) 
	  (*it).second.get_matrix()->disp(stream);
      }
      fprintf(stream, "-------------------------------------------\n");
    }

    return true;
  }

  bool Mesh_::import_file(char *full_path) {
    // empty all containers in mesh
    this->reset();
    
    FILE *fp;
    char format[32], *line;
    int i, j, n_nodes, n_elements;

    // *** file open ASCII mode 
    assert(open_file(full_path, "r", &fp));
    assert(read_line(fp, &line));
    sscanf(line, "%s", format);
    
    // *** format check
    if (strcmp("UHM", format)) {
      fprintf(stderr, "mismatch format %s\n", format);
      abort();
    } else {
      printf("format %s\n", format);
    }

    // *** read nodes
    assert(read_line(fp, &line));
    sscanf(line, "%d", &n_nodes);
    printf("Reading : n_nodes %d\n", n_nodes);

    for (i=0;i<n_nodes;i++) {
      int id, phy, ndof, p, kind;
      assert(read_line(fp, &line));
      sscanf(line, "%d %d %d %d %d", 
	     &id, &phy, &ndof, &p, &kind);
      this->add_node(std::pair<int,int>(id, phy), ndof, p, kind);
    }

    // *** read elements
    assert(read_line(fp, &line));
    sscanf(line, "%d", &n_elements);
    printf("Reading : n_elements %d\n", n_elements);

    for (i=0;i<n_elements;i++) {
      int n_nods;

      assert(read_line(fp, &line));
      sscanf(line, "%d", &n_nods);

      Element e = this->add_element();
      for (j=0;j<n_nods;j++) {
	int id, phy;

	assert(read_line(fp, &line));
	sscanf(line, "%d %d", &id, &phy);
	Node n = this->find_node(std::pair<int,int>(id, phy));
	if (n != nil_node) {
	  e->add_node(n);
	} else { 
	  fprintf(stderr, "failt to find node [ %d , %d ]\n", id, phy);
	  abort();
	}
      }
    }

    // *** clean up orphans
    this->remove_orphan_nodes();

    // *** close file
    assert(close_file(fp));

    return true;
  }

  bool Mesh_::export_graphviz_hier(char *full_path, int is_leaf2root) {
    FILE *fp;

    // *** file open ASCII mode
    assert(open_file(full_path, "w", &fp));
    linal::head_graphviz(fp, "Mesh_hierarchy");

    int max_n_dof=0;
    std::map< int, Element_ >::iterator it;
    for (it=this->elements.begin();it!=this->elements.end();++it) {
      Element e = &(it->second);
      std::pair<int,int> dof = e->get_n_dof();
      max_n_dof = max(max_n_dof, dof.first+dof.second);
    }

    for (it=this->elements.begin();it!=this->elements.end();++it) {
      Element e = &(it->second);
      e->write_graphviz_hier(fp, is_leaf2root, max_n_dof);
    }
          
    linal::tail_graphviz(fp);

    // *** close file
    assert(close_file(fp));
    return true;
  }

  bool Mesh_::export_sparse_pattern(char *full_path, char *ss, 
				    int is_fill_in) {
    FILE *fp;

    // *** file open ASCII mode
    assert(open_file(full_path, "w", &fp));

    std::map< int, Element_ >::iterator it;
    int n_dofs=0;
    for (it=this->elements.begin();it!=this->elements.end();++it) {
      Element e = &(it->second);
      n_dofs += e->get_n_dof().first;
    }

    fprintf(fp, "%% UHM export sparse pattern\n");
    fprintf(fp, "clear all; close all; \n");

    std::set< std::pair<int,int> > sparse;
    if (is_fill_in) {
      for (it=this->elements.begin();it!=this->elements.end();++it) {
	Element e = &(it->second);
	for (int i=UHM_ATL;i<UHM_ABR;++i)
	  assert(e->sparse_pattern(sparse, i));
      }
    } else {
      for (it=this->elements.begin();it!=this->elements.end();++it) {
	Element e = &(it->second);
	if (e->is_leaf()) {
	  for (int i=UHM_ATL;i<UHM_P;++i)
	    assert(e->sparse_pattern(sparse, i));
	}
      }
    }

    fprintf(fp, "i=zeros(%d,1);j=zeros(%d,1);s=zeros(%d,1);\n", 
	    (int)sparse.size(),(int)sparse.size(),(int)sparse.size());

    std::set< std::pair<int,int> >::iterator sit;
    int k=1;
    for (sit=sparse.begin();sit!=sparse.end();++sit) {
      fprintf(fp, "i(%d)=%d;j(%d)=%d;s(%d)=1;\n", 
	      k, sit->first, k, sit->second, k);
      ++k;
    }
    fprintf(fp, "%s = sparse(i,j,s,%d,%d,%d)\n", 
	    ss, n_dofs,n_dofs, (int)sparse.size());
    fprintf(fp, "spy( %s );\n", ss);
    
    // *** close file
    assert(close_file(fp));

    return true;
  } 

  bool Mesh_::export_connectivity(char *full_path, int n_rhs) {
    FILE *fp;

    // *** file open ASCII mode
    assert(open_file(full_path, "w", &fp));

    std::map< int, Element_ >::iterator it;
    int n_elts=0, n_dofs=0;
    for (it=this->elements.begin();it!=this->elements.end();it++) {
      Element e = &(it->second);
      if (e->is_leaf()) ++n_elts;
      n_dofs += e->get_n_dof().first;
    }
    fprintf(fp, "### n_elts, n_dofs, n_rhs\n");
    fprintf(fp, "%d %d %d\n", n_elts, n_dofs, n_rhs);

    for (it=this->elements.begin();it!=this->elements.end();it++) {
      Element e = &(it->second);
      if (e->is_leaf()) {
	for (int i=UHM_ATL;i<UHM_P;++i)
	  assert(e->export_matrix(fp, n_rhs, i));
	for (int i=UHM_BT;i<UHM_RT;++i)
	  assert(e->export_matrix(fp, n_rhs, i));
      }
    }

    // *** close file
    assert(close_file(fp));

    return true;
  }

  bool Mesh_::export_matrix(char *full_path, int n_rhs) {
    FILE *fp;

    // *** file open ASCII mode
    assert(open_file(full_path, "w", &fp));

    std::map< int, Element_ >::iterator it;
    int n_elts=0, n_dofs=0;
    for (it=this->elements.begin();it!=this->elements.end();it++) {
      Element e = &(it->second);
      if (e->is_leaf()) ++n_elts;
      n_dofs += e->get_n_dof().first;
    }

    fprintf(fp, "### n_elts, n_dofs, n_rhs\n");
    fprintf(fp, "%d %d %d\n", n_elts, n_dofs, n_rhs);

    for (it=this->elements.begin();it!=this->elements.end();it++) {
      Element e = &(it->second);
      if (e->is_leaf()) {
	for (int i=UHM_ATL;i<UHM_P;++i)
	  assert(e->get_matrix()->export_matrix(fp, i));
	for (int i=UHM_BT;i<UHM_RT;++i)
	  assert(e->get_matrix()->export_matrix(fp, i));
      }
    }

    // *** close file
    assert(close_file(fp));

    return true;
  }


  bool Mesh_::export_matrix(Sparse sp, int n_rhs) {
    return this->export_matrix(sp, 0, n_rhs);
  }


  bool Mesh_::export_matrix(Sparse sp, int generation, int n_rhs) {

    std::map< int, Element_ >::iterator it;
    for (it=this->elements.begin();it!=this->elements.end();it++) {
      Element e = &(it->second);
      if (e->get_generation() >= generation) 
        sp->import_element(e);
    }

    return true;
  }

  bool Mesh_::import_matrix(Sparse sp, int assemble,
                            std::vector<double> &rhs, int ldb, int n_rhs) {
    sp->reset(n_rhs);

    std::map< int, Element_ >::iterator it;
    for (it=this->elements.begin();it!=this->elements.end();it++) {
      Element e = &(it->second);
      if (e->is_leaf()) {
        sp->export_rhs(rhs, assemble, ldb, n_rhs, UHM_XT, e);
        sp->export_rhs(rhs, assemble, ldb, n_rhs, UHM_XB, e);
      }
    }
    return true;
  }
}
