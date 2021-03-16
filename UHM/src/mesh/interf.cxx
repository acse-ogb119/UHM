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

#include "uhm/matrix/uhm/fla.hxx"
#include "uhm/matrix/uhm/helper.hxx"

// for the multi-physics problem, disp = 2. Otherwise, disp = 1
// copy_in does not check buffer in the UHM
namespace uhm {
  static void copy(Mesh mesh, Element elt, int is_in,
		   int datatype, int m, int n,
		   int nod_disp, int *nods, int side, 
		   void *buffer);
  static void copy(Mesh mesh, Element elt, int is_in,
                   std::vector< std::pair<int, int> > &nods, 
		   int side,
                   linal::Flat_ A);

  // -------------------------------------------------------------------
  void Mesh_::copy_in(Element elt,
		      int datatype, int m, int n,
		      int *nods, int side,
		      void *buffer) {
    this->copy_in(elt, datatype, m, n, 
		  UHM_PHYSICS_SINGLE, nods, side,
		  buffer);
  }
  void Mesh_::copy_in(Element elt, 
		      int datatype, int m, int n,
		      int nod_disp, int *nods, int side, 
		      void *buffer) {
    copy(this, elt, true, datatype, m, n, nod_disp, nods, side, buffer);
  }
  void Mesh_::copy_out(Element elt,
		      int datatype, int m, int n,
		      int *nods, int side,
		      void *buffer) {
    this->copy_out(elt, datatype, m, n, 
		   UHM_PHYSICS_SINGLE, nods, side,
		   buffer);
  }
  void Mesh_::copy_out(Element elt, 
		       int datatype, int m, int n,
		       int nod_disp, int *nods, int side, 
		       void *buffer) {
    copy(this, elt, false, datatype, m, n, nod_disp, nods, side, buffer);
  }
  // -------------------------------------------------------------------
  void Mesh_::copy_in(Element elt,
		      std::vector< std::pair<int,int> > &nods, 
		      int side,
		      linal::Flat_ A) {
    copy(this, elt, true, nods, side, A);
  }

  void Mesh_::copy_out(Element elt,
		       std::vector< std::pair<int,int> > &nods, 
		       int side,
		       linal::Flat_ B) {
    copy(this, elt, false, nods, side, B);
  }
  // -------------------------------------------------------------------
  void copy(Mesh mesh, Element e, int is_in,
	    int datatype, int m, int n,
	    int nod_disp, int *nods, int side, 
	    void *buffer) {

    // get the element to interface matrix
    //Element e = mesh->find_element(elt);
    assert(element_valid(e) && e->is_matrix_created());

    // create_buffer is free from double malloc
    e->get_matrix()->create_buffer();

    // phony child
    Element c = new Element_;

    // create child configuration
    int offs = 0;
    for (int i=0;i<e->get_n_nodes();++i) {
      Node nod;
      if (nod_disp == UHM_PHYSICS_MULTI)
	nod = mesh->find_node(std::pair<int,int>(nods[i*2], nods[i*2+1]));
      else 
	nod = mesh->find_node(std::pair<int,int>(nods[i], 0));

      // given node should be mesh object
      if (nod == nil_node) {
	fprintf(stderr, "cannot find out node in the mesh\n");
	abort();
      }

      // child add node into schur container with offset value
      if (nod->get_n_dof()) {
	c->add_schur(nod, offs);
	offs += nod->get_n_dof();
      }
    }

    // user data sanity check
    int n_rhs;
    if (m==0 || n==0) {
      fprintf(stderr, "dimension is null\n");
      abort();
    } else {
      switch (side) {
      case UHM_LHS: assert(m == offs);n_rhs = 1;break;
      case UHM_RHS: assert(m == offs);n_rhs = n;break;
      }
    }

    // create new matrix
    Matrix hm = new Matrix_FLA_(datatype, 0, m, n_rhs); 
    hm->create_without_buffer();

    switch (side) {
    case UHM_LHS: 
      hm->create_buffer(UHM_ABR); 
      break;
    case UHM_RHS: 
      hm->create_buffer(UHM_BB);  
      hm->create_buffer(UHM_XB);  
      break;
    }

    c->set_matrix(hm);

    if (is_in) {
      switch (side) {
      case UHM_LHS: hm->copy_in(UHM_ABR, buffer); break;
      case UHM_RHS: hm->copy_in(UHM_BB,  buffer); break;
      }
    }

    // merge
    Helper_ h(e, c);

    h.set_mapper();
    switch (is_in) {
    case 0:
      if (side == UHM_LHS) { h.branch_ABR(); }
      if (side == UHM_RHS) { h.branch_rhs_x(); }
      break;
    case 1:
      if (side == UHM_LHS) { h.merge_A(); }
      if (side == UHM_RHS) { h.merge_rhs_b(); e->get_matrix()->set_rhs(true); }
      break;
    }

    if (!is_in) {
      switch (side) {
      case UHM_LHS: hm->copy_out(UHM_ABR, buffer); break;
      case UHM_RHS: hm->copy_out(UHM_XB,  buffer); break;
      }
    }
    
    // this will delete associated matrix, too
    delete c;
  }
  void copy(Mesh mesh, Element e, int is_in,
	    std::vector< std::pair<int, int> > &nods, 
	    int side, 
	    linal::Flat_ A) {
    
    assert(element_valid(e) && e->is_matrix_created());
  
    // create_buffer is free from double malloc
    e->get_matrix()->create_buffer();

    // phony child
    Element c = new Element_;

    // create child configuration
    int offs = 0;
    for (int i=0;i<e->get_n_nodes();++i) {
      Node nod;
      nod = mesh->find_node(nods.at(i));

      // given node should be mesh object
      if (nod == nil_node) {
	fprintf(stderr, "cannot find out node in the mesh\n");
	abort();
      }

      // child add node into schur container with offset value
      if (nod->get_n_dof()) {
	c->add_schur(nod, offs);
	offs += nod->get_n_dof();
      }
    }

    // user data sanity check
    int m = A.get_m(), n_rhs;
    switch (side) {
    case UHM_LHS: assert(m == offs);n_rhs = 1;break;
    case UHM_RHS: assert(m == offs);n_rhs = A.get_n();break;
    }

    // create new matrix
    Matrix hm = new Matrix_FLA_(A.get_data_type(), 0, m, n_rhs); 
    hm->create_without_buffer();

    switch (side) {
    case UHM_LHS: 
      hm->create_buffer(UHM_ABR); 
      break;
    case UHM_RHS: 
      hm->create_buffer(UHM_BB);  
      hm->create_buffer(UHM_XB);  
      break;
    }

    c->set_matrix(hm);

    if (is_in) {
      switch (side) {
      case UHM_LHS: hm->copy_in(UHM_ABR, A); break;
      case UHM_RHS: hm->copy_in(UHM_BB,  A); break;
      }
    }

    // merge
    Helper_ h(e, c);

    h.set_mapper();

    switch (is_in) {
    case 0:
      if (side == UHM_LHS) { h.branch_ABR(); }
      if (side == UHM_RHS) { h.branch_rhs_x(); }
      break;
    case 1:
      if (side == UHM_LHS) { h.merge_A(); }
      if (side == UHM_RHS) { h.merge_rhs_b();e->get_matrix()->set_rhs(true); }
      break;
    }

    if (!is_in) {
      switch (side) {
      case UHM_LHS: hm->copy_out(UHM_ABR, A); break;
      case UHM_RHS: hm->copy_out(UHM_XB,  A); break;
      }
    }
    
    // this will delete associated matrix, too
    delete c;
  }
}
