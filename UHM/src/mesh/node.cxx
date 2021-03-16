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

namespace uhm {
  // --------------------------------------------------------------
  // ** Node
  Node_::Node_()                       { _init(std::pair<int,int>(0,0),
					       0, 0, 0); }
  Node_::Node_( std::pair<int,int> id) { _init(id, 0,0,0); }
  Node_::Node_( std::pair<int,int> id, 
		int n_dof, 
		int p, int kind)       { _init(id, n_dof, p, kind); }
  Node_::~Node_()                      { }

  // --------------------------------------------------------------
  void Node_::reset_owner() { this->owner.clear(); }

  void Node_::set_kind    (int kind)   { this->kind = kind; }
  void Node_::set_n_dof   (int n_dof)  { this->n_dof = n_dof; }
  void Node_::set_p       (int p)      { this->p = p; }
  void Node_::set_offset  (int offset) { this->offset = offset; }
  void Node_::set_marker  (int marker) { this->marker = marker; }

  int  Node_::get_kind()   { return this->kind; }
  int  Node_::get_n_dof()  { return this->n_dof; }
  int  Node_::get_p()      { return this->p; }
  int  Node_::get_n_owner(){ return this->owner.size(); }
  int  Node_::get_offset() { return this->offset; }
  int  Node_::get_marker() { return this->marker; }

  bool Node_::is_owned_by(Element e) {
    assert(element_valid(e));
    return (this->owner.end() != this->owner.find(e));
  }

  bool Node_::is_same_as(Node n) {
    assert(node_valid(n));

    return ( (this->get_id()    == n->get_id()) &&
	     (this->get_n_dof() == n->get_n_dof()) &&
	     (this->get_kind()  == n->get_kind()) &&
	     (this->get_p()     == n->get_p()) );
  }

  void Node_::add_owner(Element e) {
    assert(element_valid(e));

#pragma omp critical
    this->owner.insert(e);

  }

  void Node_::remove_owner(Element e) {
    assert(element_valid(e));

#pragma omp critical
    this->owner.erase( this->owner.find(e) );

  }

  void Node_::clean_connectivity() {
    // if owner is not leaf, erase it 
    std::set< Element >::iterator it;
#pragma omp critical
    {
      for (it=this->owner.begin();it!=this->owner.end();++it) 
	if (!(*it)->is_leaf()) this->owner.erase(it);
    }
  }

  bool Node_::disp() { return this->disp(stdout); }
  bool Node_::disp(FILE *stream) {
    fprintf(stream, " - Node -\n");
    fprintf(stream, "  id [ %d , %d ] , offset [ %d ], ",
	    this->get_id().first, this->get_id().second,
	    this->get_offset() );
    fprintf(stream, "n_dof [ %d ], kind [ %d ], poly [ %d ]\n", 
	    this->get_n_dof(), this->get_kind(), this->get_p());

    if (this->get_n_owner()) {
      fprintf(stream, "  %d owners [ ", this->get_n_owner());
      
      std::set< Element >::iterator it;
      for (it=this->owner.begin();it!=this->owner.end();it++) 
	fprintf(stream, " %d ", (*it)->get_id());
      fprintf(stream, " ]\n");

    }
    return true;
  }
}
