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
#include "uhm/operation/mesh.hxx"
#include "uhm/operation/element.hxx"

#include "uhm/mesh/node.hxx"
#include "uhm/mesh/element.hxx"

#include "uhm/matrix/uhm/matrix.hxx"

#include "uhm/mesh/mesh.hxx"

namespace uhm {
  // --------------------------------------------------------------
  // ** coarsening 
  bool build_tree_var_1(Mesh m) {
    assert(m && mesh_valid(m));

    { // ** preprocess
      uhm::Scheduler_ s;
      s.load(m);
      s.execute_leaves_seq(&(op_restore_connectivity));
      s.execute_elements_seq(&(op_update_connectivity), true);
    }

    { // ** multi level coarsening mesh
      // container for building tree recursively
      std::vector< Element > *orphan = new std::vector< Element >;
      std::vector< Element > *tmp = new std::vector< Element >;

      // branch to connect two orphans
      std::set< std::pair<Element,Element> > branch;

      // branch with weight
      std::multimap< std::pair<int,int>, 
	std::pair<Element, Element> > weighted_branch;

      // ---------------------------------------------------------------
      // 1. collect all orphans for initial iteration
      // mesh iterator
      std::map< int, Element_ >::iterator mit;
      for (mit=m->elements.begin();mit!=m->elements.end();mit++) {
	// if it is orphan, push back into container
	if (mit->second.is_orphan()) orphan->push_back(&mit->second);
      }

      // 2. recursive coarsening
      while (orphan->size() > 1) {
	// 2.1 create branch

	// visit all orphans
	std::vector< Element >::iterator oit;
	for (oit=orphan->begin();oit<orphan->end();oit++) {

	  // visit all schur nodes in the orphan
	  std::map< Node, int >::iterator nit;
	  for (nit=(*oit)->nodes.begin();nit!=(*oit)->nodes.end();nit++) {
	    if (nit->second == UHM_SEPARATED_SCHUR) {
	      Node n;
	      
	      // if node has only two owner, it is branch
	      n = nit->first;
	      std::set< Element >::iterator eit;
	      int ccc = 0;
	      for (eit = n->owner.begin();eit != n->owner.end();eit++) {
		std::pair< Element, Element > neig(*oit,*eit);
		std::pair< int, int > weight(0,0);
		
		// neighbor info is unique key in branch set
		if (neig.first == neig.second) continue;
		if (neig.first > neig.second) std::swap(neig.first, neig.second);
		
		// insert neighbor into branch set
		std::pair< std::set< std::pair< Element, Element > >::iterator, bool> ret;
		ret = branch.insert(neig);
		
		// if it is new branch, calculate weight and insert weighted_branch, too
		if (ret.second) {
		  int i;
		  Element e[2];
		  
		  e[0] = (Element)ret.first->first; 
		  e[1] = (Element)ret.first->second;
		  
		  // iterate nodes in neighbors
		  // weight.first : dof connecting between neighbors
		  // weight.second : dof connecting to others beside neighbor
		  for (i=0;i<2;i++) {
		    std::map< Node, int >::iterator nnit;
		    assert(element_valid(e[i]));
		    for (nnit=e[i]->nodes.begin();nnit!=e[i]->nodes.end();nnit++) {
		      if (nnit->second == UHM_SEPARATED_SCHUR) {
			Node nn;
			nn = nnit->first;
			if (nn->get_n_owner() == 2 && 
			    nn->is_owned_by(e[(i+1)%2])) 
			  weight.first += ( (i) ? (0) : (nn->get_n_dof()));
			else 
			  weight.second -= nn->get_n_dof();
		      }
		    }
		  }

		  // insert into weighted branch - automatically sorted
		  weighted_branch.insert(std::make_pair(weight, neig));
		}
	      }
	    }
	  }
	}
	
	// 2.2 create new parent for the orphans according to its weight
	// reverse iteration big factor, small schur
	std::multimap<std::pair<int,int>, std::pair< Element, Element > >::reverse_iterator wit;
	for (wit=weighted_branch.rbegin();wit!=weighted_branch.rend();wit++) {
	  Element e[2];
	  e[0] = wit->second.first;
	  e[1] = wit->second.second;
	  if (e[0]->is_orphan() && e[1]->is_orphan()) {
	    int i, gen = min(e[0]->get_generation(), e[1]->get_generation());
	    Element p = m->add_element(--gen);
	    for (i=0;i<2;i++) { 
	      p->add_child(e[i]); 
	      e[i]->set_parent(p); 
	    }
	    tmp->push_back(p);
	    
	    p->merge_nodes_from_children();
	    p->separate_nodes();
	  }
	}
	
	// 2.3 add orphan to tmp and swap 
	for (oit=orphan->begin();oit<orphan->end();oit++)
	  if ((*oit)->is_orphan()) tmp->push_back(*oit);

	std::swap(orphan, tmp);
	tmp->clear();
      }
      delete orphan;
      delete tmp;
    }
    return true;
  }
}
