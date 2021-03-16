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
  bool orphan_graph(std::vector< Element > *orphan, 
		    Element parent,
		    //
		    std::vector< int > *xadj, 
		    std::vector< int > *adjncy,
		    std::vector< int > *adjwgt) {

    int n = orphan->size();

    { // ** mark the elements with parent id, and indexing
      double t = timer();
      for (int i=0;i<n;++i) {
	orphan->at(i)->set_marker(0, parent->get_id());
	orphan->at(i)->set_marker(1, i);
      }
      double time = timer()-t;
      //printf("orphan_graph::marking %lf\n", time);
    }

    { // ** create graph
      std::vector< std::map< int, int > > adj;
      adj.assign( n, std::map<int,int>() );
      
      { // **loop through orphan
	double t = timer();

#pragma omp parallel for schedule(static)
	for (int i=0;i<n;++i) {
	  Element elt = orphan->at(i);
	  
	  // loop through schur nodes
	  std::map< Node, int >::iterator nit;
	  for (nit=elt->nodes.begin();nit!=elt->nodes.end();++nit) {
	    
	    if (nit->second == UHM_SEPARATED_SCHUR) {
	      
	      std::set< Element >::iterator eit;
	      for (eit =(nit->first)->owner.begin();
		   eit!=(nit->first)->owner.end();
		   ++eit) {
		  
		// if owner has marker of parent id, add it in adjacency
		if ((*eit)->get_marker(0) == parent->get_id() &&
		    (*eit) != elt) { 
		    
		  // vertex point and its weight
		  std::pair< std::map< int, int >::iterator, bool > ret;
		  std::pair< int,int > in;
		    
		  in.first  = (*eit)->get_marker(1);
		  in.second = nit->first->get_n_dof();
		    
		  ret = adj[i].insert( in );
		    
		  if (ret.second==false) 
		    ret.first->second += in.second;
		}
	      }
	    }
	  }
	}
	double time = timer()-t;
	//printf("orphan_graph::adding %lf\n", time);
      }

      // remove duplicated entries
      {
	double t = timer();

	// clear input array
	if (xadj->size()) xadj->clear();
	if (adjncy->size()) adjncy->clear();
	if (adjwgt->size()) adjwgt->clear();

	xadj->reserve   ( n+1 );
	xadj->push_back ( 0 );
	
	for (int i=0;i<n;++i) {

	  std::map<int,int>::iterator ait;
	  for (ait=adj[i].begin();ait!=adj[i].end();++ait) {
	    adjncy->push_back(ait->first);
	    adjwgt->push_back(ait->second);
	  }
	  xadj->push_back(adjncy->size());
	}
	double time = timer() - t;
	//printf("orphan_graph::map to vector %lf\n", time);
      }
    }
    return true;
  }
}
