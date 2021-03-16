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

extern "C" {
  extern void METIS_PartGraphRecursive(int *, int *, int *, int *,
				       int *, int *, int *, int *, 
				       int *, int *, int *);
}

namespace uhm {
  static bool build_tree_var_2_internal(Mesh m, int nparts,
                                        std::vector< Element > *orphan,
                                        Element parent);
  // --------------------------------------------------------------
  // ** coarsening 
  bool build_tree_var_2(Mesh m) {
    return build_tree_var_2(m, 2); 
  }
  bool build_tree_var_2(Mesh m, int nparts) {
    assert(m && mesh_valid(m));

    { // ** preprocess
      // After initial processing of nodes, leaves have separated nodes
      Scheduler_ s;
      s.load(m);
      s.execute_leaves_seq(&(op_restore_connectivity));
      s.execute_elements_seq(&(op_update_connectivity), true);
    }

    { // ** multi level coarsening mesh
      // Start :: collect orphans
      std::vector< Element > *orphan = new std::vector< Element >;

      std::map< int, Element_ >::iterator mit;
      for (mit=m->elements.begin();mit!=m->elements.end();mit++) 
	if (mit->second.is_orphan()) 
	  orphan->push_back(&mit->second);

      // ---------------------------------------------------------------
      // *** bisection of orphan
      if (orphan->size() > 1) {
	Element parent = m->add_element();
#pragma omp parallel
	{

#pragma omp single nowait
	  build_tree_var_2_internal(m, nparts, orphan, parent);

	}
      }
      delete orphan;
    }

    { // ** update generation
      Scheduler_ s;
      s.load(m);
      s.execute_leaves_seq(&(op_restore_connectivity));
      s.execute_elements_seq(&(op_update_connectivity), true);
      s.execute_tree(&op_update_generation, true);
    }

    return true;
  }
  
  bool build_tree_var_2_internal(Mesh m, int nparts,
				 std::vector< Element > *orphan, 
				 Element parent) {

    std::vector< std::vector< Element > * > orphan_elt;
    for (int i=0;i<nparts;++i) 
      orphan_elt.push_back(new std::vector< Element >);
    
    { // ** bisection
      int n = orphan->size(), wgtflag = 1, numflag = 0, 
	options = 0, edgecut = 0;
      std::vector< int > xadj, adjncy, adjwgt, part;
      orphan_graph(orphan, parent,
                   &xadj, &adjncy, &adjwgt);
      part.assign    ( n, 0 );

      /*
      std::vector< int > perm, iperm;
      perm.assign    ( n, 0 );
      iperm.assign    ( n, 0 );
      METIS_EdgeND(&n, &xadj[0], &adjncy[0], &numflag, &options,
		   &perm[0], &iperm[0]);

		   
      */
	
      // run METIS
      METIS_PartGraphRecursive(&n, &xadj[0], &adjncy[0], NULL, &adjwgt[0], 
			       &wgtflag, &numflag, 
			       &nparts, &options, &edgecut, &part[0]);

      for (int i=0;i<part.size();++i) 
	orphan_elt[part.at(i)]->push_back( orphan->at(i) );
    }

    { // ** recursion
      for (int i=0;i<nparts;++i) {
	if (orphan_elt[i]->size()>1) {
	  Element elt = m->add_element();
	  elt->set_parent( parent );
	  parent->add_child( elt );

	  if (orphan_elt[i]->size() < nparts) 
	    nparts = orphan_elt[i]->size();
	  
#pragma omp task firstprivate(i, nparts, orphan_elt, elt)
	  build_tree_var_2_internal(m, nparts, orphan_elt[i], elt);
	  
	} else {
	  orphan_elt[i]->at(0)->set_parent( parent );
	  parent->add_child( orphan_elt[i]->at(0) );
	}
      }

#pragma omp taskwait
      
      for (int i=0;i<nparts;++i)
	delete orphan_elt[i];
    }

    return true;
  }
}
