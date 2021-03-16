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
#include "defs.h"
#include "struct.h"
#include "macros.h"
#include "rename.h"

  extern void SetUpGraph(GraphType *, int, int, int, 
			 idxtype *, idxtype *, idxtype *, idxtype *, int);
  extern int idxsum(int, idxtype *);
  extern void InitRandom(int);
  extern void AllocateWorkSpace(CtrlType *, GraphType *, int);
  extern void FreeWorkSpace(CtrlType *, GraphType *);
  extern void MlevelEdgeBisection(CtrlType *, GraphType *, int *, float);
  extern void SplitGraphPart(CtrlType *, 
			     GraphType *, GraphType *, GraphType *);
  extern void GKfree(void **,...); 
}

namespace uhm {
  static bool build_tree_var_4_internal(std::vector< Element > *orphan,
					Mesh m, Element parent, int is_parallel,
					CtrlType *ctrl, GraphType *graph);
  
  // --------------------------------------------------------------
  // ** coarsening 
  bool build_tree_var_4(Mesh m) {
    assert(m && mesh_valid(m));
    { // ** multi level coarsening mesh
      // Start :: collect orphans
      std::vector< Element > *orphan = new std::vector< Element >;
      { // ** preprocess                                              
	// After initial processing of nodes, leaves have separated nodes
	//double t = timer();
	Scheduler_ s;
	s.load(m);
	s.execute_elements_par(&(op_update_connectivity), true);
	//double time = timer() - t;
	//printf("time build_tree_var 4 :: preprocess %lf\n", time);

      }

      {
	//double t = timer();
	std::map< int, Element_ >::iterator mit;
	for (mit=m->elements.begin();mit!=m->elements.end();mit++) 
	  if (mit->second.is_orphan()) 
	    orphan->push_back(&mit->second);

	//double time = timer() - t;
	//printf("time build_tree_var 4 :: adding orphan %lf\n", time);
      }
      // ---------------------------------------------------------------
      // *** bisection of orphan
      if (orphan->size() > 1) {
	Element parent = m->add_element();
	int n = orphan->size(), wgtflag = 1;

	std::vector< int > xadj, adjncy, adjwgt;
	{
	  //double t = timer();

	  orphan_graph(orphan, parent,
		       &xadj, &adjncy, &adjwgt);

	  //double time = timer() - t;
	  //printf("time build_tree_var 4 :: graph orphan %lf\n", time);
	}

	int tvwgt, tpwgts2[2], *label, *where;
	float ubfactor;
	CtrlType ctrl;
	GraphType graph, lgraph, rgraph;
	
	// create graph object to interface to METIS
	{
	  //double t = timer();
	  SetUpGraph(&graph, OP_OEMETIS, n, 1, &xadj[0], &adjncy[0],
		     NULL, &adjwgt[0], wgtflag);
	  //double time = timer() - t;
	  //printf("time build_tree_var 4 :: initial metis %lf\n", time);
	}	  
	// use default control parameter                
	ctrl.CType   = PMETIS_CTYPE;
	ctrl.IType   = PMETIS_ITYPE;
	ctrl.RType   = PMETIS_RTYPE;
	ctrl.dbglvl  = PMETIS_DBGLVL;
	
	ctrl.optype    = OP_PMETIS;
	ctrl.CoarsenTo = 20;
	ctrl.maxvwgt   = 1.5*(idxsum(n, graph.vwgt)/ctrl.CoarsenTo);
 
	InitRandom(-1);

#pragma omp parallel
	{
#pragma omp single nowait
	  {
	    //double t = timer();
	    AllocateWorkSpace(&ctrl, &graph, 2);
	    build_tree_var_4_internal(orphan, m, parent, true, 
				      &ctrl, &graph);
	    FreeWorkSpace(&ctrl, &graph);	
	    //double time = timer() - t;
	    //printf("time build_tree_var 4 :: dissection %lf\n", time);
	  }
	}
      }
      delete orphan;
    }

    { // ** update generation
      //double t = timer();
      Scheduler_ s;
      s.load(m);
      //      s.execute_leaves_seq(&(op_restore_connectivity));
      //      s.execute_elements_par(&(op_update_connectivity), true);
      s.execute_tree(&op_update_generation, true);
      //double time = timer() - t;
      //printf("time build_tree_var 4 :: updating %lf\n", time);
    }
    return true;
  }
  
  bool build_tree_var_4_internal(std::vector< Element > * orphan, 
				 Mesh m, Element parent, int is_parallel,
				 CtrlType *ctrl, GraphType *graph) {

    int n=graph->nvtxs, tvwgt, tpwgts2[2];
    CtrlType  ctrl_elt[2];
    GraphType graph_elt[2];
    
    tvwgt = idxsum(n, graph->vwgt);
    tpwgts2[0] = tvwgt/2;
    tpwgts2[1] = tvwgt-tpwgts2[0];

    MlevelEdgeBisection(ctrl, graph, tpwgts2, 1.0);
    SplitGraphPart(ctrl, graph, &graph_elt[0], &graph_elt[1]);

    GKfree((void**)(&graph->gdata), 
	   (void**)(&graph->rdata), 
	   (void**)(&graph->label), LTERM);

    for (int i=0;i<2;++i) {

      if (graph_elt[i].nvtxs > 3) {
	Element elt = m->add_element();
	
	elt->set_parent( parent );
	parent->add_child( elt );

	ctrl_elt[i] = *ctrl;

	if (graph_elt[i].nvtxs > UHM_DISSECTION_TASK_SIZE && 
	    is_parallel) {
	  
#pragma omp task firstprivate(i, elt, ctrl_elt, graph_elt)
	  {
	    AllocateWorkSpace(&ctrl_elt[i], &graph_elt[i], 2);
	    build_tree_var_4_internal(orphan, m, elt, true,
				      &ctrl_elt[i], &graph_elt[i]);
	    FreeWorkSpace(&ctrl_elt[i], &graph_elt[i]);
	  }
	} else {
	  build_tree_var_4_internal(orphan, m, elt, false,
				    ctrl, &graph_elt[i]);
	}

      } else {
	if ( graph_elt[i].nvtxs ) {
	  Element elt;
	  
	  switch ( graph_elt[i].nvtxs ) {
	  case 1: 
	    elt = parent;
	    break;
	  default:
	    elt = m->add_element();
	    elt->set_parent( parent );
	    parent->add_child( elt );
	    break;
	  }
	  
	  for (int j=0;j<graph_elt[i].nvtxs;++j) {
	    orphan->at( graph_elt[i].label[j] )->set_parent( elt );
	    elt->add_child( orphan->at( graph_elt[i].label[j] ) );
	  }


	  GKfree((void**)(&graph_elt[i].gdata), 
		 (void**)(&graph_elt[i].rdata), 
		 (void**)(&graph_elt[i].label), LTERM);
	}
      }
    }
    /*
    if (is_parallel) {

#pragma omp taskwait

    }
    */
    return true;
  }
}
