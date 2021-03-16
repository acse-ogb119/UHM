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
  // * scheduler provide two types of tree traversal
  //   - tree traversal depth first  : execute_tree
  //   - tree traversal breadth first : execute_elements
  //   - tree traversal for leaves : execute_leaves
  // * execute_tree
  //   - post order traversal in tree
  //   - can create dependency between elements in tree
  //   - used int single and multi-thread parallelism using OpenMP 3.0
  //   
  // * execute_elements
  //   - used in elimination tree construction
  //   - used in Supermatrix : supermatrix take care of all dependency
  //   - routine itself does not allow multithread
  //
  // * execute_leaves
  //   - used for the operation only leaves
  //   - used for elimination tree construction
  //   - single thread only

  // --------------------------------------------------------------
  static bool op_tree_seq(int is_leaf2root, Element e, 
                          bool (*op_func)(Element));

  static bool op_leaf_to_root_seq(Element e, bool (*op_func)(Element));
  static bool op_root_to_leaf_seq(Element e, bool (*op_func)(Element));

  static bool op_tree_par(int is_leaf2root, Element e, 
                          bool (*op_func)(Element));

  static bool op_leaf_to_root_par(Element e, bool (*op_func)(Element));
  static bool op_root_to_leaf_par(Element e, bool (*op_func)(Element));
  
  // --------------------------------------------------------------

  static bool op_tree_seq(int is_leaf2root, Element e, 
                          bool(*op_func)(Element)) {
    assert(element_valid(e) && e->is_orphan());

    if (is_leaf2root) assert(op_leaf_to_root_seq(e, op_func));
    else              assert(op_root_to_leaf_seq(e, op_func));

    return true;
  }

  static bool op_leaf_to_root_seq(Element e, bool (*op_func)(Element)) {
    for (int i=0;i<e->get_n_children();i++) {
      Element c = e->get_child(i);
      assert(op_leaf_to_root_seq(c, op_func));
    }
    assert(op_func( e ));
    return true;
  }  
  
  static bool op_root_to_leaf_seq(Element e, bool (*op_func)(Element)) {
    assert(op_func( e ));
    for (int i=0;i<e->get_n_children();i++) {
      Element c = e->get_child(i);
      assert(op_root_to_leaf_seq(c, op_func));
    }
    return true;
  }  

  static bool op_tree_par(int is_leaf2root, Element e, 
                          bool(*op_func)(Element)) {
    assert(element_valid(e) && e->is_orphan());

    if (is_leaf2root) assert(op_leaf_to_root_par(e, op_func));
    else              assert(op_root_to_leaf_par(e, op_func));

    return true;
  }

  static bool op_leaf_to_root_par(Element e, bool (*op_func)(Element)) {
    for (int i=0;i<e->get_n_children();i++) {
      Element c = e->get_child(i);

#pragma omp task firstprivate(c)
      assert(op_leaf_to_root_par(c, op_func));

    }
    
#pragma omp taskwait
    assert(op_func( e ));
    return true;
  }  
  
  static bool op_root_to_leaf_par(Element e, bool (*op_func)(Element)) {
    assert(op_func( e ));

#pragma omp taskwait
    
    for (int i=0;i<e->get_n_children();i++) {
      Element c = e->get_child(i);

#pragma omp task firstprivate(c)
      assert(op_root_to_leaf_par(c, op_func));

    }
    return true;
  }  
  
  // --------------------------------------------------------------
  // ** Scheduler
  Scheduler_::Scheduler_()       { this->_init(0); }
  Scheduler_::Scheduler_(int id) { this->_init(id); }
  Scheduler_::~Scheduler_()      { }
  
  void Scheduler_::load(Mesh m) {
    assert(mesh_valid(m));
    
    { // ** mesh::elements iterator
      std::map< int, Element_ >::iterator mit;
      
      // ** push all elements into container according to its level
      for (mit=m->elements.begin();mit!=m->elements.end();mit++) {
	std::pair<std::map<int, std::vector< Element > >::iterator, bool> ret;
	std::pair< int, std::vector< Element > > dummy;
	dummy.first = mit->second.get_generation();
	ret = this->elements.insert(dummy);

	// ** push element into container
	(ret.first)->second.push_back(&mit->second);
	
	// ** if element is leaf push into leaves container, too
	if ((mit->second).is_leaf())
	  this->leaves.push_back(&mit->second);
      }
    }
  }
  
  void Scheduler_::unload() {
    this->elements.clear(); 
    this->leaves.clear();
  }
  
  int Scheduler_::get_n_elements() {
    int ret = 0;
    
    std::map< int, std::vector< Element > >::reverse_iterator sit;       
    for (sit=this->elements.rbegin();sit!=this->elements.rend();sit++) 
      ret += (int)sit->second.size();
    
    return ret;
  }
  int Scheduler_::get_n_leaves() { return (int)this->leaves.size(); }

  int Scheduler_::is_loaded() { return this->elements.size(); }

  void Scheduler_::get_orphan(std::vector<Element>& orphan) {
    orphan.clear();
    std::map< int, std::vector< Element > >::iterator sit;
    std::vector< Element >::iterator vit;
    for (sit=this->elements.begin();sit!=this->elements.end();sit++) {
      for (vit=sit->second.begin();vit<sit->second.end();vit++) {
        Element e = *vit;
        if (e->is_orphan()) orphan.push_back( e );
      }
    }
  }

  bool Scheduler_::execute_tree(bool (*op_func)(Element), 
				int is_leaf2root) { 
    
    // 1. collect all orphans
    // 
    //    for SMP level there is only one orphan which is root in tree
    std::vector< Element > orphan;
    this->get_orphan(orphan);

    std::map< int, std::vector< Element > >::iterator sit;
    std::vector< Element >::iterator vit;

    /*
    for (sit=this->elements.begin();sit!=this->elements.end();sit++) {
      for (vit=sit->second.begin();vit<sit->second.end();vit++) {
	Element e = *vit;
	if (e->is_orphan()) orphan.push_back( e );
      }
    }
    */

#ifdef UHM_MULTITHREADING_ENABLE    
    // ----------------------------------------------------------             
    // ** UHM multi thread
    // ----------------------------------------------------------  
#pragma omp parallel 
    {
#pragma omp single nowait
      {
	for (vit=orphan.begin();vit<orphan.end();vit++) {
	  Element e = *vit;

#pragma omp task firstprivate(e)
	  assert(op_tree_par(is_leaf2root,  e , op_func ));

	}
      } 
    } // end of parallel region 

#else
    // ----------------------------------------------------------             
    // ** UHM single thread
    // ----------------------------------------------------------  
    for (vit=orphan.begin();vit<orphan.end();vit++) 
      assert(op_tree_seq(is_leaf2root,  *vit , op_func ));
#endif
    
    return true;
  }

  bool Scheduler_::execute_elements_par(bool (*op_func)(Element), 
					int is_leaf2root) { 
    
    // ----------------------------------------------------------             
    // ** UHM multi thread
    // ----------------------------------------------------------  
    if (is_leaf2root) {
      std::map< int, std::vector< Element > >::reverse_iterator sit;
#pragma omp parallel
      {
#pragma omp single nowait
	{
	  for (sit=this->elements.rbegin();sit!=this->elements.rend();sit++) { 
	    std::vector< Element >::iterator vit;
	    for (vit=sit->second.begin();vit<sit->second.end();vit++) {
	      Element e = *vit;
#pragma omp task firstprivate(e)
	      assert(op_func( e ));
	    }
#pragma omp taskwait

	  }
	}
      }
    } else {
      std::map< int, std::vector< Element > >::iterator sit;
#pragma omp parallel
      {
#pragma omp single nowait
	{
	  for (sit=this->elements.begin();sit!=this->elements.end();sit++) { 
	    std::vector< Element >::iterator vit;                            
	    for (vit=sit->second.begin();vit<sit->second.end();vit++) {
	      Element e = *vit;
#pragma omp task firstprivate( e )
	      assert(op_func( e ));
	    }
#pragma omp taskwait

	  }
	}
      }
    }
    return true;
  }

  bool Scheduler_::execute_elements_seq(bool (*op_func)(Element), 
					int is_leaf2root) { 
    
    // ----------------------------------------------------------             
    // ** UHM single thread
    // ----------------------------------------------------------  
    if (is_leaf2root) {
      std::map< int, std::vector< Element > >::reverse_iterator sit;       
      for (sit=this->elements.rbegin();sit!=this->elements.rend();sit++) { 
	std::vector< Element >::iterator vit;
        for (vit=sit->second.begin();vit<sit->second.end();vit++) 
	  assert(op_func( *vit ));
      }
    } else {
      std::map< int, std::vector< Element > >::iterator sit;       
      for (sit=this->elements.begin();sit!=this->elements.end();sit++) { 
	std::vector< Element >::iterator vit;                            
	for (vit=sit->second.begin();vit<sit->second.end();vit++) 
	  assert(op_func( *vit ));
      }
    }
    return true;
  }

  bool Scheduler_::execute_leaves_seq(bool (*op_func)(Element)) {
    // ----------------------------------------------------------             
    // ** UHM single thread
    // ----------------------------------------------------------  
    // ** Do NOT parallelize, it is due to elimination tree construction 
    std::vector< Element >::iterator vit;
    for (vit=this->leaves.begin();vit<this->leaves.end();vit++) {
      assert(op_func(*vit));
    }
    return true;
  }

  bool Scheduler_::execute_leaves_par(bool (*op_func)(Element)) {
    // ----------------------------------------------------------             
    // ** UHM multi thread
    // ----------------------------------------------------------  
    std::vector< Element >::iterator vit;
#pragma omp parallel
    {
#pragma omp single nowait
      {
        for (vit=this->leaves.begin();vit<this->leaves.end();vit++) {
          Element e = *vit;

#pragma omp task firstprivate( e )
          assert(op_func(e));

        }
      }
    }
    return true;
  }
  
  bool Scheduler_::execute_elements(bool (*op_func_1)(Element),
				    bool (*op_func_2)(Element),
				    int is_leaf2root) {
    // Todo :: necessary ????
    // not yet constructed
    return true;
  }

  // for the complex variable case, datasize=16, flop = 2. 
  // for the real varialbe case, datasize=8, flop = 1 
  bool Scheduler_::disp() { return this->disp(8, 1); }
  bool Scheduler_::disp(int datasize, int flop) {
    return this->disp(stdout, datasize, flop); 
  }
  bool Scheduler_::disp(FILE *stream, int datasize, int flop) {
    // short cut
    if (!this->elements.size()) { 
      fprintf(stream, "--------------------------------------------------------\n");
      fprintf(stream, "  scheduler is empty\n");
      return true;
    }

    // mega
    const double mega = 1.0e6;

    // iterator along level < level, container >
    std::map<int, std::vector< Element > >::iterator mit;
    std::map<int, std::vector< Element > >::reverse_iterator rmit;
    mit = this->elements.begin();
    int top = mit->first;

    rmit = this->elements.rbegin();
    int bottom = rmit->first;

    // general information
    fprintf(stream, "top level [ %d ], bottom level [ %d ], diff [ %d ]\n",
	   top, bottom, bottom - top);
    fprintf(stream, "--------------------------------------------------------\n");
    fprintf(stream, "  analyze the problem for the case 'lu_nopiv' \n");
    
    // accumulator
    double all_max_hm = 0.0, all_total_mem_hm_before = 0.0, 
      all_total_mem_hm_after = 0.0, all_total_flop = 0.0;
    
    // loop through level
    for (mit=this->elements.begin();mit!=this->elements.end();mit++) {
      fprintf(stream, "  level [ %d ]\n", mit->first);

      // each level display
      double max_hm = 0.0, min_hm = 0.0,
	total_mem_hm_before = 0.0, total_mem_hm_after = 0.0,
	total_flop = 0.0;

      // iterator for the container
      std::vector< Element >::iterator vit;
      for (vit=mit->second.begin();vit<mit->second.end();vit++) {
	std::pair<int,int> hm_size = (*vit)->get_n_dof();
	double fhm = (double)(hm_size.first);
	double shm = (double)(hm_size.second);
	max_hm = max( max_hm, (fhm+shm) );

	total_mem_hm_before += (pow((fhm+shm),2)*datasize);
	total_mem_hm_after += ((fhm*fhm+fhm*shm*2)*datasize);

	// lu + trsm + trsm + gemm
	total_flop += (2/3*pow(fhm,3) + pow(fhm,2)*fhm + 
		       pow(fhm,2)*shm + 2*shm*shm*fhm)*flop;
      }
      fprintf(stream, "  max hyper matrix [ %6.0lf x %6.0lf , %6.0lf MB ]\n", 
	     max_hm, max_hm, max_hm*max_hm*datasize/mega);
      fprintf(stream, "  total memory storage before decomposition [ %6.0lf MB ]\n", 
	     total_mem_hm_before/mega);
      fprintf(stream, "  total memory storage after decomposition  [ %6.0lf MB ]\n",
	     total_mem_hm_after/mega);
      fprintf(stream, "  total flop [ %6.0lf MFLOP ]\n", total_flop/mega);
      
      all_max_hm = max(all_max_hm, max_hm);
      all_total_mem_hm_before += total_mem_hm_before;
      all_total_mem_hm_after += total_mem_hm_after;
      all_total_flop += total_flop;
    }
    fprintf(stream, "--------------------------------------------------------\n");
    fprintf(stream, "  max hyper matrix [ %6.0lf x %6.0lf , %6.0lf MB ]\n", 
	   all_max_hm, all_max_hm, all_max_hm*all_max_hm*datasize/mega);
    fprintf(stream, "  total memory storage with pre-allocation [ %6.0lf MB ]\n", 
	   all_total_mem_hm_before/mega);
    fprintf(stream, "  total memory storage after decomposition [ %6.0lf MB ]\n",
	   all_total_mem_hm_after/mega);
    fprintf(stream, "  total flop [ %6.0lf MFLOP ]\n", all_total_flop/mega);
    fprintf(stream, "--------------------------------------------------------\n");
    return true;
  }
}
