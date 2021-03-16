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

namespace uhm {
  // --------------------------------------------------------------
  // ** Mesh
  void Mesh_::set_rhs() {
    // copy b to x for leaf. clear for non-leaf
    std::map< int, Element_ >::iterator eit;
    for (eit=this->elements.begin();eit!=this->elements.end();eit++) {
      Element e = &(eit->second);
      assert(e->is_matrix_created());
      e->get_matrix()->set_rhs( e->is_leaf() );
    }
  }

  void Mesh_::create_matrix_without_buffer(int datatype, int n_rhs) {
    std::map< int, Element_ >::iterator it;
      
    for (it=this->elements.begin();it!=this->elements.end();++it) {
      Element e = &(it->second);
        
      if (e->is_matrix_reusable()) {
        assert(e->is_matrix_created()); 
      } else {
        // create matrix
        // if the matrix is not reusable, which means connectivity is 
        // changed by p or h refinement, delete matrix object and 
        // create new one, which fit to new connectivity
        std::pair<int,int> n_dof = e->get_n_dof();
          
        // TODO :: replace the routine, fix the matrix object
        Matrix hm = new Matrix_FLA_(datatype, 
                                    n_dof.first, n_dof.second, 
                                    n_rhs);
        hm->create_without_buffer();
          
        // if matrix already created, delete that object first
        if (e->is_matrix_created()) delete e->get_matrix();
        e->set_matrix(hm);
      }
    }
  }

  void Mesh_::create_leaf_matrix_buffer() {
    assert(this->get_scheduler()->is_loaded());
    Scheduler s = this->get_scheduler();
    s->execute_leaves_seq( &(op_create_matrix_buffer_with_schur) );
  }
  
  void Mesh_::create_element_matrix_buffer() { 
    this->create_element_matrix_buffer(false); 
  }
  
  void Mesh_::create_element_matrix_buffer(int is_schur) {
    assert(this->get_scheduler()->is_loaded());
    Scheduler s = this->get_scheduler();
    if (is_schur)
      s->execute_elements_seq( &(op_create_matrix_buffer_with_schur), true );
    else
      s->execute_elements_seq( &(op_create_matrix_buffer_without_schur), true );
  }
  
  void Mesh_::create_matrix_buffer() { this->create_matrix_buffer(false); }
  void Mesh_::create_matrix_buffer(int is_schur) {
    std::map< int, Element_ >::iterator it;
    for (it=this->elements.begin();it!=this->elements.end();it++) {
      Element e = &(it->second);
        
      assert(e->is_matrix_created());
        
      // default is create all buffers except for ABR
      // ABR is schur complement 
      for (int i=UHM_ATL;i<UHM_END;i++) 
        if (i!=UHM_ABR)
          e->get_matrix()->create_buffer(i);
        
      // when it has been requested, then allocte for schur
      if (is_schur) 
        e->get_matrix()->create_buffer(UHM_ABR);
    }
      
  }

  void Mesh_::free_matrix() {
    std::map< int, Element_ >::iterator it;
    for (it=this->elements.begin();it!=this->elements.end();it++) {
      Element e = &(it->second);
        
      if (e->is_matrix_created()) {
        delete e->get_matrix();
        e->set_matrix(nil_matrix);
      }
    }
  }

  void Mesh_::free_matrix_buffer() {
    
    std::map< int, Element_ >::iterator it;
    for (it=this->elements.begin();it!=this->elements.end();it++) {
      Element e = &(it->second);
      
      assert(e->is_matrix_created());
      e->get_matrix()->free_buffer();
    }
  }

  void Mesh_::random_matrix()     { this->_random_matrix( false ); }
  void Mesh_::random_spd_matrix() { this->_random_matrix( true ); }
  void Mesh_::triangularize() {
    
    std::map< int, Element_ >::iterator it;
    for (it=this->elements.begin();it!=this->elements.end();it++) {
      Element e = &(it->second);
      if ( e->is_leaf() && !e->is_matrix_reusable() )
        e->get_matrix()->triangularize(FLA_LOWER_TRIANGULAR);
    }
  }

  unsigned int Mesh_::get_n_dof() {
    unsigned int n_dof=0;
    std::map< int, Element_ >::iterator it;
    
    // collect the n_dof of factored nodes
    for (it=this->elements.begin();it!=this->elements.end();++it) {
      Element e = &(it->second);
      n_dof += e->get_n_dof().first;
    }
    return n_dof;
  }
  

  unsigned int Mesh_::get_n_nonzero_factor() {
    unsigned int n_nonzero=0;
    std::map< int, Element_ >::iterator it;
    
    // collect nonzeros which are in leaf UHM
    for (it=this->elements.begin();it!=this->elements.end();it++) {
      Element e = &(it->second);
      std::pair< int, int > n_dof = e->get_n_dof();
      n_nonzero += ( ( n_dof.first * n_dof.first  ) + 
                     ( n_dof.first * n_dof.second )*2 );
    }
    return n_nonzero;
  }

  unsigned int Mesh_::get_n_nonzero() {
    unsigned int n_nonzero=0;
    std::map< int, Element_ >::iterator it;

    // collect nonzeros which are in leaf UHM
    for (it=this->elements.begin();it!=this->elements.end();it++) {
      Element e = &(it->second);
      if (e->is_leaf()) {
        std::pair< int, int > n_dof = e->get_n_dof();
        n_nonzero += ( (n_dof.first + n_dof.second)*
                       (n_dof.first + n_dof.second) );
      }
    }
    return n_nonzero;
  }

  void Mesh_::estimate_cost(int method, int datatype, int n_rhs,
                            double &flop_decompose, 
                            double &flop_solve,
                            unsigned int &n_nonzero_factor,
                            double &buffer) {
    flop_decompose = 0.0;
    flop_solve     = 0.0;
    buffer         = 0.0;

    std::map< int, Element_ >::iterator it;

    // collect nonzeros which are in leaf UHM
    for (it=this->elements.begin();it!=this->elements.end();it++) {
      Element e = &(it->second);
      double 
        elt_decompose, elt_solve, 
        elt_buffer;
      unsigned int elt_n_nonzero_factor;

      e->estimate_cost( method, datatype, n_rhs, 
                        elt_decompose, elt_solve, 
                        elt_n_nonzero_factor, elt_buffer );

      flop_decompose   += elt_decompose;
      flop_solve       += elt_solve;
      n_nonzero_factor += elt_n_nonzero_factor;
      buffer           += elt_buffer;

    }
  }

  double Mesh_::get_residual() {
    double rval = 0.0;
    std::map< int, Element_ >::iterator it;
#ifdef UHM_MULTITHREADING_ENABLE
    // ----------------------------------------------------------
    // ** UHM multi thread
    // ----------------------------------------------------------  
    // Use reduction to collect norm in each elements
    // and should NOT be inside task parallelism

    // change container from map to vector
    std::vector< Element > elts;
    elts.reserve(this->elements.size());
    for (it=this->elements.begin();it!=this->elements.end();it++) {
      Element e = &(it->second);
      if (e->is_matrix_created()) elts.push_back(e);
    }

#pragma omp parallel for reduction(+:rval) schedule(static)
    for (int i=0;i<elts.size();++i) 
      rval += elts.at(i)->get_matrix()->get_residual();

#else
    // ----------------------------------------------------------
    // ** UHM single thread
    // ----------------------------------------------------------  
    for (it=this->elements.begin();it!=this->elements.end();it++) {
      Element e = &(it->second);
      if (e->is_matrix_created()) 
        rval += (e->get_matrix()->get_residual());
    }
#endif

    return rval;
  }
 
  double Mesh_::get_lower_triangular_norm() {
    assert(this->get_scheduler()->is_loaded());
    Scheduler s = this->get_scheduler();

    double rval = 0.0;
    std::map< int, Element_ >::iterator it;
    // ----------------------------------------------------------
    // ** UHM single thread
    // ----------------------------------------------------------  
    for (it=this->elements.begin();it!=this->elements.end();it++) {
      Element e = &(it->second);
      if (e->is_matrix_created()) 
        rval = max(rval, (e->get_matrix()->get_lower_triangular_norm()));
    }

    return rval;
  }

  void Mesh_::_random_matrix( int is_spd ) {

    std::map< int, Element_ >::iterator it;
    for (it=this->elements.begin();it!=this->elements.end();it++) {
      Element e = &(it->second);
      assert(e->is_matrix_created());
      if ( e->is_leaf() && !e->is_matrix_reusable() ) {
        e->get_matrix()->create_buffer();
        if (is_spd)
          e->get_matrix()->random();
        else
          e->get_matrix()->random_spd(FLA_LOWER_TRIANGULAR);
      }
    }
  }
}
