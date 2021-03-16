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
#include "uhm.hxx"

#include "uhm/interf/pardiso.hxx"
#include "uhm/interf/wsmp.hxx"

#include "uhm/wrapper/fort.hxx"

//----------------------------------------------------------------------
// FORTRAN INTERFACE
//----------------------------------------------------------------------

// ** global variable and util
//----------------------------------------------------------------------

void UHM_C2F(uhm_initialize_fla)              () {
  FLA_Init();
}

void UHM_C2F(uhm_finalize_fla)                () {
  FLA_Finalize();
}

void UHM_C2F(uhm_set_num_threads)             ( uhm_fort_int *n_threads ) {
  uhm::set_num_threads( *n_threads );
}

void UHM_C2F(uhm_get_num_threads)             ( uhm_fort_int *n_threads ) {
  *n_threads = uhm::get_num_threads();
}

void UHM_C2F(uhm_set_hier_block_size)         ( uhm_fort_int *blocksize ) {
  uhm::set_hier_block_size( *blocksize );
}
void UHM_C2F(uhm_get_hier_block_size)         ( uhm_fort_int *size ) {
  *size = uhm::get_hier_block_size();
}
void UHM_C2F(uhm_reset_flop)                  () {
  uhm::matrix_reset_flop();
}
void UHM_C2F(uhm_get_flop)                    ( uhm_fort_double *flop ) {
  *flop = uhm::matrix_flop();
}
void UHM_C2F(uhm_get_buffer_used)             ( uhm_fort_double *size ) {
  *size = uhm::matrix_buffer_used();
}
void UHM_C2F(uhm_get_max_buffer_used)         ( uhm_fort_double *size ) {
  *size = uhm::matrix_max_buffer_used();
}
void UHM_C2F(uhm_is_multithreading_enable)    ( uhm_fort_int *flag ) {
  *flag = uhm::is_multithreading_enable();
}
void UHM_C2F(uhm_is_hier_matrix_enable)       ( uhm_fort_int *flag ) {
  *flag = uhm::is_hier_matrix_enable();
}
void UHM_C2F(uhm_timer)                       ( uhm_fort_double *t ) {
  *t = uhm::timer();
}


// ** mesh object
//----------------------------------------------------------------------

void UHM_C2F(uhm_mesh_create)                 ( uhm_fort_p   *mesh ) {
  uhm::Mesh m = new uhm::Mesh_;
  *mesh = (uhm_fort_p)m;
}

void UHM_C2F(uhm_mesh_delete)                 ( uhm_fort_p   *mesh ) {
  uhm::Mesh m = (uhm::Mesh)( *mesh );
  delete m;
}

void UHM_C2F(uhm_add_element)                 ( uhm_fort_p   *mesh,
                                                uhm_fort_p   *elt) {
  uhm::Mesh m = (uhm::Mesh)( *mesh );
  *elt = (uhm_fort_p)m->add_element();
}

void UHM_C2F(uhm_add_node)                   ( uhm_fort_p   *mesh, 
                                               uhm_fort_int *id, 
                                               uhm_fort_int *phy,
                                               uhm_fort_int *n_dof,
                                               uhm_fort_int *p,
                                               uhm_fort_p   *nod ) {
  uhm::Mesh m = (uhm::Mesh)( *mesh );
  std::pair<int,int> in(*id, *phy);
 
  *nod = (uhm_fort_p)m->add_node(in, *n_dof, *p);
}

void UHM_C2F(uhm_element_add_node)            ( uhm_fort_p   *elt,
                                                uhm_fort_p   *nod ) {
  uhm::Element e = (uhm::Element)( *elt );
  uhm::Node    n = (uhm::Node)( *nod );
  e->add_node(n);
}

void UHM_C2F(uhm_copy_in)                     ( uhm_fort_p   *mesh,
                                                uhm_fort_p   *elt,
                                                uhm_fort_int *datatype,
                                                uhm_fort_int *mm,
                                                uhm_fort_int *nn,
                                                uhm_fort_int *node_disp,
                                                uhm_fort_int *nods,
                                                uhm_fort_int *side,
                                                uhm_fort_double *buffer ) {
  uhm::Mesh    m = (uhm::Mesh)( *mesh );
  uhm::Element e = (uhm::Element)( *elt );
  
  m->copy_in( e, *datatype, *mm, *nn, *node_disp, nods, *side, buffer );
}

void UHM_C2F(uhm_copy_out)                    ( uhm_fort_p   *mesh,
                                                uhm_fort_p   *elt,
                                                uhm_fort_int *datatype,
                                                uhm_fort_int *mm,
                                                uhm_fort_int *nn,
                                                uhm_fort_int *node_disp,
                                                uhm_fort_int *nods,
                                                uhm_fort_int *side,
                                                uhm_fort_double *buffer ) {
  uhm::Mesh m    = (uhm::Mesh)( *mesh );
  uhm::Element e = (uhm::Element)( *elt );
  
  m->copy_out( e, *datatype, *mm, *nn, *node_disp, nods, *side, buffer );
}

void UHM_C2F(uhm_mesh_disp)                   ( uhm_fort_p   *mesh ) {
  uhm::Mesh m = (uhm::Mesh)( *mesh );
  m->disp();
}

void UHM_C2F(uhm_get_n_dof)                   ( uhm_fort_p   *mesh,
                                                uhm_fort_int *n_dof ) {
  uhm::Mesh m = (uhm::Mesh)( *mesh );
  *n_dof = m->get_n_dof();
}

void UHM_C2F(uhm_get_n_nonzero)               ( uhm_fort_p   *mesh,
                                                uhm_fort_int *n_nonzero ) {
  uhm::Mesh m = (uhm::Mesh)( *mesh );
  *n_nonzero = m->get_n_nonzero();
}

void UHM_C2F(uhm_estimate_cost)               ( uhm_fort_p      *mesh,
                                                uhm_fort_int    *method,
                                                uhm_fort_int    *datatype,
                                                uhm_fort_int    *n_rhs,
                                                uhm_fort_double *flop_decompose,
                                                uhm_fort_double *flop_solve,
                                                uhm_fort_uint   *n_nonzero_factor,
                                                uhm_fort_double *buffer ) {
  uhm::Mesh m = (uhm::Mesh)( *mesh );
  m->estimate_cost( *method, *datatype, *n_rhs, 
                    *flop_decompose, *flop_solve, 
                    *n_nonzero_factor, *buffer );
}

void UHM_C2F(uhm_mesh_import)                 ( uhm_fort_p    *mesh,
                                                uhm_fort_char *filename) {
  uhm::Mesh m = (uhm::Mesh)( *mesh );
  m->import_file(filename);
}

void UHM_C2F(uhm_build_tree)                  ( uhm_fort_p   *mesh ) {
  uhm::Mesh m = (uhm::Mesh)( *mesh );
  uhm::build_tree(m);
}

void UHM_C2F(uhm_create_matrix_without_buffer)( uhm_fort_p   *mesh,
                                                uhm_fort_int *datatype,
                                                uhm_fort_int *n_rhs ) {
  uhm::Mesh m = (uhm::Mesh)( *mesh );
  m->create_matrix_without_buffer( *datatype, *n_rhs );
}



void UHM_C2F(uhm_create_leaf_matrix_buffer)   ( uhm_fort_p   *mesh ) {
  uhm::Mesh m = (uhm::Mesh)( *mesh );
  m->create_leaf_matrix_buffer();
}

void UHM_C2F(uhm_create_element_matrix_buffer)( uhm_fort_p   *mesh,
                                                uhm_fort_int *is_schur ) {
  uhm::Mesh m = (uhm::Mesh)( *mesh );
  m->create_element_matrix_buffer( *is_schur );
}

void UHM_C2F(uhm_set_rhs)                     ( uhm_fort_p  *mesh ) {
  uhm::Mesh m = (uhm::Mesh)( *mesh );
  m->set_rhs();
}

void UHM_C2F(uhm_random_spd_matrix)           ( uhm_fort_p  *mesh ) {
  uhm::Mesh m = (uhm::Mesh)( *mesh );
  m->random_spd_matrix();
}

void UHM_C2F(uhm_random_matrix)               ( uhm_fort_p  *mesh ) {
  uhm::Mesh m = (uhm::Mesh)( *mesh );
  m->random_matrix();
}
void UHM_C2F(uhm_triangularize)               ( uhm_fort_p  *mesh ) {
  uhm::Mesh m = (uhm::Mesh)( *mesh );
  m->triangularize();
}

void UHM_C2F(uhm_lock)                        ( uhm_fort_p  *mesh ) {
  uhm::Mesh m = (uhm::Mesh)( *mesh );
  m->lock();
}
void UHM_C2F(uhm_unlock)                      ( uhm_fort_p  *mesh ) {
  uhm::Mesh m = (uhm::Mesh)( *mesh );
  m->unlock();
}

// ** decomposition 
//----------------------------------------------------------------------

void UHM_C2F(uhm_chol_with_free)              ( uhm_fort_p  *mesh ) {
  uhm::Mesh m = (uhm::Mesh)( *mesh );
  m->chol_with_free();
}
void UHM_C2F(uhm_lu_nopiv_with_free)          ( uhm_fort_p  *mesh ) {
  uhm::Mesh m = (uhm::Mesh)( *mesh );
  m->lu_nopiv_with_free();
}
void UHM_C2F(uhm_lu_piv_with_free)            ( uhm_fort_p  *mesh ) {
  uhm::Mesh m = (uhm::Mesh)( *mesh );
  m->lu_piv_with_free();
}
void UHM_C2F(uhm_lu_incpiv_with_free)         ( uhm_fort_p  *mesh ) {
  uhm::Mesh m = (uhm::Mesh)( *mesh );
  m->lu_incpiv_with_free();
}
void UHM_C2F(uhm_qr_with_free)                ( uhm_fort_p  *mesh ) {
  uhm::Mesh m = (uhm::Mesh)( *mesh );
  m->qr_with_free();
}

// ** solution
//----------------------------------------------------------------------

void UHM_C2F(uhm_solve_chol)                  ( uhm_fort_p  *mesh ) {
  uhm::Mesh m = (uhm::Mesh)( *mesh );
  m->solve_chol();
}
void UHM_C2F(uhm_solve_lu_nopiv)              ( uhm_fort_p  *mesh ) {
  uhm::Mesh m = (uhm::Mesh)( *mesh );
  m->solve_lu_nopiv();
}
void UHM_C2F(uhm_solve_lu_piv)                ( uhm_fort_p  *mesh ) {
  uhm::Mesh m = (uhm::Mesh)( *mesh );
  m->solve_lu_piv();
}
void UHM_C2F(uhm_solve_qr)                    ( uhm_fort_p  *mesh ) {
  uhm::Mesh m = (uhm::Mesh)( *mesh );
  m->solve_qr();
}

// ** checking
//----------------------------------------------------------------------

void UHM_C2F(uhm_check_chol)                  ( uhm_fort_p  *mesh ) {
  uhm::Mesh m = (uhm::Mesh)( *mesh );
  m->check_chol();
}
void UHM_C2F(uhm_check_lu_nopiv)              ( uhm_fort_p  *mesh ) {
  uhm::Mesh m = (uhm::Mesh)( *mesh );
  m->check_lu_nopiv();
}
void UHM_C2F(uhm_check_lu_piv)                ( uhm_fort_p  *mesh ) {
  uhm::Mesh m = (uhm::Mesh)( *mesh );
  m->check_lu_piv();
}
void UHM_C2F(uhm_check_qr)                    ( uhm_fort_p  *mesh ) {
  uhm::Mesh m = (uhm::Mesh)( *mesh );
  m->check_qr();
}

// ** residual
//----------------------------------------------------------------------

void UHM_C2F(uhm_get_residual)                ( uhm_fort_p  *mesh,
                                                uhm_fort_double *res ) {
  uhm::Mesh m = (uhm::Mesh)( *mesh );
  *res = (uhm_fort_double)m->get_residual();
}


