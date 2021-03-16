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
#include "uhm/wrapper/pardiso_fort.hxx"

//----------------------------------------------------------------------
// FORTRAN INTERFACE
//----------------------------------------------------------------------

// ** PARDISO
//----------------------------------------------------------------------
void UHM_C2F(uhm_mesh_export_matrix_pardiso)  ( uhm_fort_p   *mesh,
                                                uhm_fort_p   *pardiso,
                                                uhm_fort_int *n_rhs,
                                                uhm_fort_int *is_sym) {
  uhm::Mesh m = (uhm::Mesh)( *mesh );
  uhm::Pardiso p = (uhm::Pardiso)( *pardiso );
  m->export_matrix(p, *n_rhs, *is_sym);
}

void UHM_C2F(uhm_pardiso_create)              ( uhm_fort_p   *pardiso,
                                                uhm_fort_int *datatype ) {
  uhm::Pardiso p = new uhm::Pardiso_( *datatype );
  *pardiso = (uhm_fort_p)p;
}

void UHM_C2F(uhm_pardiso_delete)              ( uhm_fort_p   *pardiso ) {
  uhm::Pardiso p = (uhm::Pardiso)( *pardiso );
  delete p;
}

void UHM_C2F(uhm_pardiso_is_complex)          ( uhm_fort_p   *pardiso,
                                                uhm_fort_int *flag) {
  uhm::Pardiso p = (uhm::Pardiso)( *pardiso );
  if (p->is_complex()) *flag = 1;
  else *flag = 0;
}

void UHM_C2F(uhm_pardiso_set_show_n_rhs)      ( uhm_fort_p   *pardiso,
                                                uhm_fort_int *show_n_rhs) {
  uhm::Pardiso p = (uhm::Pardiso)( *pardiso );
  p->set_show_n_rhs(*show_n_rhs);
}

void UHM_C2F(uhm_pardiso_set_phase)           ( uhm_fort_p   *pardiso,
                                                uhm_fort_int *phase) {
  uhm::Pardiso p = (uhm::Pardiso)( *pardiso );
  p->set_phase(*phase);
}

void UHM_C2F(uhm_pardiso_run)                 ( uhm_fort_p   *pardiso ) {
  uhm::Pardiso p = (uhm::Pardiso)( *pardiso );
  p->run();
}
                                                
void UHM_C2F(uhm_pardiso_set_iparm)           ( uhm_fort_p   *pardiso,
                                                uhm_fort_int *idx,
                                                uhm_fort_int *val ) {
  uhm::Pardiso p = (uhm::Pardiso)( *pardiso );
  p->set_iparm( *idx, *val );
}

void UHM_C2F(uhm_pardiso_set_dparm)           ( uhm_fort_p   *pardiso,
                                                uhm_fort_int *idx,
                                                uhm_fort_double *val ) {
  uhm::Pardiso p = (uhm::Pardiso)( *pardiso );
  p->set_dparm( *idx, *val );
}

void UHM_C2F(uhm_pardiso_init)                ( uhm_fort_p   *pardiso ) {
  uhm::Pardiso p = (uhm::Pardiso)( *pardiso );
  p->init();
}

void UHM_C2F(uhm_pardiso_analyze)             ( uhm_fort_p   *pardiso ) {
  uhm::Pardiso p = (uhm::Pardiso)( *pardiso );
  p->analyze();
}

void UHM_C2F(uhm_pardiso_decompose)           ( uhm_fort_p   *pardiso ) {
  uhm::Pardiso p = (uhm::Pardiso)( *pardiso );
  p->decompose();
}

void UHM_C2F(uhm_pardiso_solve)               ( uhm_fort_p   *pardiso ) {
  uhm::Pardiso p = (uhm::Pardiso)( *pardiso );
  p->solve();
}

void UHM_C2F(uhm_pardiso_export_matrix_uhm)   ( uhm_fort_p   *pardiso,
                                                uhm_fort_p   *mesh ) {
  uhm::Pardiso p = (uhm::Pardiso)( *pardiso );
  uhm::Mesh m = (uhm::Mesh)( *mesh );
  p->export_matrix(m);
}

void UHM_C2F(uhm_pardiso_finalize)            ( uhm_fort_p   *pardiso ) {
  uhm::Pardiso p = (uhm::Pardiso)( *pardiso );
  p->finalize();
}

