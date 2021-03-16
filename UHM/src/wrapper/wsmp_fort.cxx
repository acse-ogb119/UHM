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

#include "uhm/interf/wsmp.hxx"
#include "uhm/wrapper/wsmp_fort.hxx"

//----------------------------------------------------------------------
// FORTRAN INTERFACE
//----------------------------------------------------------------------

// ** WSMP
//----------------------------------------------------------------------
void UHM_C2F(uhm_mesh_export_matrix_wsmp)     ( uhm_fort_p   *mesh,
                                                uhm_fort_p   *wsmp,
                                                uhm_fort_int *n_rhs ) {
  uhm::Mesh m = (uhm::Mesh)( *mesh );
  uhm::WSMP w = (uhm::WSMP)( *wsmp );
  m->export_matrix(w, *n_rhs);
}

void UHM_C2F(uhm_wsmp_create)                 ( uhm_fort_p   *wsmp,
                                                uhm_fort_int *datatype ) {
  uhm::WSMP w = new uhm::WSMP_( *datatype );
  *wsmp = (uhm_fort_p)w;
}

void UHM_C2F(uhm_wsmp_delete)                 ( uhm_fort_p   *wsmp ) {
  uhm::WSMP w = (uhm::WSMP)( *wsmp );
  delete w;
}

void UHM_C2F(uhm_wsmp_set_show_n_rhs)         ( uhm_fort_p   *wsmp,
                                                uhm_fort_int *show_n_rhs) {
  uhm::WSMP w = (uhm::WSMP)( *wsmp );
  w->set_show_n_rhs(*show_n_rhs);
}

void UHM_C2F(uhm_wsmp_set_iparm)              ( uhm_fort_p   *wsmp,
                                                uhm_fort_int *idx,
                                                uhm_fort_int *val ) {
  uhm::WSMP w = (uhm::WSMP)( *wsmp );
  w->set_iparm( *idx, *val );
}

void UHM_C2F(uhm_wsmp_set_dparm)              ( uhm_fort_p   *wsmp,
                                                uhm_fort_int *idx,
                                                uhm_fort_double *val ) {
  uhm::WSMP w = (uhm::WSMP)( *wsmp );
  w->set_dparm( *idx, *val );
}

void UHM_C2F(uhm_wsmp_init)                   ( uhm_fort_p   *wsmp ) {
  uhm::WSMP w = (uhm::WSMP)( *wsmp );
  w->init();
}

void UHM_C2F(uhm_wsmp_analyze)                ( uhm_fort_p   *wsmp ) {
  uhm::WSMP w = (uhm::WSMP)( *wsmp );
  w->analyze();
}

void UHM_C2F(uhm_wsmp_decompose)              ( uhm_fort_p   *wsmp ) {
  uhm::WSMP w = (uhm::WSMP)( *wsmp );
  w->decompose();
}

void UHM_C2F(uhm_wsmp_solve)                  ( uhm_fort_p   *wsmp ) {
  uhm::WSMP w = (uhm::WSMP)( *wsmp );
  w->solve();
}

void UHM_C2F(uhm_wsmp_refine)                 ( uhm_fort_p   *wsmp ) {
  uhm::WSMP w = (uhm::WSMP)( *wsmp );
  w->refine();
}

void UHM_C2F(uhm_wsmp_export_matrix_uhm)      ( uhm_fort_p   *wsmp,
                                                uhm_fort_p   *mesh ) {
  uhm::WSMP w = (uhm::WSMP)( *wsmp );
  uhm::Mesh m = (uhm::Mesh)( *mesh );
  w->export_matrix(m);
}

void UHM_C2F(uhm_wsmp_finalize)               ( uhm_fort_p   *wsmp ) {
  uhm::WSMP w = (uhm::WSMP)( *wsmp );
  w->finalize();
}
