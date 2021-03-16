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

#include "dmumps_c.h"
#include "zmumps_c.h"

#include "uhm/interf/mumps.hxx"
#include "uhm/wrapper/mumps_fort.hxx"

//----------------------------------------------------------------------
// FORTRAN INTERFACE
//----------------------------------------------------------------------

// ** MUMPS
//----------------------------------------------------------------------
void UHM_C2F(uhm_mesh_export_matrix_mumps)  ( uhm_fort_p   *mesh,
                                              uhm_fort_p   *mumps,
                                              uhm_fort_int *n_rhs,
                                              uhm_fort_int *is_sym ) {
  uhm::Mesh m = (uhm::Mesh)( *mesh );
  uhm::Mumps mp = (uhm::Mumps)( *mumps );
  m->export_matrix(mp, *n_rhs, *is_sym);
}

void UHM_C2F(uhm_mumps_create)              ( uhm_fort_p   *mumps,
                                              uhm_fort_int *datatype ) {
  uhm::Mumps mp = new uhm::Mumps_( *datatype );
  *mumps = (uhm_fort_p)mp;
}

void UHM_C2F(uhm_mumps_delete)              ( uhm_fort_p   *mumps ) {
  uhm::Mumps mp = (uhm::Mumps)( *mumps );
  delete mp;
}

void UHM_C2F(uhm_mumps_is_complex)          ( uhm_fort_p   *mumps,
                                              uhm_fort_int *flag) {
  uhm::Mumps mp = (uhm::Mumps)( *mumps );
  if (mp->is_complex()) *flag = 1;
  else *flag = 0;
}

void UHM_C2F(uhm_mumps_set_show_n_rhs)      ( uhm_fort_p   *mumps,
                                              uhm_fort_int *show_n_rhs) {
  uhm::Mumps mp = (uhm::Mumps)( *mumps );
  mp->set_show_n_rhs(*show_n_rhs);
}

void UHM_C2F(uhm_mumps_set_par)             ( uhm_fort_p   *mumps,
                                              uhm_fort_int *par) {
  uhm::Mumps mp = (uhm::Mumps)( *mumps );
  mp->set_par(*par);
}

void UHM_C2F(uhm_mumps_set_sym)             ( uhm_fort_p   *mumps,
                                              uhm_fort_int *sym) {
  uhm::Mumps mp = (uhm::Mumps)( *mumps );
  mp->set_sym(*sym);
}

void UHM_C2F(uhm_mumps_set_comm)            ( uhm_fort_p   *mumps,
                                              uhm_fort_int *comm) {
  uhm::Mumps mp = (uhm::Mumps)( *mumps );
  mp->set_comm(*comm);
}

void UHM_C2F(uhm_mumps_set_job)             ( uhm_fort_p   *mumps,
                                              uhm_fort_int *job) {
  uhm::Mumps mp = (uhm::Mumps)( *mumps );
  mp->set_job(*job);
}

void UHM_C2F(uhm_mumps_run)                 ( uhm_fort_p   *mumps ) {
  uhm::Mumps mp = (uhm::Mumps)( *mumps );
  mp->run();
}

void UHM_C2F(uhm_mumps_set_icntl)           ( uhm_fort_p   *mumps,
                                              uhm_fort_int *idx,
                                              uhm_fort_int *val ) {
  uhm::Mumps mp = (uhm::Mumps)( *mumps );
  mp->set_icntl( *idx, *val );
}

void UHM_C2F(uhm_mumps_set_cntl)             ( uhm_fort_p   *mumps,
                                                uhm_fort_int *idx,
                                                uhm_fort_double *val ) {
  uhm::Mumps mp = (uhm::Mumps)( *mumps );
  mp->set_cntl( *idx, *val );
}

void UHM_C2F(uhm_mumps_get_icntl)             ( uhm_fort_p   *mumps,
                                                uhm_fort_int *idx,
                                                uhm_fort_int *val ) {
  uhm::Mumps mp = (uhm::Mumps)( *mumps );
  *val = mp->get_icntl( *idx );
}

void UHM_C2F(uhm_mumps_get_cntl)              ( uhm_fort_p   *mumps,
                                                uhm_fort_int *idx,
                                                uhm_fort_double *val ) {
  uhm::Mumps mp = (uhm::Mumps)( *mumps );
  *val = mp->get_cntl( *idx );
}

void UHM_C2F(uhm_mumps_get_info)              ( uhm_fort_p   *mumps,
                                                uhm_fort_int *idx,
                                                uhm_fort_int *val ) {
  uhm::Mumps mp = (uhm::Mumps)( *mumps );
  *val = mp->get_info( *idx );
}

void UHM_C2F(uhm_mumps_get_infog)             ( uhm_fort_p   *mumps,
                                                uhm_fort_int *idx,
                                                uhm_fort_int *val ) {
  uhm::Mumps mp = (uhm::Mumps)( *mumps );
  *val = mp->get_infog( *idx );
}

void UHM_C2F(uhm_mumps_get_rinfo)             ( uhm_fort_p   *mumps,
                                                uhm_fort_int *idx,
                                                uhm_fort_double *val ) {
  uhm::Mumps mp = (uhm::Mumps)( *mumps );
  *val = mp->get_rinfo( *idx );
}

void UHM_C2F(uhm_mumps_get_rinfog)            ( uhm_fort_p   *mumps,
                                                uhm_fort_int *idx,
                                                uhm_fort_double *val ) {
  uhm::Mumps mp = (uhm::Mumps)( *mumps );
  *val = mp->get_rinfog( *idx );
}



void UHM_C2F(uhm_mumps_init)                ( uhm_fort_p   *mumps ) {
  uhm::Mumps mp = (uhm::Mumps)( *mumps );
  mp->init();
}

void UHM_C2F(uhm_mumps_analyze)             ( uhm_fort_p   *mumps ) {
  uhm::Mumps mp = (uhm::Mumps)( *mumps );
  mp->analyze();
}

void UHM_C2F(uhm_mumps_decompose)           ( uhm_fort_p   *mumps ) {
  uhm::Mumps mp = (uhm::Mumps)( *mumps );
  mp->decompose();
}

void UHM_C2F(uhm_mumps_solve)               ( uhm_fort_p   *mumps ) {
  uhm::Mumps mp = (uhm::Mumps)( *mumps );
  mp->solve();
}

void UHM_C2F(uhm_mumps_export_matrix_uhm)   ( uhm_fort_p   *mumps,
                                              uhm_fort_p   *mesh ) {
  uhm::Mumps mp = (uhm::Mumps)( *mumps );
  uhm::Mesh m = (uhm::Mesh)( *mesh );
  mp->export_matrix(m);
}

void UHM_C2F(uhm_mumps_finalize)            ( uhm_fort_p   *mumps ) {
  uhm::Mumps mp = (uhm::Mumps)( *mumps );
  mp->finalize();
}

