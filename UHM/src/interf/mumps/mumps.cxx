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

#include "uhm/mesh/element.hxx"

#include "uhm/matrix/uhm/matrix.hxx"

#include "uhm/mesh/mesh.hxx"

#ifdef UHM_INTERF_MUMPS_ENABLE
#include "dmumps_c.h"
#include "zmumps_c.h"
#endif

#include "uhm/interf/sparse.hxx"
#include "uhm/interf/mumps.hxx"

namespace uhm {
  bool Mesh_::export_matrix(Mumps mumps, int n_rhs, int is_sym) {

    // sparse enviroments : fortran index, compressed row;
    int fort=1, cs=0;

    // sparse object
    Sparse sp;

    // easy setting
    if (mumps->is_complex()) 
      sp = new ZSparse_(fort, cs, n_rhs);
    else 
      sp = new DSparse_(fort, cs, n_rhs);
    
    this->export_matrix(sp, n_rhs);
    
    mumps->set_sparse_matrix(sp, is_sym);

    printf("- MUMPS INTERFACED -\n");
    printf("n_dof       : %d\n", (int)mumps->get_n_dof());
    printf("n_nonzero   : %d\n", (int)mumps->get_n_nonzero());
    printf("n_rhs       : %d\n", (int)mumps->get_n_rhs());
    printf("is complex  : %d\n", (int)mumps->is_complex());

    return true;
  }

  bool Mumps_::export_matrix(Mesh m) {
    // sparse enviroments : fortran index, compressed row;
    int fort=1, cs=0;

    // sparse object
    Sparse sp;

    // easy setting
    if (this->is_complex())
      sp = new ZSparse_(fort, cs, n_rhs);
    else
      sp = new DSparse_(fort, cs, n_rhs);

    m->import_matrix(sp, UHM_ASSEMBLED, this->x, this->n_dof, this->n_rhs);

    return true;
  }

  
  // ----------------------------------------------------------------
  // ** mumps control
  bool Mumps_:: is_complex() {
    return (this->datatype == UHM_COMPLEX);
  }

  void Mumps_::set_show_n_rhs(int show_n_rhs) {
    this->show_n_rhs = show_n_rhs;
  }

  void Mumps_::set_par(int par) {
#ifdef UHM_INTERF_MUMPS_ENABLE
    if (is_complex())
      this->zmumps.par = par;
    else
      this->dmumps.par = par;
#endif
  }

  void Mumps_::set_sym(int sym) {
#ifdef UHM_INTERF_MUMPS_ENABLE
    if (is_complex())
      this->zmumps.sym = sym;
    else
      this->dmumps.sym = sym;
#endif
  }

  void Mumps_::set_comm(int comm) {
#ifdef UHM_INTERF_MUMPS_ENABLE
    if (is_complex())
      this->zmumps.comm_fortran = comm;
    else
      this->dmumps.comm_fortran = comm;
#endif
  }

  void Mumps_::set_job(int job) {
#ifdef UHM_INTERF_MUMPS_ENABLE
    if (is_complex())
      this->zmumps.job = job;
    else
      this->dmumps.job = job;
#endif
  }

  void Mumps_::set_icntl(int idx, int val) {
#ifdef UHM_INTERF_MUMPS_ENABLE
    if (is_complex())
      this->zmumps.icntl[idx-1] = val;
    else
      this->dmumps.icntl[idx-1] = val;
#endif
  }

  void Mumps_::set_cntl(int idx, double val) {
#ifdef UHM_INTERF_MUMPS_ENABLE
    if (is_complex())
      this->zmumps.cntl[idx-1] = val;
    else
      this->dmumps.cntl[idx-1] = val;
#endif
  }

  void Mumps_::set_sparse_matrix(Sparse sp, int is_sym)   {
    sp->triplet( is_sym,
                 this->n_dof, this->n_rhs, this->n_nz,
                 this->ia, this->ja,
                 this->a, this->b );

    this->x.clear();
    if (this->is_complex())
      this->x.reserve( this->b.size()*2 );
    else
      this->x.reserve( this->b.size() );

    for (int i=0;i<this->b.size();++i)
      this->x.push_back( this->b.at(i) );

    if (this->show_n_rhs > (int)this->b.size())
      this->show_n_rhs = (int)this->b.size();

#ifdef UHM_INTERF_MUMPS_ENABLE
    if (this->is_complex()) {
      this->zmumps.n    = this->n_dof;
      this->zmumps.nz   = this->n_nz;
      this->zmumps.irn  = (int *)(&(this->ia[0]));
      this->zmumps.jcn  = (int *)(&(this->ja[0]));
      this->zmumps.a    = (mumps_double_complex*)(&(this->a[0]));
      this->zmumps.nrhs = this->n_rhs;
      this->zmumps.lrhs = this->n_dof;
      this->zmumps.rhs  = (mumps_double_complex*)(&(this->x[0]));

    } else {
      this->dmumps.n    = this->n_dof;
      this->dmumps.nz   = this->n_nz;
      this->dmumps.irn  = (int *)(&(this->ia[0]));
      this->dmumps.jcn  = (int *)(&(this->ja[0]));
      this->dmumps.a    = (double*)(&(this->a[0]));
      this->dmumps.nrhs = this->n_rhs;
      this->dmumps.lrhs = this->n_dof;
      this->dmumps.rhs  = (double *)(&(this->x[0]));
    }
#endif

    for (int i = 0; i < this->show_n_rhs; ++i)
      printf("x [ %d ] = % E, b [ %d ] = % E\n",
             i, this->x.at(i),
             i, this->b.at(i) );
  }

  int Mumps_::get_par() {
    int r_val = 0;
#ifdef UHM_INTERF_MUMPS_ENABLE
    if (is_complex())
      r_val = this->zmumps.par;
    else
      r_val = this->dmumps.par;
#endif
    return r_val;
  }

  int Mumps_::get_sym() {
    int r_val = 0;
#ifdef UHM_INTERF_MUMPS_ENABLE
    if (is_complex())
      r_val =  this->zmumps.sym;
    else
      r_val =  this->dmumps.sym;
#endif
    return r_val;
  }

  int Mumps_::get_comm() {
    int r_val = 0;
#ifdef UHM_INTERF_MUMPS_ENABLE
    if (is_complex())
      r_val =  this->zmumps.comm_fortran;
    else
      r_val =  this->dmumps.comm_fortran;
#endif
    return r_val;
  }

  int Mumps_::get_job() {
    int r_val = 0;
#ifdef UHM_INTERF_MUMPS_ENABLE
    if (is_complex())
      r_val =  this->zmumps.job;
    else
      r_val =  this->dmumps.job;
#endif
    return r_val;
  }

  int Mumps_::get_icntl(int idx) {
    int r_val = 0;
#ifdef UHM_INTERF_MUMPS_ENABLE
    if (is_complex())
      r_val =  this->zmumps.icntl[idx-1];
    else
      r_val =  this->dmumps.icntl[idx-1];
#endif
    return r_val;
  }

  double Mumps_::get_cntl(int idx) {
    int r_val = 0;
#ifdef UHM_INTERF_MUMPS_ENABLE
    if (is_complex())
      r_val =  this->zmumps.cntl[idx-1];
    else
      r_val =  this->dmumps.cntl[idx-1];
#endif
    return r_val;
  }

  int Mumps_::get_info(int idx) {
    int r_val = 0;
#ifdef UHM_INTERF_MUMPS_ENABLE
    if (is_complex())
      r_val =  this->zmumps.info[idx-1];
    else
      r_val =  this->dmumps.info[idx-1];
#endif
    return r_val;
  }

  int Mumps_::get_infog(int idx) {
    int r_val = 0;
#ifdef UHM_INTERF_MUMPS_ENABLE
    if (is_complex())
      r_val =  this->zmumps.infog[idx-1];
    else
      r_val =  this->dmumps.infog[idx-1];
#endif
    return r_val;
  }

  double Mumps_::get_rinfo(int idx) {
    int r_val = 0;
#ifdef UHM_INTERF_MUMPS_ENABLE
    if (is_complex())
      r_val =  this->zmumps.rinfo[idx-1];
    else
      r_val =  this->dmumps.rinfo[idx-1];
#endif
    return r_val;
  }

  double Mumps_::get_rinfog(int idx) {
    int r_val = 0;
#ifdef UHM_INTERF_MUMPS_ENABLE
    if (is_complex())
      r_val =  this->zmumps.rinfog[idx-1];
    else
      r_val =  this->dmumps.rinfog[idx-1];
#endif
    return r_val;
  }

  int  Mumps_::get_n_dof()     { return this->n_dof; }
  int  Mumps_::get_n_rhs()     { return this->n_rhs; }
  int  Mumps_::get_n_nonzero() { return this->n_nz; }

  bool Mumps_::validate_sparse_matrix() {
    return true;
  }

  int  Mumps_::run() {
    int ierr = 0;
#ifdef UHM_INTERF_MUMPS_ENABLE
    if (is_complex()) {
      zmumps_c(&(this->zmumps));
    } else {
      dmumps_c(&(this->dmumps));
    }

    ierr = this->get_infog(1);

    if (ierr)
      printf("[MUMPS] ERROR during run: %d\n", ierr);
#endif
    return ierr;
  }

  bool Mumps_::init() {
    int ierr=0;

    this->set_job(-1);
    ierr = this->run();

    if (ierr) {
      printf("[MUMPS] ERROR during init: %d\n", ierr);
      printf("[MUMPS] INIT FAIL INFO[1] %d, INFO[2] %d\n",
             this->get_infog(1), this->get_infog(2));
      return false;
    }
    printf("[MUMPS] Complete init ...\n");
    return true;
  }
  
  bool Mumps_::analyze() {
    int ierr = 0;
    this->set_job(1);

    ierr = this->run();

    if (ierr) {
      printf("[MUMPS] ERROR during analysis: %d\n", ierr);
      printf("[MUMPS] ANALYSIS FAIL INFO[1] %d, INFO[2] %d\n",
             this->get_infog(1), this->get_infog(2));
      return false;
    }
    
    printf( "- Analysis -\n");
    printf( "RINFOG(1)::TOTAL FLOP          %E\n",
            this->get_rinfog(1));
    printf( "INFOG(5) ::MAX FRONT SIZE      %6d\n",
            this->get_infog(5));
    printf( "INFOG(7) ::ORDERING METHOD     %6d\n",
            this->get_infog(7));
    printf( "INFOG(16)::MAX MEMORY FOR CORE %6d [MB]\n",
            this->get_infog(16));
    printf( "INFOG(17)::MAX MEMORY FOR ALL  %6d [MB]\n",
            this->get_infog(17));
    printf( "\n\n");

    printf("[MUMPS] Complete analysis ...\n");
    return true;
  }

  bool Mumps_::decompose() {
    int ierr = 0;
    this->set_job(2);
    ierr = this->run();

    if (ierr) {
      printf("[MUMPS] ERROR during decompose: %d\n", ierr);
      printf("[MUMPS] DECOMPOSITION FAIL INFO[1] %d, INFO[2] %d\n",
             this->get_infog(1), this->get_infog(2));
      return false;
    }

    printf( "- Decomposition -\n");
    printf( "RINFOG(2)::FLOP ASSEM          %E\n",
            this->get_rinfog(2));
    printf( "RINFOG(3)::FLOP DECOMPOSITION  %E\n",
            this->get_rinfog(3));
    printf( "INFOG(11)::MAX FRONT SIZE     %6d\n",
            this->get_infog(11));
    printf( "INFOG(18)::MAX MEMORY FOR CORE %6d [MB]\n",
            this->get_infog(18));
    printf( "INFOG(19)::MAX MEMORY FOR ALL  %6d [MB]\n",
            this->get_infog(19));
    printf( "\n\n");

    printf("[MUMPS] Complete decompose ...\n");
    return true;
  }

  bool Mumps_::solve() {
    int ierr = 0;
    this->set_job(3);
    ierr = this->run();

    if (ierr) {
      printf("[MUMPS] ERROR during solve: %d\n", ierr);
      printf("[MUMPS] SOLUTION FAIL INFO[1] %d, INFO[2] %d\n",
             this->get_infog(1), this->get_infog(2));
      return false;
    }

    printf( "- Solution phase -\n");
    printf( "INFOG(30)::MAX MEMORY FOR CORE %6d [MB]\n",
            this->get_infog(30));
    printf( "INFOG(31)::MAX MEMORY FOR ALL  %6d [MB]\n",
            this->get_infog(31));
    printf( "\n\n");

    printf("[MUMPS] Complete solve ...\n");

    for (int i = 0; i < this->show_n_rhs; ++i)
      printf("x [ %d ] = % E, b [ %d ] = % E\n",
             i, this->x.at(i),
             i, this->b.at(i) );

    return true;
  }

  bool Mumps_::finalize() {
    int ierr = 0;
    this->set_job(-2);
    ierr = this->run();

    if (ierr) {
      printf("[MUMPS] ERROR during finalization: %d\n", ierr);
      printf("[MUMPS] FINALIZATION FAIL INFO[1] %d, INFO[2] %d\n",
             this->get_infog(1), this->get_infog(2));
      return false;
    }
    
    printf("[MUMPS] Complete finalize ...\n");
    return true;
  }
}
