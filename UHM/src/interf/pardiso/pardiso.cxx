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

#include "uhm/interf/sparse.hxx"
#include "uhm/interf/pardiso.hxx"

extern "C" {
  void UHM_C2F(pardisoinit) (void   *, int    *,   int *, int *, double *, int *);
  void UHM_C2F(pardiso)     (void   *, int    *,   int *, int *,    int *, int *,
                             double *, int    *,    int *, int *,   int *, int *,
                             int *, double *, double *, int *, double *);
}


namespace uhm {
  // ----------------------------------------------------------------
  // ** Interface to Pardiso
  bool Mesh_::export_matrix(Pardiso pardiso, int n_rhs, int is_sym) {

    // sparse enviroments : fortran index, compressed row;
    int fort=1, cs=0;

    // environements
    int mtype, solver, maxfct, mnum, msglvl;

    // sparse object
    Sparse sp;

    // easy setting
    if (pardiso->is_complex()) {
      sp = new ZSparse_(fort, cs, n_rhs);

      // here I do not assume Hermitian PD
      if (is_sym)
        mtype  = -4; // complex and hermitian indefinite
      else
        mtype  =  3; // complex and sym structure

      solver =  0; // direct
      maxfct =  1; // do not care
      mnum   =  1; // do not care
      msglvl =  1; // show statistics

    } else {
      sp = new DSparse_(fort, cs, n_rhs);

      // here I do not assume SPD case
      if (is_sym)
        mtype  = -2; // real and sym indefinite
      else
        mtype  =  1; // real and sym structure 
      
      solver =  0; // direct
      maxfct =  1; // do not care
      mnum   =  1; // do not care
      msglvl =  1; // show statistics

    }

    this->export_matrix(sp, n_rhs);

    // default setting
    pardiso->set_env( mtype, solver, maxfct, mnum, msglvl);
    pardiso->set_iparm( 1, 0);
    
    pardiso->set_sparse_matrix(sp, is_sym);

    printf("- PARDISO INTERFACED -\n");
    printf("mtype       : % d\n", mtype);
    printf("-2 real    sym indef, 1 real    sym struct\n");
    printf("-4 complex sym indef, 3 complex sym struct\n");
    printf("\n\n");
    printf("n_dof       : % d\n", (int)pardiso->get_n_dof() );
    printf("n_nonzero   : % d\n", (int)pardiso->get_n_nonzero() );
    printf("n_rhs       : % d\n", (int)pardiso->get_n_rhs() );
    printf("is complex  : % d\n", (int)pardiso->is_complex() );

    return true;
  }

  bool Pardiso_::export_matrix(Mesh m) {
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
  // ** pardiso control
  bool Pardiso_:: is_complex() {
    return (this->datatype == UHM_COMPLEX);
  }

  void Pardiso_::set_show_n_rhs(int show_n_rhs) { 
    this->show_n_rhs = show_n_rhs; 
  }

  void Pardiso_::set_env(int mtype, int solver, int maxfct, int mnum, int msglvl) {
    this->mtype  = mtype;
    this->solver = solver;
    this->maxfct = maxfct;
    this->mnum   = mnum;
    this->msglvl = msglvl;
  }

  void Pardiso_::set_phase(int phase)           { this->phase = phase; }
  void Pardiso_::set_iparm(int idx, int val)    { this->iparm[idx-1] = val; }
  void Pardiso_::set_dparm(int idx, double val) { this->dparm[idx-1] = val; }
  void Pardiso_::set_sparse_matrix(Sparse sp, int is_sym)   {
    sp->compress( is_sym,
                  this->n_dof, this->n_rhs, this->n_nz,
		  this->ia, this->ja,
		  this->a, this->b );
    this->x.clear();
    if (this->is_complex()) 
      this->x.reserve( this->b.size()*2 );
    else
      this->x.reserve( this->b.size() );

    for (int i=0;i<this->b.size();++i)
      this->x.push_back(0.0);

    if (this->show_n_rhs > (int)this->b.size())
      this->show_n_rhs = (int)this->b.size();

    for (int i = 0; i < this->show_n_rhs; ++i) 
      printf("x [ %d ] = % E, b [ %d ] = % E\n", 
	     i, this->x.at(i),
	     i, this->b.at(i) );
  }

  int    Pardiso_::get_iparm(int idx)  { return this->iparm[idx-1]; }
  double Pardiso_::get_dparm(int idx)  { return this->dparm[idx-1]; }

  int  Pardiso_::get_n_dof()     { return this->n_dof; }
  int  Pardiso_::get_n_rhs()     { return this->n_rhs; }
  int  Pardiso_::get_n_nonzero() { return this->n_nz; }

  bool Pardiso_::validate_sparse_matrix() {
    return true;
  }

  int  Pardiso_::run() {
    int idum, ierr=0;
    double ddum;

#ifdef UHM_INTERF_PARDISO_ENABLE
    UHM_C2F(pardiso) (this->pt, &(this->maxfct), &(this->mnum), &(this->mtype),
                      &(this->phase),
                      &(this->n_dof), &(this->a[0]), &(this->ia[0]), &(this->ja[0]),
                      &idum, &(this->n_rhs),
                      this->iparm, &(this->msglvl),
                      &(this->b[0]), &(this->x[0]), &ierr,  this->dparm);
    
    if (ierr !=0) 
      printf("[PARDISO] ERROR during run: %d\n", ierr);
#endif
    return ierr;
  }

  bool Pardiso_::init() {
    int ierr=0;

    // this initialization is must-do
    for (int i=0;i<64;++i) {
      this->pt[i]    = 0;
      this->iparm[i] = 0;
      this->dparm[i] = 0.0;
    }

#ifdef UHM_INTERF_PARDISO_ENABLE
    pardisoinit_ (this->pt,  &this->mtype, &this->solver, 
                 this->iparm, this->dparm, &ierr);
#endif

    switch (ierr) {
    case -10: printf("No license file found \n");break;
    case -11: printf("License is expired \n");break;
    case -12: printf("Wrong username or hostname \n");break;
    default:
      printf("[PARDISO]: License check was successful ... \n");
    }

    return true;
  }

  bool Pardiso_::analyze() {
    this->set_phase(11);

    int ierr = this->run();
    if (ierr != 0) {
      printf("[PARDISO] ERROR during symbolic factorization: %d\n", ierr);
      abort();
    }

    printf("[PARDISO] Reordering completed ... \n");
    printf("[PARDISO] Number of nonzeros in factors  = %d\n", this->get_iparm(18));
    printf("[PARDISO] Number of factorization MFLOPS = %d\n", this->get_iparm(19));

    return true;
  }

  bool Pardiso_::decompose() {
    this->set_phase(22);
    
    int ierr = this->run();
    if (ierr != 0) {
      printf("[PARDISO] ERROR during numerical factorization: %d\n", ierr);
      abort();
    }
    printf("[PARDISO] Factorization completed ...\n");
    
    return true;
  }

  bool Pardiso_::solve() {
    this->set_phase(33);

    int ierr = this->run();
    if (ierr != 0) {
      printf("[PARDISO] ERROR during solution: %d\n", ierr);
      abort();
    }

    printf("[PARDISO] Solve completed ... \n");
    printf("The solution of the system is: \n");

    for (int i = 0; i < this->show_n_rhs; ++i) 
      printf("x [ %d ] = % E, b [ %d ] = % E\n", 
	     i, this->x.at(i),
	     i, this->b.at(i) );

    return true;
  }

  bool Pardiso_::finalize() {
    this->set_phase(-1);
    int ierr = this->run();
    if (ierr != 0) {
      printf("[PARDISO] Error during finalization: %d\n", ierr);
      abort();
    }
    
    printf("[PARDISO] Finalized\n");
    
    return true;
  }

  void Pardiso_::disp() { this->disp(stdout); }
  void Pardiso_::disp( FILE *stream) { 
    
  }

}
