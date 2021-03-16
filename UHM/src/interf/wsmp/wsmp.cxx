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
#include "uhm/interf/wsmp.hxx"

extern "C" {
  void UHM_C2F(wsetmaxthrds) ( int * );
  void UHM_C2F(wsmp_initialize) ();
  void UHM_C2F(wsmp_clear) ();
  void UHM_C2F(wgffree) ();
  void UHM_C2F(wgsfree)();
  void UHM_C2F(wgsmp) ( int *, int *, int *, double *, double *,
                int *, int *, double *, int *, double * );
}


namespace uhm {
  // ----------------------------------------------------------------
  // ** Interface to WSMP

  bool Mesh_::export_matrix(WSMP wsmp, int n_rhs) {

    // sparse enviroments : fortran index, compressed row;
    int fort=1, cs=0;

    // sparse object
    Sparse sp;

    // easy setting
    if (wsmp->is_complex()) 
      sp = new ZSparse_(fort, cs, n_rhs);
    else 
      sp = new DSparse_(fort, cs, n_rhs);

    this->export_matrix(sp, n_rhs);

    // default setting
    for (int i = 1;i<65;++i) {
      wsmp->set_iparm( i, 0);
      wsmp->set_dparm( i, 0.0);
    }

    wsmp->set_sparse_matrix(sp);

    printf("- WSMP INTERFACED -\n");
    printf("n_dof       : %d\n", (int)wsmp->get_n_dof() );
    printf("n_nonzero   : %d\n", (int)wsmp->get_n_nonzero() );
    printf("n_rhs       : %d\n", (int)wsmp->get_n_rhs() );
    printf("is complex  : %d\n", (int)wsmp->is_complex() );

    delete sp;

    return true;
  }


  bool WSMP_::export_matrix(Mesh m) {
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
  // ** wsmp control
  bool WSMP_:: is_complex() {
    return (this->datatype == UHM_COMPLEX);
  }
  void WSMP_::set_show_n_rhs(int show_n_rhs) { this->show_n_rhs = show_n_rhs; }

  void WSMP_::set_iparm(int idx, int val)    { this->iparm[idx-1] = val; }
  void WSMP_::set_dparm(int idx, double val) { this->dparm[idx-1] = val; }
  void WSMP_::set_sparse_matrix(Sparse sp)   {

    sp->compress( false,
                  this->n_dof, this->n_rhs, this->n_nz,
		  this->ia, this->ja,
		  this->a, this->b );

    this->x.clear();
    this->r.clear();

    if (this->is_complex()) {
      this->x.reserve( this->b.size()*2 );
      this->r.reserve( this->b.size()*2 );

    } else {
      this->x.reserve( this->b.size() );
      this->r.reserve( this->b.size() );

    }

    for (int i=0;i<this->b.size();++i) {
      this->x.push_back( this->b.at(i) );
      this->r.push_back( 0.0 );
    }

    if (this->show_n_rhs > (int)this->b.size())
      this->show_n_rhs = (int)this->b.size();

    for (int i = 0; i < this->show_n_rhs; ++i) 
      printf("x [ %d ] = % E, b [ %d ] = % E, r [ %d ] = % E\n", 
             i, this->x.at(i),
             i, this->b.at(i),
             i, this->r.at(i) );

  }

  int    WSMP_::get_iparm(int idx)  { return this->iparm[idx-1]; }
  double WSMP_::get_dparm(int idx)  { return this->dparm[idx-1]; }

  int  WSMP_::get_n_dof()     { return this->n_dof; }
  int  WSMP_::get_n_rhs()     { return this->n_rhs; }
  int  WSMP_::get_n_nonzero() { return this->n_nz; }

  bool WSMP_::validate_sparse_matrix() {
    return true;
  }

  int  WSMP_::run() {

#ifdef UHM_INTERF_WSMP_ENABLE
    UHM_C2F(wgsmp) (&(this->n_dof), 
                    &(this->ia[0]), &(this->ja[0]), &(this->a[0]), &(this->x[0]),
                    &(this->n_dof), &(this->n_rhs), 
                    &(this->r[0]), this->iparm, this->dparm);
    
    if (this->get_iparm(64) !=0) 
      printf("[WSMP] ERROR during run: %d\n", this->get_iparm(64));
#endif

    return this->get_iparm(64);
  }

  bool WSMP_::init() {
    int n_thread = get_num_threads();

#ifdef UHM_INTERF_WSMP_ENABLE
    UHM_C2F(wsetmaxthrds) ( &n_thread );

    this->set_iparm(1, 0);
    this->set_iparm(2, 0);
    this->set_iparm(3, 0);

    UHM_C2F(wsmp_initialize)();

    this->run();

    if (this->get_iparm(64) !=0) {
      printf("[WSMP] ERROR during initialization: %d\n", this->get_iparm(64));
      abort();
    }
#endif
    return true;
  }

  bool WSMP_::analyze() {
    this->set_iparm(2, 1);
    this->set_iparm(3, 1);
    

    this->run();

    if (this->get_iparm(64) != 0) {
      printf("[WSMP] ERROR during symbolic factorization: %d\n", this->get_iparm(64));
      abort();
    }

    printf("[WSMP] Reordering completed ... \n");
    printf("[WSMP] Number of nonzeros in factors  = %d\n", 1000*this->get_iparm(24));
    printf("[WSMP] Number of factorization MFLOPS = %E\n", this->get_dparm(24));

    return true;
  }

  bool WSMP_::decompose() {
    this->set_iparm(2, 2);
    this->set_iparm(3, 2);

    this->run();

    if (this->get_iparm(64) != 0) {
      printf("[WSMP] ERROR during numerical factorization: %d\n", this->get_iparm(64));
      abort();
    }
    printf("[WSMP] Factorization completed ...\n");
    
    return true;
  }

  bool WSMP_::solve() {
    this->set_iparm(2, 3);
    this->set_iparm(3, 3);

    this->run();

    if (this->get_iparm(64) != 0) {
      printf("[WSMP] ERROR during solution: %d\n", this->get_iparm(64));
      abort();
    }

    printf("[WSMP] Solve completed ... \n");
    printf("The solution of the system is: \n");
    for (int i = 0; i < this->show_n_rhs; ++i) 
      printf("x [ %d ] = % E, b [ %d ] = % E, r [ %d ] = % E\n", 
             i, this->x.at(i),
             i, this->b.at(i),
             i, this->r.at(i) );

    return true;
  }

  bool WSMP_::refine() {
    this->set_iparm(2, 4);
    this->set_iparm(3, 4);

    this->run();

    if (this->get_iparm(64) != 0) {
      printf("[WSMP] ERROR during iterative refinement: %d\n", this->get_iparm(64));
      abort();
    }

    printf("[WSMP] Iterative refinement is completed ... \n");
    printf("The solution of the system is: \n");

    for (int i = 0; i < this->show_n_rhs; ++i) 
      printf("x [ %d ] = % E, b [ %d ] = % E, r [ %d ] = % E\n", 
             i, this->x.at(i),
             i, this->b.at(i),
             i, this->r.at(i) );

    return true;
  }

  bool WSMP_::finalize() {

#ifdef UHM_INTERF_WSMP_ENABLE
    UHM_C2F(wsmp_clear) ();
    UHM_C2F(wgffree) ();
    UHM_C2F(wgsfree) ();

    printf("[WSMP] Finalized\n");
#endif    
    return true;
  }

  void WSMP_::disp() { this->disp(stdout); }
  void WSMP_::disp( FILE *stream) { 
    
  }

}
