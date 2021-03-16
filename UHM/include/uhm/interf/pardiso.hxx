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
#ifndef UHM_INTERF_PARDISO_HXX
#define UHM_INTERF_PARDISO_HXX

namespace uhm {
  typedef class Sparse_*  Sparse;
  typedef class Pardiso_* Pardiso;

  // ----------------------------------------------------------------
  // ** PARDISO
  class Pardiso_ {
  private:
  protected:

    int datatype, show_n_rhs;

    // pardiso control parameter

    // execution step
    // 11 - analysis
    // 12 - analysis + factorization
    // 13 - analysis + factorization + solve + iterative refinement
    // 22 - factorization
    // 23 - factorization + solve + iterative refinement
    // 33 - solve + iterative refinement
    // 0  - release internal memory for LU matrix number MNUM
    // -1 - release all internal memory for all matrices
    int phase;

    // internal data address pointer never touch
    void  *pt[64];

    //  1 - real structrue sym
    //  2 - real sym pos def
    // -2 - real sym indef
    //  3 - complex structrue sym
    //  4 - complex hermitian pos def
    // -4 - complex and hermitian indef
    //  6 - complex and sym
    // 11 - real and nonsym
    // 13 - complex and nonsym
    int    mtype;      

    //  0 - direct , 1 - multi recursive iterative solver
    int    solver;

    // parameters - see manual
    int    iparm[64];
    double dparm[64];
    
    int    maxfct, mnum, msglvl;

    int    num_procs;

    // matrix CSR format
    int n_dof, n_rhs, n_nz;
    std::vector<int> ia, ja;
    std::vector<double> a, b, x;

    void _init(int datatype) {
      assert(datatype == UHM_REAL ||
             datatype == UHM_COMPLEX);
      this->datatype = datatype;
      this->show_n_rhs = 0;
    }
    
  public:

    // methods
    Pardiso_() { _init(UHM_REAL); }
    Pardiso_(int datatype) { _init(datatype); }
    virtual ~Pardiso_() { }

    bool export_matrix(Mesh m);

    void disp();
    void disp( FILE *stream );


    bool is_complex();
    void set_show_n_rhs(int show_n_rhs);
    void set_env(int mtype, int solver, int maxfct, int mnum, int msglvl);
    void set_phase(int phase);
    void set_iparm(int idx, int val);
    void set_dparm(int idx, double val);
    void set_sparse_matrix(Sparse sp, int is_sym);

    int    get_iparm(int idx);
    double get_dparm(int idx);

    int  get_n_dof();    
    int  get_n_rhs();    
    int  get_n_nonzero();

    bool validate_sparse_matrix();
    int  run();

    bool init();
    bool analyze();
    bool decompose();
    bool solve();
    bool finalize();
  };

  // ----------------------------------------------------------------
}

#endif
