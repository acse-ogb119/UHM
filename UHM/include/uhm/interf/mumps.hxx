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
#ifndef UHM_INTERF_MUMPS_HXX
#define UHM_INTERF_MUMPS_HXX

namespace uhm {
  typedef class  Mumps_*         Mumps;

  // ----------------------------------------------------------------
  // ** MUMPS
  class Mumps_ {
  private:
  protected:
    int datatype, show_n_rhs;

    DMUMPS_STRUC_C dmumps;
    ZMUMPS_STRUC_C zmumps;
    
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
    Mumps_() { _init(UHM_REAL); }
    Mumps_(int datatype) { _init(datatype); }
    virtual ~Mumps_() { }

    bool   export_matrix(Mesh m);

    void   disp();
    void   disp( FILE *stream );

    bool   is_complex();

    void   set_show_n_rhs(int show_n_rhs);

    void   set_par    (int par);
    void   set_sym    (int sym);
    void   set_comm   (int comm);
    void   set_job    (int job);
    
    void   set_icntl  (int idx, int val);
    void   set_cntl   (int idx, double val);

    void   set_sparse_matrix(Sparse sp, int is_sym);

    int    get_par();
    int    get_sym();
    int    get_comm();
    int    get_job();

    int    get_icntl  (int idx);
    double get_cntl   (int idx);

    int    get_info   (int idx);
    int    get_infog  (int idx);
    double get_rinfo  (int idx);
    double get_rinfog (int idx);
    
    int    get_n_dof();
    int    get_n_rhs();
    int    get_n_nonzero();

    bool   validate_sparse_matrix();
    int    run();

    bool   init();
    bool   analyze();
    bool   decompose();
    bool   solve();
    bool   finalize();
  };

  // ----------------------------------------------------------------
}


#endif
