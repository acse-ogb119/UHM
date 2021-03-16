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
#ifndef UHM_INTERF_SPARSE_HXX
#define UHM_INTERF_SPARSE_HXX

namespace uhm {
  typedef class Element_* Element;

  typedef class Sparse_*  Sparse;
  typedef class DSparse_* DSparse;
  typedef class ZSparse_* ZSparse;

  // ----------------------------------------------------------------
  // ** Sparse matrix ij triplet and CSR format
  class Sparse_ {
  private:
  protected:
    int fort, datatype, cs, n_dof, n_rhs, axpy;
    void _init(int datatype, int fort, int cs,
               int n_dof, int n_rhs);
      
  public:
    Sparse_() { }
    virtual ~Sparse_() { }

    bool is_complex();
    bool is_fort_numbering();
    bool is_column_major();
    bool is_axpy();

    virtual void set_axpy(int axpy);

    virtual bool import_element(Element e)=0;


    virtual bool import_rhs(std::vector<double> &rhs, int assemble,
                            int ldb, int n_rhs, Element e);
    virtual bool import_rhs(std::vector<double> &rhs, int assemble,
                            int ldb, int n_rhs,
                            int mat, Element e);
    virtual bool export_rhs(std::vector<double> &rhs, int assemble,
                            int ldb, int n_rhs, Element e);
    virtual bool export_rhs(std::vector<double> &rhs, int assemble,
                            int ldb, int n_rhs,
                            int mat, Element e);

        
    virtual void reset(int n_rhs)=0;

    virtual bool triplet(int is_sym, 
                         int &n_dof, int &n_rhs, int &n_nz,
                         std::vector<int> &ia, 
                         std::vector<int> &ja,
                         std::vector<double> &a,
                         std::vector<double> &b)=0;
    
    virtual bool compress(int is_sym,
                          int &n_dof, int &n_rhs, int &n_nz,
                          std::vector<int> &ia, 
                          std::vector<int> &ja,
                          std::vector<double> &a,
                          std::vector<double> &b)=0;

    virtual void disp(FILE *stream, int ab)=0;
    virtual void disp(int ab)=0;
  };

  // ----------------------------------------------------------------
  class DSparse_ : public Sparse_ {
  private:
  protected:
    std::map< std::pair<int,int>, double > a_ij_val, b_ij_val;
  public:
    DSparse_() { }
    DSparse_(int fort, int cs, int n_rhs) { _init(UHM_REAL, fort, cs, 0, n_rhs); }
    virtual ~DSparse_() { }
    
    virtual void reset(int n_rhs);

    virtual bool import_element(Element e);

    virtual bool triplet(int is_sym,
                         int &n_dof, int &n_rhs, int &n_nz,
                         std::vector<int> &ia, 
                         std::vector<int> &ja,
                         std::vector<double> &a,
                         std::vector<double> &b);
    
    virtual bool compress(int is_sym,
                          int &n_dof, int &n_rhs, int &n_nz,
                          std::vector<int> &ia, 
                          std::vector<int> &ja,
                          std::vector<double> &a,
                          std::vector<double> &b);

    virtual void disp(FILE *stream, int ab);
    virtual void disp(int ab);
  };

  // ----------------------------------------------------------------
  class ZSparse_ : public Sparse_ {
  private:
  protected:
    std::map< std::pair<int,int>, std::complex<double> > a_ij_val, b_ij_val;
  public:
    ZSparse_() { }
    ZSparse_(int fort, int cs, int n_rhs) { _init(UHM_COMPLEX, fort, cs, 0, n_rhs); }
    virtual ~ZSparse_() { }
    
    virtual void reset(int n_rhs);
    
    virtual bool import_element(Element e);

    virtual bool triplet(int is_sym,
                         int &n_dof, int &n_rhs, int &n_nz,
                         std::vector<int> &ia, 
                         std::vector<int> &ja,
                         std::vector<double> &a,
                         std::vector<double> &b);
    
    virtual bool compress(int is_sym, 
                          int &n_dof, int &n_rhs, int &n_nz,
                          std::vector<int> &ia, 
                          std::vector<int> &ja,
                          std::vector<double> &a,
                          std::vector<double> &b);

    virtual void disp(FILE *stream, int ab);
    virtual void disp(int ab);
  };
}


#endif
