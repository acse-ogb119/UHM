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

#include "uhm/object.hxx"
#include "uhm/mesh/node.hxx"
#include "uhm/mesh/element.hxx"

#include "uhm/matrix/uhm/matrix.hxx"

#include "uhm/util.hxx"

#include "uhm/interf/sparse.hxx"

namespace uhm {
  // ----------------------------------------------------------------
  // ** Sparse base object
  // This is for bench mark. I neither consider nor expect the performance.
  //
  //
  void Sparse_::_init(int datatype, int fort, int cs, 
                      int n_dof, int n_rhs) {
    assert(datatype==UHM_REAL ||
           datatype==UHM_COMPLEX);
    this->datatype  = datatype;
    this->fort      = fort;
    this->cs        = cs;
    this->n_dof     = n_dof;
    this->n_rhs     = n_rhs;
    this->axpy      = 1;
  }

  bool Sparse_::is_complex()        { return (this->datatype == UHM_COMPLEX); }
  bool Sparse_::is_fort_numbering() { return (bool)(this->fort); }
  bool Sparse_::is_column_major()   { return (bool)(this->cs);   }
  bool Sparse_::is_axpy()           { return (bool)(this->axpy);   }

  void Sparse_::set_axpy(int axpy)  { this->axpy = axpy; }

  bool Sparse_::import_rhs(std::vector<double> &rhs, int assemble,
                           int ldb, int n_rhs, Element e) { 
    assert(e->is_leaf());
    this->import_rhs(rhs, assemble, ldb, n_rhs, UHM_XT, e);
    this->import_rhs(rhs, assemble, ldb, n_rhs, UHM_XB, e);
    return true;
  }

  bool Sparse_::import_rhs(std::vector<double> &rhs, int assemble,
                           int ldb, int n_rhs, 
                           int mat, Element e) { 
    assert(mat == UHM_XT || mat == UHM_XB);

    int m1, m2, n1, n2;
    std::vector<double> val;
    std::vector< std::pair<int,int> > ij;
    
    e->export_matrix(n_rhs, m1, n1, ij, mat);
    e->get_matrix()->export_matrix(m2, n2, val, mat);

    // sanity check
    if ( m1 != m2 || n1 != n2 || m1*n1 != ij.size() ) {
      printf("import_rhs error size mismatch");
      abort();
    }

    if (this->is_complex()) { 
      for (int j=0;j<n1;++j) {
        for (int i=0;i<m1;++i) {
          int re_loc_from = i*2 +     j*2*m1;
          int im_loc_from = i*2 + 1 + j*2*m1;

          int re_loc_to   = ij.at(i+j*m1).first*2 +     ij.at(i+j*m1).second*2*ldb;
          int im_loc_to   = ij.at(i+j*m1).first*2 + 1 + ij.at(i+j*m1).second*2*ldb;

          switch (assemble) {
          case UHM_UNASSEMBLED:
            rhs[ re_loc_to ] += val.at( re_loc_from );
            rhs[ im_loc_to ] += val.at( im_loc_from );
            break;
          case UHM_ASSEMBLED:
            rhs[ re_loc_to ]  = val.at( re_loc_from );
            rhs[ im_loc_to ]  = val.at( im_loc_from );
            break;
          }
        }
      }
    } else {
      for (int j=0;j<n1;++j) {
        for (int i=0;i<m1;++i) {
          int re_loc_from = i + j*m1;
          int re_loc_to   =  ij.at(i+j*m1).first + ij.at(i+j*m1).second*ldb;

          switch (assemble) {
          case UHM_UNASSEMBLED:
            rhs[ re_loc_to ] += val.at( re_loc_from );
            break;
          case UHM_ASSEMBLED:
            rhs[ re_loc_to ]  = val.at( re_loc_from );
            break;
          }
        }
      }
    }
    
    return true;
  }

  bool Sparse_::export_rhs(std::vector<double> &rhs, int assemble, 
                           int ldb, int n_rhs, Element e) {
    assert(e->is_leaf());
    this->export_rhs(rhs, assemble, ldb, n_rhs, UHM_XT, e);
    this->export_rhs(rhs, assemble, ldb, n_rhs, UHM_XB, e);
    return true;
  }

  bool Sparse_::export_rhs(std::vector<double> &rhs, int assemble, 
                           int ldb, int n_rhs,
                           int mat, Element e) {
    assert(mat == UHM_XT || mat == UHM_XB);

    int m, n;
    std::vector<double> val;
    std::vector< std::pair<int,int> > ij;
    
    e->export_matrix(n_rhs, m, n, ij, mat);

    if (this->is_complex()) 
      for (int i=0;i<ij.size();++i) {
        double re =  rhs.at( ij.at(i).first*2 +     ij.at(i).second*2*ldb );
        double im =  rhs.at( ij.at(i).first*2 + 1 + ij.at(i).second*2*ldb );
        val.push_back( re );
        val.push_back( im );
        
        // if unassembled mode, solution is redistribute as unassembled format
        // for the iterative procedure.
        if (assemble == UHM_UNASSEMBLED) {
          val[ re ] = 0.0;
          val[ im ] = 0.0;
        }
      }
    else 
      for (int i=0;i<ij.size();++i) {
        double re =  rhs.at( ij.at(i).first + ij.at(i).second*ldb );
        val.push_back( re );
        
        if (assemble == UHM_UNASSEMBLED)
          val[ re ] = 0.0;
      }
    
    e->get_matrix()->import_matrix(m, n, ldb, val, mat);

    return true;
  }

  // ----------------------------------------------------------------
  // ** Sparse double format
  void DSparse_::reset(int n_rhs) {
    this->n_rhs = n_rhs;
    this->n_dof = 0;
    this->a_ij_val.clear();
    this->b_ij_val.clear();
    this->axpy = 1;
  }

  bool DSparse_::import_element(Element e) {

    // workspace for the interface to UHM
    std::vector< std::pair<int,int> > ij_a, ij_b;
    std::vector< double > a, b;

    // leaf element only and real datatype
    if (e->is_leaf() &&
        !e->get_matrix()->is_complex_datatype()) {
      int m, n;
      for (int i=UHM_ATL;i<UHM_P;++i) {
        if (e->get_matrix()->is_buffer(i)) {
          e->export_matrix(this->n_rhs, m, n, ij_a, i);
          e->get_matrix()->export_matrix(m, n, a, i);
        }
      }
      for (int i=UHM_XT;i<UHM_BT;++i) {
        if (e->get_matrix()->is_buffer(i)) {
          e->export_matrix(this->n_rhs, m, n, ij_b, i);
          e->get_matrix()->export_matrix(m, n, b, i);
        }
      }
    }

    // fortran indexing 
    int shift = this->is_fort_numbering();
    
    // LHS assemble
    for (int i=0;i<ij_a.size();++i) {

      // change the index if it use fortran indexing
      std::pair< int, int > idx;
      if (!this->is_column_major()) {
        idx.first  = ij_a.at(i).first  + shift;
        idx.second = ij_a.at(i).second + shift;
      } else {
        idx.first  = ij_a.at(i).second + shift;
        idx.second = ij_a.at(i).first  + shift;
      }
      
      // insert the entree 
      std::pair< std::map< std::pair<int,int>, double >::iterator, bool > ret;
      ret = this->a_ij_val.insert( std::pair< std::pair<int,int> , double >
                                   (idx, a.at(i)) );
      
      // if exist, add it
      if (ret.second == false && this->is_axpy()) 
        ret.first->second += a.at(i);

    }

    // RHS assemble - column numbering only
    for (int i=0;i<ij_b.size();++i) {

      // change the index if it use fortran indexing
      std::pair< int, int > idx;
      idx.first  = ij_b.at(i).second + shift;
      idx.second = ij_b.at(i).first  + shift;
      
      // insert the entree 
      std::pair< std::map< std::pair<int,int>, double >::iterator, bool > ret;
      ret = this->b_ij_val.insert( std::pair< std::pair<int,int> , double >
                                   (idx, b.at(i)) );
      
      // if exist, add it
      if (ret.second == false && this->is_axpy()) 
        ret.first->second += b.at(i);

    }

    return true;
  }

  bool DSparse_::triplet(int is_sym,
                         int &n_dof, int &n_rhs, int &n_nz,
                         std::vector<int> &ia,
                         std::vector<int> &ja,
                         std::vector<double> &a,
                         std::vector<double> &b) {
    // clear the vector container
    ia.clear();  a.clear();   
    ja.clear();  b.clear();

    std::map< std::pair<int,int>, double >::iterator it;
    if (this->is_column_major()) {

      // ia - entree row number
      // ja - start and end
      if (is_sym)  // symmetric
        for (it=this->a_ij_val.begin();it!=this->a_ij_val.end();++it) {

          // store upper triangular only
          if (it->first.second >= it->first.first) {
            ia.push_back( it->first.second );
            ja.push_back( it->first.first  );
            a.push_back( it->second );
          }

        }
      else         // non-symmetric
        for (it=this->a_ij_val.begin();it!=this->a_ij_val.end();++it) {
          ia.push_back( it->first.second );
          ja.push_back( it->first.first  );
          a.push_back( it->second );
        }
      
      this->n_dof = ja.at( ja.size() - 1 );
      
    } else {
      
      // ia - start and end
      // ja - entree column number
      if (is_sym)  // symmetric
        for (it=this->a_ij_val.begin();it!=this->a_ij_val.end();++it) {
          
          // store upper triangular only
          if (it->first.second >= it->first.first) {
            ia.push_back( it->first.first  );
            ja.push_back( it->first.second );
            a.push_back( it->second );
          }

        }
      else         // non-symmetric
        for (it=this->a_ij_val.begin();it!=this->a_ij_val.end();++it) {
          ia.push_back( it->first.first  );
          ja.push_back( it->first.second );
          a.push_back( it->second );
        }
      
      this->n_dof = ia.at( ia.size() - 1 );

    }

    // columnmajor vector
    for (it=this->b_ij_val.begin();it!=this->b_ij_val.end();++it)
      b.push_back( it->second );

    assert(this->n_rhs == (b.size() / this->n_dof));

    // assign variable
    n_dof = this->n_dof;
    n_rhs = this->n_rhs;
    n_nz  = a.size();

    printf("ndof %d, nrhs %d, nz %d\n", n_dof, n_rhs, n_nz);

    return true;
  }

  bool DSparse_::compress(int is_sym, 
                          int &n_dof, int &n_rhs, int &n_nz,
                          std::vector<int> &ia,
                          std::vector<int> &ja,
                          std::vector<double> &a,
                          std::vector<double> &b) {
    // clear the vector container
    ia.clear();  a.clear();   
    ja.clear();  b.clear();

    std::map< std::pair<int,int>, double >::iterator it;
    int prev, cnt;

    if (this->is_column_major()) {

      // ia - entree row number
      // ja - start and end
      prev = this->a_ij_val.begin()->first.first;

      cnt = this->is_fort_numbering();
      ja.push_back( cnt );

      if (is_sym) { // symmetric
        for (it=this->a_ij_val.begin();it!=this->a_ij_val.end();++it) {
          
          if (it->first.second >= it->first.first) {
            
            if (it->first.first != prev) {
              ja.push_back( cnt );
              prev = it->first.first;
            }
            
            ++cnt;
            ia.push_back( it->first.second );
            a.push_back( it->second );
            
          }
        }
        this->n_dof = ja.size();
        ja.push_back( cnt );
      } else {     // non-symmetric
        for (it=this->a_ij_val.begin();it!=this->a_ij_val.end();++it) {
          if (it->first.first != prev) {
            ja.push_back( cnt );
            prev = it->first.first;
          }
          
          ++cnt;
          ia.push_back( it->first.second );
          a.push_back( it->second );
        }
        this->n_dof = ja.size();
        ja.push_back( cnt );
      }
      
    } else {

      // ia - start and end
      // ja - entree column number
      prev = this->a_ij_val.begin()->first.first;

      cnt = this->is_fort_numbering();
      ia.push_back( cnt );
      if (is_sym) {
        for (it=this->a_ij_val.begin();it!=this->a_ij_val.end();++it) {
          
          if (it->first.second >= it->first.first) {

            if (it->first.first != prev) {
              ia.push_back( cnt );
              prev = it->first.first;
            }
            
            ++cnt;
            ja.push_back( it->first.second );
            a.push_back( it->second );

          }
        }
        this->n_dof = ia.size();
        ia.push_back( cnt );
      } else {     // non-sym
        for (it=this->a_ij_val.begin();it!=this->a_ij_val.end();++it) {
          
          if (it->first.first != prev) {
            ia.push_back( cnt );
            prev = it->first.first;
          }
          
          ++cnt;
          ja.push_back( it->first.second );
          a.push_back( it->second );
        }
        this->n_dof = ia.size();
        ia.push_back( cnt );
      }
    }

    // columnmajor vector
    for (it=this->b_ij_val.begin();it!=this->b_ij_val.end();++it)
      b.push_back( it->second );

    assert(this->n_rhs == (b.size() / this->n_dof));

    // assign variable
    n_dof = this->n_dof;
    n_rhs = this->n_rhs;
    n_nz  = a.size();

    printf("ndof %d, nrhs %d, nz %d\n", n_dof, n_rhs, n_nz);
    
    return true;
  }
                          


  void DSparse_::disp(int ab) {
    this->disp(stdout, ab);
  }
  void DSparse_::disp(FILE *stream, int ab) {
    std::map< std::pair<int,int>, double >::iterator it;

    // fortran index :: matlab format
    int cnt = 1;

    switch (ab) {
    case UHM_LHS:
      fprintf(stream, "i_A = zeros( %d, 1 ); j_A = zeros( %d, 1 ); s_A = zeros( %d , 1);",
              (int)a_ij_val.size(), (int)a_ij_val.size(), (int)a_ij_val.size());
      fprintf(stream, "m_A = %d; n_A = %d;", this->n_dof, this->n_dof);
             
      if (this->is_column_major()) 
        for (it=this->a_ij_val.begin();it!=this->a_ij_val.end();++it) {
          fprintf(stream, "i_A(%d)= %6d; j_A(%d) = %6d; s_A(%d) = % E ; \n", 
                  cnt, it->first.second, cnt, it->first.first, 
                  cnt, it->second);
          ++cnt;
        }
      else 
        for (it=this->a_ij_val.begin();it!=this->a_ij_val.end();++it) {
          fprintf(stream, "i_A(%d)= %6d; j_A(%d)= %6d; s_A(%d)= % E ; \n", 
                  cnt, it->first.first, cnt, it->first.second, 
                  cnt, it->second);
          ++cnt;
        }
      break;
    case UHM_RHS:
      fprintf(stream, "i_B = zeros( %d, 1 ); j_B = zeros( %d, 1 ); s_B = zeros( %d , 1);",
              (int)b_ij_val.size(), (int)b_ij_val.size(), (int)b_ij_val.size());
      fprintf(stream, "m_B = %d; n_B = %d;", this->n_dof, this->n_dof);

      for (it=this->b_ij_val.begin();it!=this->b_ij_val.end();++it) {
        fprintf(stream, "i_B(%d)= %6d; j_B(%d) = %6d; s_B(%d)= % E; \n", 
                cnt, it->first.second, cnt, it->first.first, 
                cnt, it->second);
        ++cnt;
      }
      break;
    }
  }

  // ----------------------------------------------------------------
  // ** Sparse double format
  void ZSparse_::reset(int n_rhs) {
    this->n_rhs = n_rhs;
    this->n_dof = 0;
    this->a_ij_val.clear();
    this->b_ij_val.clear();
    this->axpy = 1;
  }

  bool ZSparse_::import_element(Element e) {

    // increase the number of dofs
    std::pair< int, int > n_dof_elt = e->get_n_dof();
    this->n_dof += n_dof_elt.first;
    
    // workspace for the interface to UHM
    std::vector< std::pair<int,int> > ij_a, ij_b;
    std::vector< double > a, b;

    // leaf element only and real datatype
    if (e->is_leaf() &&
        e->get_matrix()->is_complex_datatype()) {
      int m, n;
      for (int i=UHM_ATL;i<UHM_P;++i) {
        e->export_matrix(this->n_rhs, m, n, ij_a, i);
        e->get_matrix()->export_matrix(m, n, a, i);
      }
      for (int i=UHM_BT;i<UHM_RT;++i) {
        e->export_matrix(n_rhs, m, n, ij_b, i);
        e->get_matrix()->export_matrix(m, n, b, i);
      }
    }

    // fortran indexing 
    int shift = this->is_fort_numbering();
    
    // LHS assemble
    for (int i=0;i<ij_a.size();++i) {

      // change the index if it use fortran indexing
      std::pair< int, int > idx;
      if (!this->is_column_major()) {
        idx.first  = ij_a.at(i).second + shift;
        idx.second = ij_a.at(i).first  + shift;

      } else {
        idx.first  = ij_a.at(i).first  + shift;
        idx.second = ij_a.at(i).second + shift;

      }

      // insert the entree 
      std::pair< std::map< std::pair<int,int>, std::complex<double> >::iterator, bool > ret;
      ret = this->a_ij_val.insert( std::pair< std::pair<int,int> , std::complex<double> >
                                   (idx, std::complex<double>(a.at(2*i), a.at(2*i+1))));
      
      // if exist, add it
      if (ret.second ==false) 
        ret.first->second += a.at(i);

    }
    
    // RHS : assemble
    for (int i=0;i<ij_b.size();++i) {
      
      // change the index if it use fortran indexing
      std::pair< int, int > idx;
      idx.first  = ij_b.at(i).second + shift;
      idx.second = ij_b.at(i).first  + shift;
      
      // insert the entree 
      std::pair< std::map< std::pair<int,int>, std::complex<double> >::iterator, bool > ret;
      ret = this->b_ij_val.insert( std::pair< std::pair<int,int> , std::complex<double> >
                                   (idx, std::complex<double>(b.at(2*i), b.at(2*i+1))));
      
      // if exist, add it
      if (ret.second ==false) 
        ret.first->second += b.at(i);
      
    }
    return true;
  }
  

  bool ZSparse_::triplet(int is_sym,
                         int &n_dof, int &n_rhs, int &n_nz,
                         std::vector<int> &ia,
                         std::vector<int> &ja,
                         std::vector<double> &a,
                         std::vector<double> &b) {
    // clear the vector container
    ia.clear();  a.clear();   
    ja.clear();  b.clear();

    std::map< std::pair<int,int>, std::complex<double> >::iterator it;
    if (this->is_column_major()) {

      // ia - entree row number
      // ja - start and end
      if (is_sym) // symmetric
        for (it=this->a_ij_val.begin();it!=this->a_ij_val.end();++it) {
          
          // store upper triangular only
          if (it->first.second >= it->first.first) {
            ia.push_back( it->first.second );
            ja.push_back( it->first.first  );
          
            a.push_back( it->second.real() );
            a.push_back( it->second.imag() );
          }
        }
      else       // non-symmetric
        for (it=this->a_ij_val.begin();it!=this->a_ij_val.end();++it) {
          
          ia.push_back( it->first.second );
          ja.push_back( it->first.first  );
          
          a.push_back( it->second.real() );
          a.push_back( it->second.imag() );
          
        }

      this->n_dof = ja.at( ja.size() - 1 );
      
    } else {
      
      // ia - start and end
      // ja - entree column number
      if (is_sym)
        for (it=this->a_ij_val.begin();it!=this->a_ij_val.end();++it) {

          if (it->first.second >= it->first.first) {
            ia.push_back( it->first.first  );
            ja.push_back( it->first.second );
            
            a.push_back( it->second.real() );
            a.push_back( it->second.imag() );
          }
        }
      else
        for (it=this->a_ij_val.begin();it!=this->a_ij_val.end();++it) {

          ia.push_back( it->first.first  );
          ja.push_back( it->first.second );
          
          a.push_back( it->second.real() );
          a.push_back( it->second.imag() );
        }

      this->n_dof = ia.at( ia.size() - 1 );

    }

    // columnmajor vector
    for (it=this->b_ij_val.begin();it!=this->b_ij_val.end();++it) {
      b.push_back( it->second.real() );
      b.push_back( it->second.imag() );
    }

    assert(this->n_rhs == (b.size()/2/this->n_dof));

    // assign variable
    n_dof = this->n_dof;
    n_rhs = this->n_rhs;
    n_nz  = a.size() / 2;

    return true;
  }

  bool ZSparse_::compress(int is_sym, 
                          int &n_dof, int &n_rhs, int &n_nz,
                          std::vector<int> &ia,
                          std::vector<int> &ja,
                          std::vector<double> &a,
                          std::vector<double> &b) {
    // clear the vector container
    ia.clear();  a.clear();   
    ja.clear();  b.clear();

    std::map< std::pair<int,int>, std::complex<double> >::iterator it;
    int prev, cnt;

    if (this->is_column_major()) {

      // ia - entree row number
      // ja - start and end
      prev = this->a_ij_val.begin()->first.first;

      cnt = this->is_fort_numbering();
      ja.push_back( cnt );

      if (is_sym) { // symmetric
        for (it=this->a_ij_val.begin();it!=this->a_ij_val.end();++it) {

          if (it->first.second >= it->first.first) {

            if (it->first.first != prev) {
              ja.push_back( cnt );
              prev = it->first.first;
            }

            ++cnt;
            ia.push_back( it->first.second );
            
            a.push_back( it->second.real() );
            a.push_back( it->second.imag() );
          }
        }
        this->n_dof = ja.size();
        ja.push_back( cnt );

      } else {     // non-symmetric

        for (it=this->a_ij_val.begin();it!=this->a_ij_val.end();++it) {

          if (it->first.first != prev) {
            ja.push_back( cnt );
            prev = it->first.first;
          }
          
          ++cnt;
          ia.push_back( it->first.second );
          
          a.push_back( it->second.real() );
          a.push_back( it->second.imag() );
        }
        this->n_dof = ja.size();
        ja.push_back( cnt );

      }
    } else {

      // ia - start and end
      // ja - entree column number
      prev = this->a_ij_val.begin()->first.first;

      cnt = this->is_fort_numbering();
      ia.push_back( cnt );

      if (is_sym) {
        for (it=this->a_ij_val.begin();it!=this->a_ij_val.end();++it) {
         
          if (it->first.second >= it->first.first) {

            if (it->first.first != prev) {
              ia.push_back( cnt );
              prev = it->first.first;
            }
            
            ++cnt;
            ja.push_back( it->first.second );
            a.push_back( it->second.real() );
            a.push_back( it->second.imag() );
          }

        }
        this->n_dof = ia.size();
        ia.push_back( cnt );

      } else {
        for (it=this->a_ij_val.begin();it!=this->a_ij_val.end();++it) {
          
          if (it->first.first != prev) {
            ia.push_back( cnt );
            prev = it->first.first;
          }
          
          ++cnt;
          ja.push_back( it->first.second );
          a.push_back( it->second.real() );
          a.push_back( it->second.imag() );
        }
        this->n_dof = ia.size();
        ia.push_back( cnt );

      }
    }


    // columnmajor vector
    for (it=this->b_ij_val.begin();it!=this->b_ij_val.end();++it) {
      b.push_back( it->second.real() );
      b.push_back( it->second.imag() );
    }

    assert(this->n_rhs == (b.size()/2/this->n_dof));

    n_dof = this->n_dof;
    n_rhs = this->n_rhs;
    n_nz  = a.size()/2;

    return true;
  }


  
  void ZSparse_::disp(int ab) {
    this->disp(stdout, ab);
  }

  void ZSparse_::disp(FILE *stream, int ab) {
    std::map< std::pair<int,int>, std::complex<double> >::iterator it;
    
    switch (ab) {
    case UHM_LHS:
      if (this->is_column_major()) 
        for (it=this->a_ij_val.begin();it!=this->a_ij_val.end();++it) 
          fprintf(stream, "( %6d, %6d )::val ( % E, % E ) \n", 
                  it->first.second, it->first.first, 
                  it->second.real(), it->second.imag() );
      else
        for (it=this->a_ij_val.begin();it!=this->a_ij_val.end();++it) 
          fprintf(stream, "( %6d, %6d )::val ( % E, % E ) \n", 
                  it->first.first, it->first.second, 
                  it->second.real(), it->second.imag() );
      break;
    case UHM_RHS:
      for (it=this->b_ij_val.begin();it!=this->b_ij_val.end();++it) 
        fprintf(stream, "( %6d, %6d )::val ( %8.3lf, %8.3lf ) \n", 
                it->first.second, it->first.first, 
                it->second.real(), it->second.imag() );
      break;
    }
  }
}

