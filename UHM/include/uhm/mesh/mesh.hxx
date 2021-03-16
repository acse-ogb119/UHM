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
#ifndef UHM_MESH_MESH_HXX
#define UHM_MESH_MESH_HXX

namespace uhm {
  typedef class Node_*      Node;
  typedef class Element_*   Element;
  typedef class Mesh_*      Mesh;
  typedef class Matrix_*    Matrix;
  typedef class Scheduler_* Scheduler;

  typedef class Sparse_*    Sparse;
  typedef class Mumps_*     Mumps;
  typedef class Pardiso_*   Pardiso;
  typedef class WSMP_*      WSMP;
  


  //extern void mesh_new( Mesh &m );
  //extern void mesh_delete( Mesh &m );

  // internal inline functions
  bool   mesh_valid( Mesh m );

  // ----------------------------------------------------------------
  // ** Mesh class
  class Mesh_ : public Object_<int> {
  private:
    int id_element, locker;
  protected:

    std::map < std::pair<int,int>, Node_ > nodes;
    std::map < int, Element_ > elements;

    Scheduler_ scheduler;

    void _init( int id, int id_element );
    void _random_matrix( int is_spd );

  public:
    Mesh_();
    Mesh_( int id );
    Mesh_( int id, int id_element );
    virtual ~Mesh_();

    virtual bool disp();
    virtual bool disp( int mode );
    virtual bool disp( FILE *stream );
    virtual bool disp( FILE *stream, int mode );

    bool import_file(char *full_path);
    bool export_graphviz_hier(char *full_path, int is_leaf2root);
    bool export_sparse_pattern(char *full_path, char *ss, int is_fill_in);
    bool export_connectivity(char *full_path, int n_rhs);
    bool export_matrix(char *full_path, int n_rhs);

    bool export_matrix(Sparse sp, int n_rhs);
    bool export_matrix(Sparse sp, int generation, int n_rhs);
    bool import_matrix(Sparse sp, int assemble,
                       std::vector<double> &rhs, int ldb, int n_rhs);

    // external solver interface
    bool export_matrix(Mumps mumps, int n_rhs, int is_sym);
    bool export_matrix(Pardiso pardiso, int n_rhs, int is_sym);
    bool export_matrix(WSMP wsmp, int n_rhs);

    bool operator<(const Mesh_ &b) const;

    void reset();

    bool is_elements_separated();
    bool is_nodes_numbered();

    void backup( Mesh backup );
    void check_reuse( Mesh backup );

    Scheduler get_scheduler();

    int  get_n_nodes();
    int  get_n_elements();

    void remove_all_nodes();
    void remove_orphan_nodes();

    void remove_all_elements();
    void remove_lower_elements( int gen );

    // single physics
    Node add_node ( int id, int ndof );
    Node add_node ( int id, int ndof, int p );

    // multi physics
    Node add_node ( std::pair<int,int> id, int ndof );
    Node add_node ( std::pair<int,int> id, int ndof, int p );

    // advanced user
    Node add_node ( std::pair<int,int> id, int ndof, int p, int kind );

    Node find_node( int id );
    Node find_node( std::pair<int, int> id );

    bool remove_node( int id );
    bool remove_node( std::pair<int,int> id );

    Element get_root();
    Element add_element();
    Element add_element( int gen );

    Element find_element   ( int id );
    Element insert_element ( int id );
    bool    remove_element ( int id );
    
    void adjust_element_numbering();

    Element refine_element  ( int id, int n_children );
    Element refine_element  ( int id, int is_binary, int n_children );
    Element unrefine_element( int id );

    // solving sequence
    void lock();
    void unlock();
    int  is_locked();

    void create_matrix_without_buffer( int datatype, int n_rhs );

    // i want to separate create buffer between leaf and element
    void create_leaf_matrix_buffer();
    void create_element_matrix_buffer();
    void create_element_matrix_buffer(int is_schur);
    
    // old one same as element buffer creation
    void create_matrix_buffer();
    void create_matrix_buffer( int is_schur );

    void free_matrix();
    void free_matrix_buffer();
 
    void random_matrix();
    void random_spd_matrix();
    void triangularize();

    void copy_in  ( Element elt,
		    int datatype, int m, int n,
		    int *nods, int side,
		    void *buffer );
    void copy_in  ( Element elt,
		    int datatype, int m, int n,
		    int nod_disp, int *nods, int side,
		    void *buffer );
    void copy_in  ( Element elt,
		    std::vector< std::pair<int,int> > &nods, 
		    int side,
		    linal::Flat_ A );

    void copy_out ( Element elt,
		    int datatype, int m, int n,
		    int *nods, int side,
		    void *buffer );
    void copy_out ( Element elt,
		    int datatype, int m, int n,
		    int nod_disp, int *nods, int side,
		    void *buffer );
    void copy_out ( Element elt,
		    std::vector< std::pair<int,int> > &nods, 
		    int side,
		    linal::Flat_ B );

    void         set_rhs();
    double       get_residual();
    double       get_lower_triangular_norm();
    unsigned int get_n_dof();
    unsigned int get_n_nonzero_factor();
    unsigned int get_n_nonzero();

    void   estimate_cost( int method, int datatype, int n_rhs,
                          double &flop_decompose, double &flop_solve,
                          unsigned int &n_nonzero_factor, double &buffer );
    
    void export_graphviz_lu_nopiv(char *full_path, int bmn);


    // ** I am not going to complete all ooc part. 
    //    It is not necessary now.
    
    // ---------------------
    void chol_with_free();
    void chol_without_free();
    void chol_with_ooc();       // not done
    // ---------------------
    void solve_chol_1();
    void solve_chol_2();
    void solve_chol();

    void solve_chol_1_ooc();   // not done
    void solve_chol_2_ooc();   // not done
    void solve_chol_ooc();     // not done

    void check_chol_1();
    void check_chol_2();
    void check_chol();

    void improve_chol();       // not done
    // ---------------------
    void lu_nopiv_with_free();
    void lu_nopiv_without_free();
    void lu_nopiv_with_ooc();   // not done
    // ---------------------
    void solve_lu_nopiv_1();
    void solve_lu_nopiv_2();
    void solve_lu_nopiv();

    void solve_lu_nopiv_1_ooc(); // not done
    void solve_lu_nopiv_2_ooc(); // not done
    void solve_lu_nopiv_ooc();   // not done

    void check_lu_nopiv_1();
    void check_lu_nopiv_2();
    void check_lu_nopiv();

    void improve_lu_nopiv();     // not done
    // ---------------------
    void lu_piv_with_free();
    void lu_piv_without_free();
    void lu_piv_with_ooc();
    // ---------------------
    void lu_incpiv_with_free();
    void lu_incpiv_without_free();
    void lu_incpiv_with_ooc();
    // ---------------------
    void solve_lu_piv_1();
    void solve_lu_piv_2();
    void solve_lu_piv();

    void solve_lu_piv_1_ooc();
    void solve_lu_piv_2_ooc();
    void solve_lu_piv_ooc();

    void check_lu_piv_1();
    void check_lu_piv_2();
    void check_lu_piv();

    void check_lu_piv_1_ooc();
    void check_lu_piv_2_ooc();
    void check_lu_piv_ooc();

    void improve_lu_piv();      // not done
    void improve_lu_piv_ooc();  // not done
    // ---------------------
    void qr_with_free();
    void qr_without_free();
    void qr_with_ooc();       // not done
    // ---------------------
    void solve_qr_1();    // done
    void solve_qr_2();    // done
    void solve_qr();      // done :  not verified

    void check_qr_1();    // not done
    void check_qr_2();    // not done
    void check_qr();      // not done

    void improve_qr();    // not done
    // ---------------------
    friend bool mesh_valid( Mesh m );
    friend bool build_tree_var_1( Mesh m );
    friend bool build_tree_var_2( Mesh m, int nparts );
    friend bool build_tree_var_3( Mesh m );
    friend bool build_tree_var_4( Mesh m );
    friend bool build_tree_var_5( Mesh m );
    friend class Scheduler_;
  };
  // ----------------------------------------------------------------
  // ** Definition
  inline void Mesh_::_init( int id, int id_element ) {
    this->cookie     = UHM_MESH_COOKIE;
    this->id         = id;
    this->id_element = id_element;
  }
  inline bool Mesh_::operator<(const Mesh_ &b) const { 
    return (this->id < b.id); 
  }
  inline bool mesh_valid( Mesh m ) { 
    return (m && m->cookie == UHM_MESH_COOKIE); 
  }

}

#endif
