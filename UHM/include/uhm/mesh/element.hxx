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
#ifndef UHM_MESH_ELEMENT_HXX
#define UHM_MESH_ELEMENT_HXX


namespace uhm {
  typedef class Mesh_*    Mesh;
  typedef class Node_*    Node;
  typedef class Element_* Element;
  typedef class Matrix_*  Matrix;
  typedef class Comm_*    Comm;

  bool element_valid(Element e);

  // ----------------------------------------------------------------
  // ** Element class
  class Element_ : public Object_< int > {
  protected:
    // generation of elements
    int generation; 
    
    // tree information :: support only binary tree
    Element parent;
    std::vector< Element > children;
    
    // node, 0 - not separated, 1 factor, 2 schur
    std::map   < Node, int > nodes;

    // node, offset
    std::vector< std::pair<Node, int> > factor, schur;
    
    // matrix objects
    Matrix hm;

    bool   reuse;         // reuse flag
    int    marker[2];     // build_tree_var_2 need marker

    void _init(int id, int gen);
    
  public:
    Element_();
    Element_(int id);
    Element_(int id, int gen);
    virtual ~Element_();

    virtual bool disp();
    virtual bool disp(FILE *stream);

    bool operator<(const Element_ &b) const;
    
    //void inherit_comm(MPI_Comm &comm_in);
    
    // reset operation does not affect to realted objects
    // use reset for global operation
    void reset_parent();
    void reset_children();
    void reset_nodes();
    void reset_factor();
    void reset_schur();
    
    void update_generation();

    void set_matrix(Matrix hm);
    void set_reuse(int flag);
    void set_marker(int index, int marker);
    void set_parent(Element p);

    int  get_generation();
    int  get_height();

    Element get_parent();
    Element get_child(int loc);
    Matrix  get_matrix();
    int     get_marker(int index);

    int  get_n_children();
    int  get_n_nodes();
    int  get_n_factor_nodes();
    int  get_n_schur_nodes();

    std::pair<int,int> get_n_dof();

    void estimate_cost( int method, int datatype, int n_rhs,
                        double &flop_decompose, double &flop_solve,
                        unsigned int &n_nonzero_factor, double &buffer );

    bool is_orphan();
    bool is_leaf();
    bool is_nodes_separated();
    bool is_nodes_arranged();
    bool is_matrix_created();
    bool is_matrix_reusable();

    void collect_leaf_children( int n_max, int &n_leaves, Element *leaves );
    void collect_leaf_children( std::vector< Element > &leaves );

    void add_child (Element c);
    void add_node  (Node n);
    void add_node  (Node n, int separated);
    void add_schur (Node n, int offs);
    void add_factor(Node n, int offs);
        
    void separate_nodes();
    void arrange_nodes();
    void merge_nodes_from_children();
    void merge_nodes(Element c);

    void restore_connectivity();
    void numbering();

    bool write_graphviz_hier(FILE *fp, int is_leaf2root, double max_n_dof);
    bool sparse_pattern(std::set< std::pair<int,int> > &s, int mat);
    bool export_matrix(FILE *stream, int n_rhs, int mat);
    bool export_matrix(int n_rhs, 
		       int &m, int &n, 
		       std::vector< std::pair<int,int> > &ij, 
		       int mat);
    void write_lu_nopiv(FILE *fp, int bmn, int &start, int &end);
    
    friend class Mesh_;
    friend class Helper_;
    friend bool orphan_graph(std::vector< Element > *orphan,
			     Element parent,
			     //
			     std::vector< int > *xadj,
			     std::vector< int > *adjncy,
			     std::vector< int > *adjwgt);

    friend bool build_tree_var_1(Mesh m);

    friend bool element_valid(Element e);
  };

  // ----------------------------------------------------------------
  // ** Definition
  inline void Element_::_init(int id, int gen) {
    //    std::pair<int,int> zero(0,0);
    this->id         = id;
    this->cookie     = UHM_ELEMENT_COOKIE;
    this->generation = gen;
    this->parent     = nil_element;
    this->hm         = nil_matrix;
    this->reuse      = 0;

    for (int i=0;i<2;++i) {
      this->marker[i] = 0;
    }
  }
  inline bool Element_::operator<(const Element_ &b) const { 
    return (this->id < b.id); 
  }
  inline bool element_valid(Element e) { 
    return (e && e->cookie == UHM_ELEMENT_COOKIE); 
  }

}

#endif
