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
#ifndef UHM_MESH_NODE_HXX
#define UHM_MESH_NODE_HXX

namespace uhm {
  typedef class Node_*    Node;
  typedef class Element_* Element;
  typedef class Mesh_*    Mesh;

  bool node_valid(Node n);

  // ----------------------------------------------------------------
  // ** Node class
  class Node_ : public Object_< std::pair <int,int> > {
  protected:
    // store nodal connectivity according to its generation
    std::set< Element > owner;

    // n_dof  - actual dof that the node has
    // p      - order from FE
    // kind   - boundary node and the others
    // offset - offset in assembled matrix, where this node places in mesh
    int n_dof, p, kind, offset;
    int marker;

    void _init(std::pair<int,int> id, int n_dof, int p,int kind);

  public:
    Node_();
    Node_(std::pair<int,int> id);
    Node_(std::pair<int,int> id, int n_dof, int p, int kind);
    virtual ~Node_();
    virtual bool disp();
    virtual bool disp(FILE *stream);

    bool operator<(const Node_ &b) const;

    void reset_owner();
   
    void set_kind   (int kind);
    void set_n_dof  (int n_dof);
    void set_p      (int p);
    void set_offset (int offset);
    void set_marker (int marker);
    
    int  get_kind();
    int  get_n_dof();
    int  get_p();
    int  get_n_owner();
    int  get_offset();
    int  get_marker();

    bool is_owned_by  (Element e);
    bool is_same_as   (Node n);

    void add_owner    (Element e);
    void remove_owner (Element e);

    void clean_connectivity();

    friend class Mesh_;
    friend bool orphan_graph(std::vector< Element > *orphan,
			     Element parent,
			     //                                
			     std::vector< int > *xadj,
			     std::vector< int > *adjncy,
			     std::vector< int > *adjwgt);
    friend bool build_tree_var_1(Mesh m);

    friend bool node_valid(Node n);
  };
  // ----------------------------------------------------------------
  // ** Definition
  inline void Node_::_init(std::pair<int,int> id, int n_dof, int p,int kind) {
    this->cookie = UHM_NODE_COOKIE;
    this->id = id;
    this->n_dof = n_dof;
    this->p = p;
    this->kind = kind;
    this->offset = 0;
    this->marker = -1;
  }
  inline bool Node_::operator<(const Node_ &b) const { 
    return (this->id < b.id); 
  }
  inline bool node_valid(Node n) {
    return (n && n->cookie == UHM_NODE_COOKIE); 
  }
}

#endif
