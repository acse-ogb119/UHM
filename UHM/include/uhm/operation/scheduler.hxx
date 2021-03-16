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
#ifndef UHM_SCHEDULER_HXX
#define UHM_SCHEDULER_HXX

namespace uhm {
  typedef class Node_*      Node;
  typedef class Element_*   Element;
  typedef class Mesh_*      Mesh;
  typedef class Scheduler_* Scheduler;

  bool scheduler_valid(Scheduler s);

  // ----------------------------------------------------------------
  // ** Scheduler class
  class Scheduler_: public Object_<int> {
  protected:
    std::map    < int, std::vector< Element > > elements;
    std::vector < Element >                     leaves;

    void _init(int);
  public:
    Scheduler_();
    Scheduler_(int id);
    virtual ~Scheduler_();
    virtual bool disp();
    virtual bool disp(int datasize, int flop);
    virtual bool disp(FILE *stream, int datasize, int flop);

    void load(Mesh m);
    void unload();
        
    int  get_n_elements();
    int  get_n_leaves();

    int  is_loaded();

    void get_orphan(std::vector<Element>& orphan);

    bool execute_tree(bool (*op_func)(Element), int is_leaf2root);
    bool execute_elements_seq(bool (*op_func)(Element), int is_leaf2root);
    bool execute_elements_par(bool (*op_func)(Element), int is_leaf2root);
    bool execute_leaves_seq(bool (*op_func)(Element));
    bool execute_leaves_par(bool (*op_func)(Element));


    // following is not implemented yet : actually it is not necessary
    bool execute_elements(bool (*op_func_1)(Element), 
			  bool (*op_func_2)(Element),
			  int is_leaf2root);
  
    // friends
    friend bool scheduler_valid(Scheduler s);
  };
  // ----------------------------------------------------------------
  // ** Definition
  inline void Scheduler_::_init(int) {
    this->cookie = UHM_SCHEDULER_COOKIE;
    this->id = id;
  }
  inline bool scheduler_valid(Scheduler s) { 
    return (s && s->cookie == UHM_SCHEDULER_COOKIE);
  }
}


#endif
