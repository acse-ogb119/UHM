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
#ifndef UHM_OPERATION_MESH_HXX
#define UHM_OPERATION_MESH_HXX

namespace uhm {
  // Var 1 :: Heavy edge matching
  // Var 2 :: Metis graph partitioning
  // Var 3 :: Metis nested dissection
  extern bool orphan_graph(std::vector< Element > *orphan,
			   Element parent,
			   //                         
			   std::vector< int > *xadj,
			   std::vector< int > *adjncy,
			   std::vector< int > *adjwgt);

  // ** front end
  extern bool build_tree(Mesh m);

  // ** various implementations
  extern bool build_tree_var_1(Mesh m);

  extern bool build_tree_var_2(Mesh m);
  extern bool build_tree_var_2(Mesh m, int nparts);

  extern bool build_tree_var_3(Mesh m);

  extern bool build_tree_var_4(Mesh m);

  extern bool build_tree_var_5(Mesh m);
}


#endif
