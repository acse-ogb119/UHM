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

#include "uhm/mesh/node.hxx"
#include "uhm/mesh/element.hxx"

#include "uhm/matrix/uhm/matrix.hxx"

#include "uhm/mesh/mesh.hxx"



namespace uhm {
  // --------------------------------------------------------------
  // ** LU nopiv
  void Mesh_::export_graphviz_lu_nopiv(char *full_path, int bmn) {

    FILE *fp;
    assert(open_file(full_path, "w", &fp));
    linal::head_graphviz(fp, "Mesh_tasks");

    int start=0, end=0;
    std::map< int, Element_ >::iterator it;
    for (it=this->elements.begin();it!=this->elements.end();++it) {
      it->second.write_lu_nopiv(fp, bmn, start, end);
      start = end;
    }

    for (it=this->elements.begin();it!=this->elements.end();++it) {
      Element e = &(it->second);
      if (!e->is_orphan())
	fprintf(fp, "%d -> {%d;};\n", 
		e->get_marker(1),e->get_parent()->get_marker(0));
    }

    linal::tail_graphviz(fp);
    assert(close_file(fp));
  }
  
  void Element_::write_lu_nopiv(FILE *fp, int bmn, int &start, int &end) {
    
    // write on the file
    linal::Flat_ fA;
    linal::Hier_ hA;

    std::pair<int,int> n_dof = this->get_n_dof();
    int m   = n_dof.first + n_dof.second;
    fA.create(LINAL_REAL, m, m);
    hA.create(fA, bmn, bmn);
    
    linal::lu_nopiv(fp, start, end, hA);

    printf("start %d, end %d\n", start, end-1);
    this->set_marker(0,start);
    this->set_marker(1,end-1);

    fA.free();
    hA.free();
  }
}
