/*
  Copyright © 2011, Kyungjoo Kim
  All rights reserved.
  
  This file is part of LINAL.
  
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
#include "linal/common.hxx"
#include "linal/const.hxx"
#include "linal/util.hxx"
#include "linal/matrix.hxx"
#include "linal/flat.hxx"
#include "linal/hier.hxx"
#include "linal/operation.hxx"

namespace linal {
  void head_graphviz(FILE *fp, char *name) { 
    fprintf(fp, "digraph %s { \n", name);  
  }
  void tail_graphviz(FILE *fp) {
    fprintf(fp, "} \n");
  }
  
  int push_graphviz(Hier_ A, int begin, 
		    std::vector< std::pair<int,double> > &viz) {
    viz.clear();
    for (int k2=0;k2<A.get_n();++k2) {
      for (int k1=0;k1<A.get_m();++k1) {
	std::pair<int,double> in;
	in.first  = begin++;
	in.second = sqrt(A(k1,k2).m * A(k1,k2).n);
	viz.push_back(in);
      }
    }
    return begin;
  }

  void write_graphviz(FILE *fp, char *task, char *style, char *color,
		      std::vector< std::pair<int,double> > &viz_in,
		      std::vector< std::pair<int,double> > &viz_out) {
    std::vector< std::pair<int,double> >::iterator in, out;
    for (in=viz_in.begin();in<viz_in.end();++in) {
      fprintf(fp, "%d [label=\"%s\",style=%s,color=%s]; %d -> {",
	     in->first, task, style, color, in->first);

      for (out=viz_out.begin();out<viz_out.end();++out) {
	fprintf(fp,"%d;",out->first);

      }
      fprintf(fp,"};\n");
    }      
  }
}
