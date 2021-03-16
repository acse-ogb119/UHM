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
  static double inv_norm1_var1(int seed, int uplo, int diag, Flat_ A);
  static double inv_norm1_var2(int uplo, int diag, Flat_ A);

  double norm1(int uplo, int diag, Flat_ A) {
    double rval;

    Flat_ norm;
    norm.create(LINAL_DOUBLE_REAL, 1, 1);

    if (uplo && diag) {
      Flat_ A_trig;
      A_trig.create(A.get_data_type(), A.get_m(), A.get_n());
      FLA_Copy(~A, ~A_trig);
      FLA_Triangularize(uplo, diag, ~A_trig);
      FLA_Norm1(~A_trig, ~norm);
      A_trig.free();
    } else {
      FLA_Norm1(~A, ~norm);
    }
    rval = norm(0,0);
    norm.free();

    return rval;
  }

  double inv_norm1(int uplo, int diag, Flat_ A) {
    double m = (double)A.get_m();
    return inv_norm1(0.0, uplo, diag, A);
  }
  double inv_norm1(double sample, int uplo, int diag, Flat_ A) {
    assert(sample>=0.0 && sample<=1.0);
    int n = (int)(A.get_m()*sample);
    int s = A.get_m()/(n ? n : 1);   

    double val = 0.0;
    for (int i=-s;i<A.get_m();i+=s) {
      double tmp = inv_norm1_var1(i, uplo, diag, A);
      //      printf("seed %d, tmp %E, val %E\n", i, tmp, val);
      if (tmp > val) val = tmp;
    }
    
    return val;
  }

  static double inv_norm1_var1(int seed, int uplo, int diag, Flat_ A) {
    // square matrix
    assert(A.get_m() == A.get_n());
    double rval;
    int datatype = A.get_data_type();
    int m        = A.get_m();

    linal::Flat_ x, y, z, norm, rho;
    x.create     (datatype, m, 1);
    y.create     (datatype, m, 1);
    z.create     (datatype, m, 1);

    norm.create (LINAL_DOUBLE_REAL, 1, 1);
    rho.create  (LINAL_DOUBLE_REAL, 1, 1);

    // initial set up
    if (seed < 0) {
      for (int i=0;i<m;++i) 
        x(i,0) = 1.0/(double)m;
    } else {
      FLA_Obj_set_to_scalar(FLA_ZERO, ~x);
      x(seed,0) = 1.0;
    }

    // iteration until loop out condition satisfied
    for (int k=0;k<(m+1);++k) {
      FLA_Copy( ~x, ~y );
      FLA_Trsm( FLA_LEFT, uplo, FLA_NO_TRANSPOSE,
		diag, FLA_ONE, ~A, ~y );

      for (int i=0;i<m;++i)
	z(i,0) = y(i,0) >= 0.0 ? 1.0 : -1.0;

      FLA_Trsm( FLA_LEFT, uplo, FLA_TRANSPOSE,
		diag, FLA_ONE, ~A, ~z );

      FLA_Norm_inf( ~z, ~norm );
      FLA_Dot( ~z, ~x, ~rho );

      if ( (norm(0,0)-rho(0,0)) < LINAL_NORM_ERROR_TOL ) 
	break;

      // find the location of max value in z
      int index;      
      double val = 0.0;
      for (int i=0;i<m;++i) 
	if (val*val < z(i,0)*z(i,0)) {
	  index = i; 
	  val   = z(i,0);
	}
      
      // update x by setting x(index, 0) = 1.0;
      FLA_Obj_set_to_scalar(FLA_ZERO, ~x);
      x(index,0) = 1.0;
    }

    // norm
    FLA_Norm1( ~y, ~norm );
    rval = norm(0,0);

    rho.free();
    norm.free();

    z.free();
    y.free();
    x.free();

    return rval;
  }

  static double inv_norm1_var2(int uplo, int diag, Flat_ A) {
    // square matrix
    assert(A.get_m() == A.get_n());

    double rval;
    int datatype = A.get_data_type();
    int m        = A.get_m();

    linal::Flat_ x, y, z, r, norm, rho;
    x.create     (datatype, m, 1);
    y.create     (datatype, m, 1);
    z.create     (datatype, m, 1);

    norm.create (LINAL_DOUBLE_REAL, 1, 1);
    rho.create  (LINAL_DOUBLE_REAL, 1, 1);

    // ------------------------------------------------
    // initial set up 
    for (int i=0;i<m;++i) 
      y(i,0) = 1.0/(double)m;
    
    FLA_Trsm( FLA_LEFT, uplo, FLA_NO_TRANSPOSE,
	      diag, FLA_ONE, ~A, ~y );
    FLA_Norm1( ~y, ~norm );
    rval = norm(0,0);

    for (int i=0;i<m;++i)
      z(i,0) = y(i,0) >= 0.0 ? 1.0 : -1.0;

    FLA_Copy( ~z, ~x );
    FLA_Trsm( FLA_LEFT, uplo, FLA_TRANSPOSE,
	      diag, FLA_ONE, ~A, ~z );

    // ------------------------------------------------
    if (m>1) {

      // iteration until loop out condition satisfied
      for (int k=0;k<(m+1);++k) {

	// index = min( i : | z_index | == ||z||_inf )
	int    index = 0, index_back, flag1, flag2;      
	double val   = 0.0;
	for (int i=(m-1);i>=0;--i) 
	  if (val*val < z(i,0)*z(i,0)) {
	    index = i; 
	    val   = z(i,0);
	  }

	// update x by setting y(index, 0) = 1.0;
	FLA_Obj_set_to_scalar(FLA_ZERO, ~y);
	y(index,0) = 1.0;
		
	FLA_Trsm( FLA_LEFT, uplo, FLA_NO_TRANSPOSE,
		  diag, FLA_ONE, ~A, ~y );
	FLA_Norm1( ~y, ~norm );

	// check rval decreasing
	flag1 = norm(0,0) <= rval;
	rval  = norm(0,0);
	
	// check cyclic
	double sum = 0.0;
	for (int i=0;i<m;++i)
	  sum += abs( ( (y(i,0)>=0.0) ? 1.0 : -1.0 ) - x(i,0) );
	flag2 = sum < 1.0; 
	
	// modify vector
	if (flag1 && flag2) {
	  for (int i=0;i<m;++i)
	    z(i,0) = pow(-1.0, i)*(1.0 + (double)i/(m-1.0));

	  FLA_Trsm( FLA_LEFT, uplo, FLA_NO_TRANSPOSE,
		    diag, FLA_ONE, ~A, ~z );
	  FLA_Norm1( ~z, ~norm );

	  if ((2.0*norm(0,0) / (3.0*m)) > rval ) {
	    FLA_Copy( ~z, ~y );
	    rval = 2*norm(0,0) / (3.0*m);
	  }
	}
		
	// sign 
	for (int i=0;i<m;++i)
	  z(i,0) = y(i,0) >= 0.0 ? 1.0 : -1.0;
	
	FLA_Trsm( FLA_LEFT, uplo, FLA_TRANSPOSE,
		  diag, FLA_ONE, ~A, ~z );
	
	// breaking condition
	index_back = index;
        for (int i=(m-1);i>=0;--i)
          if (val*val < z(i,0)*z(i,0)) {
            index = i;
            val   = z(i,0);
          }
	
	if ( index == index_back )
	  break;
      }
    }

    rho.free();
    norm.free();

    z.free();
    y.free();
    x.free();

    return rval;
  }
}
