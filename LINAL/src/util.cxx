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

namespace linal {
  static int computing_model = 0;



  // --------------------------------------------------------------
  // ** Computing model setup
  void set_computing_model(int model) {
    computing_model = model;
  }
  int  get_computing_model() {
    return computing_model;
  }

  // --------------------------------------------------------------
  // ** FLOP and Memory counting 
  /*!
    Flop counting is from LAWN 41
  */
#define LINAL_MUL(is_complex) ( (is_complex) ?  (6.0) : (1.0) )
#define LINAL_ADD(is_complex) ( (is_complex) ?  (2.0) : (1.0) )
#define LINAL_MEM(is_complex) ( (is_complex) ? (16.0) : (8.0) )

  double get_flop_gemm( int is_complex, int mm, int nn, int kk ) {
    double m = (double)mm;    double n = (double)nn;    double k = (double)kk;
    return ( LINAL_MUL(is_complex)*(m*n*k) + 
             LINAL_ADD(is_complex)*(m*n*k) );
  }
  double get_flop_syrk( int is_complex, int kk, int nn ) {
    double k = (double)kk;    double n = (double)nn;    
    return ( LINAL_MUL(is_complex)*(0.5*k*n*(n+1.0)) + 
             LINAL_ADD(is_complex)*(0.5*k*n*(n+1.0)) );
  }
  double get_flop_trsm_lower( int is_complex, int mm, int nn ) {
    double m = (double)mm;    double n = (double)nn;    
    return ( LINAL_MUL(is_complex)*(0.5*n*m*(m+1.0)) + 
             LINAL_ADD(is_complex)*(0.5*n*m*(m-1.0)) );
  }
  double get_flop_trsm_upper( int is_complex, int mm, int nn ) {
    double m = (double)mm;    double n = (double)nn;    
    return ( LINAL_MUL(is_complex)*(0.5*m*n*(n+1.0)) + 
             LINAL_ADD(is_complex)*(0.5*m*n*(n-1.0)) );
  }
  double get_flop_trmm_lower( int is_complex, int mm, int nn ) {
    double m = (double)mm;    double n = (double)nn;    
    return get_flop_trsm_lower( is_complex, m, n);
  }
  double get_flop_trmm_upper( int is_complex, int mm, int nn ) {
    double m = (double)mm;    double n = (double)nn;    
    return get_flop_trsm_upper( is_complex, m, n);
  }
  double get_flop_lu( int is_complex, int mm, int nn ) {
    double m = (double)mm;    double n = (double)nn;    
    if (m > n) 
      return ( LINAL_MUL(is_complex)*(0.5*m*n*n-(1.0/6.0)*n*n*n+0.5*m*n-0.5*n*n+(2.0/3.0)*n) +
               LINAL_ADD(is_complex)*(0.5*m*n*n-(1.0/6.0)*n*n*n-0.5*m*n+        (1.0/6.0)*n) );
    else
      return ( LINAL_MUL(is_complex)*(0.5*n*m*m-(1.0/6.0)*m*m*m+0.5*n*m-0.5*m*m+(2.0/3.0)*m) +
               LINAL_ADD(is_complex)*(0.5*n*m*m-(1.0/6.0)*m*m*m-0.5*n*m+        (1.0/6.0)*m) );
  }
  double get_flop_chol( int is_complex, int nn ) {
    double n = (double)nn;
    return ( LINAL_MUL(is_complex)*((1.0/6.0)*n*n*n+0.5*n*n+(1.0/3.0)*n) +
             LINAL_ADD(is_complex)*((1.0/6.0)*n*n*n-        (1.0/6.0)*n) );

  }
  double get_flop_qr( int is_complex, int mm, int nn ) {
    double m = (double)mm;    double n = (double)nn;    
    if (m > n)
      return ( LINAL_MUL(is_complex)*(m*n*n-(1.0/3.0)*n*n*n+m*n+0.5*n*n+(23.0/6.0)*n) +
               LINAL_ADD(is_complex)*(m*n*n-(1.0/3.0)*n*n*n+   +0.5*n*n+( 5.0/6.0)*n) );
    else
      return ( LINAL_MUL(is_complex)*(n*m*m-(1.0/3.0)*m*m*m+2.0*n*m-0.5*m*m+(23.0/6.0)*m) +
               LINAL_ADD(is_complex)*(n*m*m-(1.0/3.0)*m*m*m+    n*m-0.5*m*m+( 5.0/6.0)*m) );
  }
  double get_flop_q( int is_complex, int mm, int nn, int kk ) {
    double m = (double)mm;    double n = (double)nn;    double k = (double)kk;
    return ( LINAL_MUL(is_complex)*(2.0*m*n*k-(m+n)*k*k+(2.0/3.0)*k*k*k+2.0*n*k-k*k-(5.0/3.0)*k) +
             LINAL_ADD(is_complex)*(2.0*m*n*k-(m+n)*k*k+(2.0/3.0)*k*k*k+    n*k-m*k+(1.0/3.0)*k) );
  }

  double get_memory( int is_complex, int mm, int nn ) {
    double m = (double)mm;    double n = (double)nn;
    return LINAL_MEM(is_complex)*(m*n);
  }
  
  // --------------------------------------------------------------
  // ** Basic File IO
  /*! 
    Read one line from file pointer. Max number of characters in
    one line is 10000.
  */
  bool read_line(FILE *fp, char **line_out) {
    char line_in[10000];
    char *c;

    while (fgets(line_in, 1000, fp) != NULL) {
      c = line_in;
      while (*c && *c <= ' ') c++;

      // skip the comment line
      if (*c && *c != '#') {
        *line_out = line_in;
        return true;
      }
    }
    *line_out = NULL;
    return false;
  }

  /*!
    Open file. If it succeed, it return "true" value
  */
  bool open_file(char *fullpath, char *mode, FILE **fp) {
    *fp = fopen(fullpath, mode);
    if (*fp == NULL) {
      fprintf(stderr, "fail to open file : %s\n", fullpath);
      return false;
    }
    return true;
  }

  /*!
    Close file. If it succeed, it return "true" value
  */
  bool close_file(FILE *fp) {
    if (fclose(fp)) {
      fprintf(stderr,"fail to close file");
      return false;
    }
    return true;
  }
}
