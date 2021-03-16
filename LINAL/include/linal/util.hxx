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
#ifndef LINAL_UTIL_HXX
#define LINAL_UTIL_HXX

namespace linal {
  // ----------------------------------------------------------------
  /*! set computing model
   */
  extern void set_computing_model(int model);
  extern int  get_computing_model();

  // ----------------------------------------------------------------
  /*! FLOP counting - From LAPACK working note #41
   */
  extern double get_flop_gemm       ( int is_complex, int m, int n, int k );
  extern double get_flop_syrk       ( int is_complex, int k, int n );
  extern double get_flop_trsm_lower ( int is_complex, int m, int n );
  extern double get_flop_trsm_upper ( int is_complex, int m, int n );
  extern double get_flop_trmm_lower ( int is_complex, int m, int n );
  extern double get_flop_trmm_upper ( int is_complex, int m, int n );
  extern double get_flop_lu         ( int is_complex, int m, int n );
  extern double get_flop_chol       ( int is_complex, int n );
  extern double get_flop_qr         ( int is_complex, int m, int n );
  extern double get_flop_q          ( int is_complex, int m, int n, int k );
  
  extern double get_memory          ( int is_complex, int m, int n );

  // ----------------------------------------------------------------
  /*! Basic File IO.
   */
  extern bool read_line(FILE *fp, char **line_out);
  extern bool open_file(char *fullpath, char *mode, FILE **fp);
  extern bool close_file(FILE *fp);

  // ----------------------------------------------------------------
  /*! Error handler
   */
#define LINAL_ERROR(eval,error_msg) linal::error((eval), __FILE__,__LINE__,(error_msg));  
  extern inline void error(const int eval,
                           const char *file, const int line, 
			   const char *error_msg) {
    // if eval fails, then thru error
    if (!eval) {
      std::ostringstream msg;
      msg << "[LINAL ERROR], " << file << " , " << line << std::endl
          << error_msg << std::endl;
      throw std::logic_error( msg.str() );
    }
  }

  // ----------------------------------------------------------------
  /*! Check util
   */
  extern inline bool exist(FLA_Obj &A) { return (A.base != NULL); }
  extern inline bool exist(FLA_Obj &A, FLA_Obj &B) { 
    return (A.base != NULL && B.base != NULL);
  }
  extern inline bool exist(FLA_Obj &A, FLA_Obj &B, FLA_Obj &C) {
    return (A.base != NULL && B.base != NULL && C.base != NULL);
  }

  extern inline bool check_all_scalar_type(int type) {
    return ( type == LINAL_INT ||
             type == LINAL_SINGLE_REAL ||
             type == LINAL_DOUBLE_REAL ||
             type == LINAL_SINGLE_COMPLEX ||
             type == LINAL_DOUBLE_COMPLEX );
  }

  extern inline bool check_double_scalr_type(int type) {
    return ( type == LINAL_INT ||
             type == LINAL_DOUBLE_REAL ||
             type == LINAL_DOUBLE_COMPLEX );
  }

  extern inline bool check_single_scalar_type(int type) {
    return ( type == LINAL_INT ||
             type == LINAL_SINGLE_REAL ||
             type == LINAL_SINGLE_COMPLEX );
  }

}

#endif
