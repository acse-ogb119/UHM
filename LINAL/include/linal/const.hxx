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
#ifndef LINAL_CONST_HXX_
#define LINAL_CONST_HXX_

/*!
  Define the constants and null objects.
*/
#define LINAL_INT               FLA_INT

#define LINAL_REAL              FLA_DOUBLE
#define LINAL_SINGLE_REAL       FLA_FLOAT
#define LINAL_DOUBLE_REAL       FLA_DOUBLE
 
#define LINAL_COMPLEX           FLA_DOUBLE_COMPLEX
#define LINAL_SINGLE_COMPLEX    FLA_COMPLEX
#define LINAL_DOUBLE_COMPLEX    FLA_DOUBLE_COMPLEX

#define LINAL_MATRIX            FLA_MATRIX

#define LINAL_TRANSPOSE         FLA_TRANSPOSE
#define LINAL_NO_TRANSPOSE      FLA_NO_TRANSPOSE

#define LINAL_CONJ_TRANSPOSE    FLA_CONJ_TRANSPOSE
#define LINAL_CONJ_NO_TRANSPOSE FLA_CONJ_NO_TRANSPOSE

#define LINAL_UPPER_TRIANGULAR  FLA_UPPER_TRIANGULAR
#define LINAL_LOWER_TRIANGULAR  FLA_LOWER_TRIANGULAR

#define LINAL_NORM_ERROR_TOL    1.0e-12

#define LINAL_CPU               1
#define LINAL_GPU               100
#define LINAL_CPU_GPU           200


// Null matrix
namespace linal {
  typedef class Matrix_* Matrix;
  
  const Matrix nil_matrix = NULL;
}

#endif



// FLAME constants
// 00072 // FLA_Side
// 00073 #define FLA_TOP               200
// 00074 #define FLA_BOTTOM            201
// 00075 #define FLA_LEFT              210
// 00076 #define FLA_RIGHT             211
// 00077 #define FLA_SIDE_MASK         0x1
// 00078 
// 00079 // FLA_Uplo
// 00080 #define FLA_LOWER_TRIANGULAR  300
// 00081 #define FLA_UPPER_TRIANGULAR  301
// 00082 #define FLA_UPLO_MASK         0x1
// 00083 
// 00084 // FLA_Trans
// 00085 #define FLA_NO_TRANSPOSE      400
// 00086 #define FLA_TRANSPOSE         401
// 00087 #define FLA_CONJ_TRANSPOSE    402
// 00088 #define FLA_CONJ_NO_TRANSPOSE 403
// 00089 #define FLA_TRANS_MASK        0x3
// 00090 
// 00091 // FLA_Conj
// 00092 #define FLA_NO_CONJUGATE      450
// 00093 #define FLA_CONJUGATE         451
