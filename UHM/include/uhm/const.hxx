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
#ifndef UHM_CONST_HXX
#define UHM_CONST_HXX

// ** cookies
#define UHM_OBJECT_COOKIE          10
#define UHM_NODE_COOKIE           100
#define UHM_ELEMENT_COOKIE        110
#define UHM_MESH_COOKIE           120
#define UHM_SCHEDULER_COOKIE      200
#define UHM_HELPER_COOKIE         300
#define UHM_MATRIX_FLA_COOKIE    1000
#define UHM_MATRIX_EL_COOKIE     2000


// control variable
#define UHM_DISSECTION_TASK_SIZE 2000
#define UHM_UNROLL_N                8

// should be re-defined 
#define UHM_INT            LINAL_INT

#define UHM_REAL           LINAL_REAL
#define UHM_SINGLE_REAL    LINAL_SINGLE_REAL
#define UHM_DOUBLE_REAL    LINAL_DOUBLE_REAL

#define UHM_COMPLEX        LINAL_COMPLEX
#define UHM_SINGLE_COMPLEX LINAL_SINGLE_COMPLEX
#define UHM_DOUBLE_COMPLEX LINAL_DOUBLE_COMPLEX


#define UHM_FULL_MATRIX       100
#define UHM_LOWER_TRIANGULAR  200
#define UHM_UPPER_TRIANGULAR  300

#define UHM_LEAF_TO_ROOT      1
#define UHM_ROOT_TO_LEAF      0

enum { UHM_PHYSICS_SINGLE=1,     UHM_PHYSICS_MULTI };
enum { UHM_NODE_KIND_DEFAULT=0, UHM_NODE_KIND_BOUNDARY };
enum { UHM_DISP_ALL=1, UHM_DISP_NODE, UHM_DISP_ELEMENT, UHM_DISP_MATRIX };
enum { UHM_LHS=1, UHM_RHS };
enum { UHM_ATL=1, UHM_ATR, UHM_ABL, UHM_ABR, UHM_P, UHM_T,
       UHM_XT, UHM_XB, UHM_BT, UHM_BB, UHM_RT, UHM_RB, UHM_END };
enum { UHM_CHOL=1,   UHM_LU_NOPIV,     UHM_LU_PIV,     UHM_LU_INCPIV,     
       UHM_QR,       UHM_QRLQ,         UHM_RRQRLQ,     UHM_SVD,
       UHM_CHOL_OOC, UHM_LU_NOPIV_OOC, UHM_LU_PIV_OOC, UHM_LU_INCPIV_OOC, 
       UHM_QR_OOC,   UHM_QRLQ_OOC,     UHM_RRQRLQ_OOC, UHM_SVD_OOC};
enum { UHM_NOT_SEPARATED=0, UHM_SEPARATED_FACTOR, UHM_SEPARATED_SCHUR };
enum { UHM_UNASSEMBLED=0, UHM_ASSEMBLED };

// ** uhm name space
namespace uhm {
  typedef class Element_*   Element;
  typedef class Node_*      Node;
  typedef class Mesh_*      Mesh;
  typedef class Matrix_*    Matrix;
  typedef class Scheduler_* Scheduler;

  //  const void*     nil_buffer     = NULL;
  const Element   nil_element    = NULL;
  const Node      nil_node       = NULL;
  const Mesh      nil_mesh       = NULL;
  const Matrix    nil_matrix     = NULL;
  const Scheduler nil_scheduler  = NULL;
}


#endif
