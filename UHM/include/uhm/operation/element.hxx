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
#ifndef UHM_OPERATION_ELEMENT_HXX
#define UHM_OPERATION_ELEMENT_HXX


namespace uhm {
  // ---------------------------------------------------
  extern bool op_create_matrix_buffer_without_schur(Element e);
  extern bool op_create_matrix_buffer_with_schur   (Element e);

  extern bool op_restore_connectivity              (Element e);
  extern bool op_numbering                         (Element e);
  extern bool op_update_connectivity               (Element e);
  extern bool op_arrange_nodes                     (Element e);
  extern bool op_update_generation                 (Element e);
  // ---------------------------------------------------
  extern bool op_merge_full_with_free              (Element e);
  extern bool op_merge_full_without_free           (Element e);
  extern bool op_merge_full_with_ooc               (Element e);

  extern bool op_merge_rhs_x                       (Element e);
  extern bool op_merge_rhs_r                       (Element e);

  extern bool op_branch_rhs_x                      (Element e);
  extern bool op_branch_rhs_r                      (Element e);
  // ---------------------------------------------------
  extern bool op_chol_without_merge                (Element e);
  extern bool op_chol_with_merge_and_free          (Element e);
  extern bool op_chol_with_merge_and_no_free       (Element e);

  extern bool op_solve_chol_1_x_without_merge      (Element e);
  extern bool op_solve_chol_1_x_with_merge         (Element e);
  extern bool op_solve_chol_2_x_without_branch     (Element e);
  extern bool op_solve_chol_2_x_with_branch        (Element e);

  extern bool op_solve_chol_1_r_without_merge      (Element e);
  extern bool op_solve_chol_1_r_with_merge         (Element e);
  extern bool op_solve_chol_2_r_without_branch     (Element e);
  extern bool op_solve_chol_2_r_with_branch        (Element e);

  extern bool op_check_chol_1                      (Element e);
  extern bool op_check_chol_2                      (Element e);
  // ---------------------------------------------------
  extern bool op_lu_nopiv_without_merge            (Element e);
  extern bool op_lu_nopiv_with_merge_and_free      (Element e);
  extern bool op_lu_nopiv_with_merge_and_no_free   (Element e);

  extern bool op_solve_lu_nopiv_1_x_without_merge  (Element e);
  extern bool op_solve_lu_nopiv_1_x_with_merge     (Element e);
  extern bool op_solve_lu_nopiv_2_x_without_branch (Element e);
  extern bool op_solve_lu_nopiv_2_x_with_branch    (Element e);

  extern bool op_solve_lu_nopiv_1_r_without_merge  (Element e);
  extern bool op_solve_lu_nopiv_1_r_with_merge     (Element e);
  extern bool op_solve_lu_nopiv_2_r_without_branch (Element e);
  extern bool op_solve_lu_nopiv_2_r_with_branch    (Element e);

  extern bool op_check_lu_nopiv_1                  (Element e);
  extern bool op_check_lu_nopiv_2                  (Element e);
  // ---------------------------------------------------
  extern bool op_lu_incpiv_without_merge           (Element e);
  extern bool op_lu_incpiv_with_merge_and_free     (Element e);
  extern bool op_lu_incpiv_with_merge_and_no_free  (Element e);
  // ---------------------------------------------------
  extern bool op_lu_piv_without_merge              (Element e);
  extern bool op_lu_piv_with_merge_and_free        (Element e);
  extern bool op_lu_piv_with_merge_and_no_free     (Element e);


  extern bool op_solve_lu_piv_1_x_without_merge    (Element e);
  extern bool op_solve_lu_piv_1_x_with_merge       (Element e);
  extern bool op_solve_lu_piv_2_x_without_branch   (Element e);
  extern bool op_solve_lu_piv_2_x_with_branch      (Element e);

  extern bool op_solve_lu_piv_1_r_without_merge    (Element e);
  extern bool op_solve_lu_piv_1_r_with_merge       (Element e);
  extern bool op_solve_lu_piv_2_r_without_branch   (Element e);
  extern bool op_solve_lu_piv_2_r_with_branch      (Element e);

  extern bool op_check_lu_piv_1                    (Element e);
  extern bool op_check_lu_piv_2                    (Element e);
  // ---------------------------------------------------

  extern bool op_qr_without_merge                  (Element e);
  extern bool op_qr_with_merge_and_free            (Element e);
  extern bool op_qr_with_merge_and_no_free         (Element e);
  // ---------------------------------------------------
  extern bool op_solve_qr_1_x_without_merge        (Element e);
  extern bool op_solve_qr_1_x_with_merge           (Element e);
  extern bool op_solve_qr_2_x_without_branch       (Element e);
  extern bool op_solve_qr_2_x_with_branch          (Element e);

  extern bool op_solve_qr_1_r_without_merge        (Element e);
  extern bool op_solve_qr_1_r_with_merge           (Element e);
  extern bool op_solve_qr_2_r_without_branch       (Element e);
  extern bool op_solve_qr_2_r_with_branch          (Element e);

  extern bool op_check_qr_1                        (Element e);
  extern bool op_check_qr_2                        (Element e);
  // ---------------------------------------------------

  // ---------------------------------------------------

  extern bool op_chol_with_merge_and_ooc           (Element e);
  extern bool op_lu_nopiv_with_merge_and_ooc       (Element e);
  extern bool op_lu_piv_with_merge_and_ooc         (Element e);
  extern bool op_lu_incpiv_with_merge_and_ooc      (Element e);
  extern bool op_qr_with_merge_and_ooc             (Element e);

  extern bool op_solve_lu_piv_1_x_without_merge_ooc(Element e);
  extern bool op_solve_lu_piv_1_x_with_merge_ooc   (Element e);
  extern bool op_solve_lu_piv_2_x_without_branch_ooc(Element e);
  extern bool op_solve_lu_piv_2_x_with_branch_ooc  (Element e);

  extern bool op_solve_lu_piv_1_r_without_merge_ooc(Element e);
  extern bool op_solve_lu_piv_1_r_with_merge_ooc   (Element e);
  extern bool op_solve_lu_piv_2_r_without_branch_ooc(Element e);
  extern bool op_solve_lu_piv_2_r_with_branch_ooc  (Element e);

  extern bool op_check_lu_piv_1_ooc                (Element e);
  extern bool op_check_lu_piv_2_ooc                (Element e);
  // ---------------------------------------------------
  extern bool op_check_solution                    (Element e);
  extern bool op_improve_solution                  (Element e);
  // ---------------------------------------------------
  // ---------------------------------------------------
  extern bool op_test                              (Element e);
}


#endif
