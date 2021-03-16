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

#include "uhm/object.hxx"

#include "uhm/operation/scheduler.hxx"
#include "uhm/operation/element.hxx"

#include "uhm/matrix/uhm/matrix.hxx"

#include "uhm/mesh/mesh.hxx"


namespace uhm {
  // --------------------------------------------------------------
  // ** LU nopiv
  void Mesh_::chol_with_free() {
    assert(this->get_scheduler()->is_loaded());
    Scheduler s = this->get_scheduler();
#ifdef UHM_MULTITHREADING_ENABLE
    s->execute_tree(&op_chol_with_merge_and_free, true);
#else
    s->execute_elements_seq(&op_chol_with_merge_and_free, true);
#endif
  }

  void Mesh_::chol_without_free() {
    assert(this->get_scheduler()->is_loaded());
    Scheduler s = this->get_scheduler();
#ifdef UHM_MULTITHREADING_ENABLE
    s->execute_tree(&op_chol_with_merge_and_no_free, true);
#else
    s->execute_elements_seq(&op_chol_with_merge_and_no_free, true);
#endif
  }

  void Mesh_::solve_chol_1() {
    assert(this->get_scheduler()->is_loaded());
    Scheduler s = this->get_scheduler();
#ifdef UHM_MULTITHREADING_ENABLE
    s->execute_tree(&op_solve_chol_1_x_with_merge, 
		    true);
#else
    s->execute_elements_seq(&op_solve_chol_1_x_with_merge, 
			    true);
#endif
  }



  void Mesh_::solve_chol_2() {
    assert(this->get_scheduler()->is_loaded());
    Scheduler s = this->get_scheduler();
#ifdef UHM_MULTITHREADING_ENABLE
    s->execute_tree(&op_solve_chol_2_x_with_branch, 
		    false);
#else
    s->execute_elements_seq(&op_solve_chol_2_x_with_branch, 
			    false);
#endif
  }

  void Mesh_::solve_chol() {
    this->solve_chol_1();
    this->solve_chol_2();
  }


  void Mesh_::check_chol_1() {
    assert(this->get_scheduler()->is_loaded());
    Scheduler s = this->get_scheduler();
#ifdef UHM_MULTITHREADING_ENABLE
    s->execute_tree(&op_check_chol_1, false);
#else
    s->execute_elements_seq(&op_check_chol_1, false);
#endif
  }

  void Mesh_::check_chol_2() {
    assert(this->get_scheduler()->is_loaded());
    Scheduler s = this->get_scheduler();
#ifdef UHM_MULTITHREADING_ENABLE
    s->execute_tree(&op_check_chol_2, true);
#else
    s->execute_elements_seq(&op_check_chol_2, true);
#endif
  }


  void Mesh_::check_chol() {
    assert(this->get_scheduler()->is_loaded());
    Scheduler s = this->get_scheduler();

    this->check_chol_1();
    this->check_chol_2();

#ifdef UHM_MULTITHREADING_ENABLE
    s->execute_elements_seq(&op_check_solution, true);
#else
    s->execute_elements_seq(&op_check_solution, true);
#endif
  }
  
  void Mesh_::improve_chol() {
    // ** Not working... i don't know why...it is supposed to be working
    assert(this->get_scheduler()->is_loaded());
    Scheduler s = this->get_scheduler();
#ifdef UHM_MULTITHREADING_ENABLE

#else
    s->execute_elements_seq(&op_solve_chol_1_r_with_merge, 
			    true);
    s->execute_elements_seq(&op_solve_chol_2_r_with_branch, 
			    false);
    s->execute_elements_seq(&op_improve_solution, true);
#endif
  }
}
