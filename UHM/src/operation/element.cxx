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
#include "uhm/operation/mesh.hxx"
#include "uhm/operation/element.hxx"

#include "uhm/mesh/node.hxx"
#include "uhm/mesh/element.hxx"

#include "uhm/matrix/uhm/matrix.hxx"

#include "uhm/mesh/mesh.hxx"

#include "uhm/matrix/uhm/helper.hxx"

namespace uhm {
  // --------------------------------------------------------------
  static bool op_create_matrix_buffer( Element e, int is_schur );

  static bool op_write_ooc ( Element e );
  static bool op_read_ooc  ( Element e );
  static bool op_merge     ( Element e, int a, int x, int b, int r, 
			     int is_create_buffer, int is_buffer_free );
  static bool op_branch    ( Element e, int x, int b, int r );
  static bool op_decompose ( Element e, int type, int is_merge, int free_option );
  static bool op_solve     ( Element e, int type, int level, int x, int r,
			     int is_merge, int is_branch );
  static bool op_check     ( Element e, int type, int level );

  // --------------------------------------------------------------
  static bool op_create_matrix_buffer( Element e, int is_schur ) {
    assert(element_valid(e));
    assert(e->is_matrix_created());
    
    // default is ABR which is schur complement is not allocated
    for (int i=UHM_ATL;i<UHM_END;i++)
      if (i!=UHM_ABR)
        e->get_matrix()->create_buffer(i);
    
    // when it has been requested, then allocte for schur
    if (is_schur)
      e->get_matrix()->create_buffer(UHM_ABR);

    return true;
  }

  static bool op_create_buffer_ooc(Element e) {
    assert(element_valid(e));
    for (int i=UHM_ATL;i<UHM_P;++i) 
      e->get_matrix()->create_buffer(i);
    return true;
  }

  static bool op_free_buffer_ooc(Element e) {
    assert(element_valid(e));
    for (int i=UHM_ATL;i<UHM_P;++i) 
      e->get_matrix()->free_buffer(i);
    return true;
  }

  static bool op_write_ooc(Element e) {
    assert(element_valid(e));
    for (int i=UHM_ATL;i<UHM_P;++i) {
      FILE *fp;
      char fullpath[256], tmp[64];

      // name the file
      strcpy(fullpath, get_ooc_dir());
      sprintf(tmp, "/_uhm_%d_mat_%d_", e->get_id(), i);
      strncat(fullpath, tmp, strlen(tmp));

      // write on the file
      assert(open_file(fullpath, "wb", &fp));
      e->get_matrix()->write_to_ooc(fp, i);
      assert(close_file(fp));
    }
    return true;
  }

  static bool op_read_ooc(Element e) {
    assert(element_valid(e));
    for (int i=UHM_ATL;i<UHM_P;++i) {
      FILE *fp;
      char fullpath[256], tmp[64];

      // name the file
      strcpy(fullpath, get_ooc_dir());
      sprintf(tmp, "/_uhm_%d_mat_%d_", e->get_id(), i);
      strncat(fullpath, tmp, strlen(tmp));

      // write on the file
      assert(open_file(fullpath, "rb", &fp));
      e->get_matrix()->read_from_ooc(fp, i);
      assert(close_file(fp));
    }
    return true;
  }

  static bool op_merge(Element e, int a, int x, int b, int r, 
		       int is_create_buffer, 
		       int free_option) {
    assert(element_valid(e) && e->is_matrix_created());
    { // create_matrix buffer
      if (is_create_buffer) {
	if (a) 
	  for (int i=UHM_ATL;i<UHM_XT;++i) 
	    e->get_matrix()->create_buffer(i);
	if (b || x || r) 
	  for (int i=UHM_XT;i<UHM_END;++i) 
	    e->get_matrix()->create_buffer(i);
      }
    }
    
    { // collect children's matrices
      for (int i=0;i<e->get_n_children();++i) {
	Element c = e->get_child(i);
	assert(element_valid(c));

	Helper_ h(e, c);
	h.set_mapper();

	if (a) h.merge_A();
	if (b) h.merge_rhs_b();
	if (x) h.merge_rhs_x();
	if (r) h.merge_rhs_r();
      }

      switch (free_option) {
      case 0: {	break; }
      case 1: { // without ooc
	for (int i=0;i<e->get_n_children();++i) {
	  Element c = e->get_child(i);
	  assert(element_valid(c));
	  c->get_matrix()->free_buffer(UHM_ABR);
	}
	break; 
      }
      case 2: { // with ooc
	for (int i=0;i<e->get_n_children();++i) {
	  Element c = e->get_child(i);
	  assert(op_write_ooc(c));
	  assert(op_free_buffer_ooc(c));
	}
	break; 
      }
      }	  
    }
    return true;
  }

  static bool op_branch(Element e, int x, int b, int r) {
    assert(element_valid(e) && e->is_matrix_created());
    for (int i=0;i<e->get_n_children();++i) {
      Element c = e->get_child(i);
      assert(element_valid(c));
      
      Helper_ h(e, c);
      h.set_mapper();
      if (x) h.branch_rhs_x();
      if (b) h.branch_rhs_b();
      if (r) h.branch_rhs_r();
    }
    return true;
  }    

  static bool op_decompose(Element e, int type, int is_merge, int free_option) {
    assert(element_valid(e) && e->is_matrix_created());
    if (e->is_matrix_reusable()) {
      // do nothing
    } else {
      if (is_merge) {
	switch (free_option) {
	case 0:
	  op_merge_full_without_free(e);
	  break;
	case 1:
	  op_merge_full_with_free(e);
	  break;
	case 2:
	  op_merge_full_with_ooc(e);
	  break;
	}
      }
      switch (type) {
      case UHM_CHOL     : e->get_matrix()->chol();      break;
      case UHM_LU_NOPIV : e->get_matrix()->lu_nopiv();  break;
      case UHM_LU_INCPIV: e->get_matrix()->lu_incpiv(); break;
      case UHM_LU_PIV   : e->get_matrix()->lu_piv();    break;
      case UHM_QR       : e->get_matrix()->qr();        break;
      }
      if (e->is_orphan()) {
	switch(free_option) {
	case 0:
	case 1:
	  break;
	case 2:
          // QR is not considered for OOC yet.
	  assert(op_write_ooc(e));
	  assert(op_free_buffer_ooc(e));
	  break;
	}
      }
    }
    return true;
  }
  
  static bool op_solve(Element e, int type, int level, 
		       int x, int r,
		       int is_merge, int is_branch) {
    assert(element_valid(e) && e->is_matrix_created() &&
	   !(is_merge && is_branch) && !(x && r));

    if (is_merge) {
      op_merge(e, 0, x, 0, r, 0, 0);
    }
    if (level == 1) {
      if (x) {
	switch (type) {
	case UHM_CHOL       : e->get_matrix()->solve_chol_1_x();      break;
	case UHM_LU_NOPIV   : e->get_matrix()->solve_lu_nopiv_1_x();  break;
	case UHM_LU_PIV     : e->get_matrix()->solve_lu_piv_1_x();    break;
        case UHM_QR         : e->get_matrix()->solve_qr_1_x();        break;
	}
      }
      if (r) {
	switch (type) {
	case UHM_CHOL     : e->get_matrix()->solve_chol_1_r();      break;
	case UHM_LU_NOPIV : e->get_matrix()->solve_lu_nopiv_1_r();  break;
	case UHM_LU_PIV   : e->get_matrix()->solve_lu_piv_1_r();    break;
        case UHM_QR       : e->get_matrix()->solve_qr_1_r();        break;
	}
      }
    } 
    if (level == 2) {
      if (x) {
	switch (type) {
	case UHM_CHOL     : e->get_matrix()->solve_chol_2_x();      break;
	case UHM_LU_NOPIV : e->get_matrix()->solve_lu_nopiv_2_x();  break;
	case UHM_LU_PIV   : e->get_matrix()->solve_lu_piv_2_x();    break;
        case UHM_QR       : e->get_matrix()->solve_qr_2_x();        break;
	}
      }
      if (r) {
	switch (type) {
	case UHM_CHOL     : e->get_matrix()->solve_chol_2_r();      break;
	case UHM_LU_NOPIV : e->get_matrix()->solve_lu_nopiv_2_r();  break;
	case UHM_LU_PIV   : e->get_matrix()->solve_lu_piv_2_r();    break;
        case UHM_QR       : e->get_matrix()->solve_qr_2_r();        break;
	}
      }
    } 

    if (is_branch) {
      op_branch(e, x, 0, r);
    }
    
    return true;
  }

  static bool op_check(Element e, int type, int level) {
    // Todo :: for the test purpose check routine is not parallelized
    // I want coorect residual right now
    assert(element_valid(e) && e->is_matrix_created());

    if (level == 1) {
      switch (type) {
      case UHM_CHOL     : e->get_matrix()->check_chol_1();      break;
      case UHM_LU_NOPIV : e->get_matrix()->check_lu_nopiv_1();  break;
      case UHM_LU_PIV   : e->get_matrix()->check_lu_piv_1();    break;
      case UHM_QR       : e->get_matrix()->check_qr_1();        break;
      }
    }
    if (level == 2) {
      switch (type) {
      case UHM_CHOL     : e->get_matrix()->check_chol_2();      break;
      case UHM_LU_NOPIV : e->get_matrix()->check_lu_nopiv_2();  break;
      case UHM_LU_PIV   : e->get_matrix()->check_lu_piv_2();    break;
      case UHM_QR       : e->get_matrix()->check_qr_2();        break;
      }
      op_merge_rhs_r(e);
    }
    return true;
  }
    
  // --------------------------------------------------------------
  // ** operation set 

  // create matrix buffer
  bool op_create_matrix_buffer_without_schur (Element e) { return op_create_matrix_buffer( e, false ); }
  bool op_create_matrix_buffer_with_schur    (Element e) { return op_create_matrix_buffer( e, true  ); }

  // for leaves
  bool op_restore_connectivity(Element e) {
    assert(element_valid(e));
    e->restore_connectivity();
    return true;
  }

  bool op_numbering(Element e) {
    assert(element_valid(e));
    e->numbering();
    return true;
  }

  // for elements
  bool op_update_connectivity(Element e) {
    assert(element_valid(e));
    e->set_marker(0,-1);
    e->set_marker(1,-1);
    e->merge_nodes_from_children();
    e->separate_nodes();
    return true;
  }

  bool op_arrange_nodes(Element e) {
    assert(element_valid(e));
    e->arrange_nodes();
    return true;
  }

  bool op_update_generation(Element e) {
    assert(element_valid(e));
    e->update_generation();
    return true;
  }

  // --------------------------------------------------------------
  bool op_merge_full_with_free               (Element e) { return op_merge(e, 1, 1, 1, 0, 1, 1); }
  bool op_merge_full_without_free            (Element e) { return op_merge(e, 1, 1, 1, 0, 1, 0); }

  bool op_merge_full_with_ooc                (Element e) { return op_merge(e, 1, 1, 1, 0, 1, 2); }

  bool op_merge_rhs_x                        (Element e) { return op_merge(e, 0, 1, 0, 0, 0, 0); }
  bool op_merge_rhs_r                        (Element e) { return op_merge(e, 0, 0, 0, 1, 0, 0); }

  bool op_branch_rhs_x                       (Element e) { return op_branch(e, 1, 0, 0); }
  bool op_branch_rhs_r                       (Element e) { return op_branch(e, 0, 0, 1); }
  // --------------------------------------------------------------
  bool op_chol_without_merge                 (Element e) { return op_decompose(e, UHM_CHOL, 0, 0); }
  bool op_chol_with_merge_and_free           (Element e) { return op_decompose(e, UHM_CHOL, 1, 1); }
  bool op_chol_with_merge_and_no_free        (Element e) { return op_decompose(e, UHM_CHOL, 1, 0); }
  bool op_chol_with_merge_and_ooc            (Element e) { return op_decompose(e, UHM_CHOL, 1, 2); }
  // --------------------------------------------------------------
  bool op_solve_chol_1_x_without_merge       (Element e) { return op_solve(e, UHM_CHOL, 1, 1, 0, 0, 0); }
  bool op_solve_chol_1_x_with_merge          (Element e) { return op_solve(e, UHM_CHOL, 1, 1, 0, 1, 0); }
  bool op_solve_chol_2_x_without_branch      (Element e) { return op_solve(e, UHM_CHOL, 2, 1, 0, 0, 0); }
  bool op_solve_chol_2_x_with_branch         (Element e) { return op_solve(e, UHM_CHOL, 2, 1, 0, 0, 1); }

  bool op_solve_chol_1_r_without_merge       (Element e) { return op_solve(e, UHM_CHOL, 1, 0, 1, 0, 0); }
  bool op_solve_chol_1_r_with_merge          (Element e) { return op_solve(e, UHM_CHOL, 1, 0, 1, 1, 0); }
  bool op_solve_chol_2_r_without_branch      (Element e) { return op_solve(e, UHM_CHOL, 2, 0, 1, 0, 0); }
  bool op_solve_chol_2_r_with_branch         (Element e) { return op_solve(e, UHM_CHOL, 2, 0, 1, 0, 1); }
  // --------------------------------------------------------------
  bool op_check_chol_1                       (Element e) { return op_check(e, UHM_CHOL, 1); }
  bool op_check_chol_2                       (Element e) { return op_check(e, UHM_CHOL, 2); } 
  // --------------------------------------------------------------
  bool op_lu_nopiv_without_merge             (Element e) { return op_decompose(e, UHM_LU_NOPIV, 0, 0); }
  bool op_lu_nopiv_with_merge_and_free       (Element e) { return op_decompose(e, UHM_LU_NOPIV, 1, 1); }
  bool op_lu_nopiv_with_merge_and_no_free    (Element e) { return op_decompose(e, UHM_LU_NOPIV, 1, 0); }
  bool op_lu_nopiv_with_merge_and_ooc        (Element e) { return op_decompose(e, UHM_LU_NOPIV, 1, 2); }
  // --------------------------------------------------------------
  bool op_solve_lu_nopiv_1_x_without_merge   (Element e) { return op_solve(e, UHM_LU_NOPIV, 1, 1, 0, 0, 0); }
  bool op_solve_lu_nopiv_1_x_with_merge      (Element e) { return op_solve(e, UHM_LU_NOPIV, 1, 1, 0, 1, 0); }
  bool op_solve_lu_nopiv_2_x_without_branch  (Element e) { return op_solve(e, UHM_LU_NOPIV, 2, 1, 0, 0, 0); }
  bool op_solve_lu_nopiv_2_x_with_branch     (Element e) { return op_solve(e, UHM_LU_NOPIV, 2, 1, 0, 0, 1); }

  bool op_solve_lu_nopiv_1_r_without_merge   (Element e) { return op_solve(e, UHM_LU_NOPIV, 1, 0, 1, 0, 0); }
  bool op_solve_lu_nopiv_1_r_with_merge      (Element e) { return op_solve(e, UHM_LU_NOPIV, 1, 0, 1, 1, 0); }
  bool op_solve_lu_nopiv_2_r_without_branch  (Element e) { return op_solve(e, UHM_LU_NOPIV, 2, 0, 1, 0, 0); }
  bool op_solve_lu_nopiv_2_r_with_branch     (Element e) { return op_solve(e, UHM_LU_NOPIV, 2, 0, 1, 0, 1); }
  // --------------------------------------------------------------
  bool op_check_lu_nopiv_1                   (Element e) { return op_check(e, UHM_LU_NOPIV, 1); }
  bool op_check_lu_nopiv_2                   (Element e) { return op_check(e, UHM_LU_NOPIV, 2); } 
  // --------------------------------------------------------------
  bool op_lu_incpiv_without_merge            (Element e) { return op_decompose(e, UHM_LU_INCPIV, 0, 0); }
  bool op_lu_incpiv_with_merge_and_free      (Element e) { return op_decompose(e, UHM_LU_INCPIV, 1, 1); }
  bool op_lu_incpiv_with_merge_and_no_free   (Element e) { return op_decompose(e, UHM_LU_INCPIV, 1, 0); }
  bool op_lu_incpiv_with_merge_and_ooc       (Element e) { return op_decompose(e, UHM_LU_INCPIV, 1, 2); }
  // --------------------------------------------------------------
  bool op_lu_piv_without_merge               (Element e) { return op_decompose(e, UHM_LU_PIV, 0, 0); }
  bool op_lu_piv_with_merge_and_free         (Element e) { return op_decompose(e, UHM_LU_PIV, 1, 1); }
  bool op_lu_piv_with_merge_and_no_free      (Element e) { return op_decompose(e, UHM_LU_PIV, 1, 0); }
  bool op_lu_piv_with_merge_and_ooc          (Element e) { return op_decompose(e, UHM_LU_PIV, 1, 2); }
  // --------------------------------------------------------------
  bool op_solve_lu_piv_1_x_without_merge     (Element e) { return op_solve(e, UHM_LU_PIV, 1, 1, 0, 0, 0); }
  bool op_solve_lu_piv_1_x_with_merge        (Element e) { return op_solve(e, UHM_LU_PIV, 1, 1, 0, 1, 0); }
  bool op_solve_lu_piv_2_x_without_branch    (Element e) { return op_solve(e, UHM_LU_PIV, 2, 1, 0, 0, 0); }
  bool op_solve_lu_piv_2_x_with_branch       (Element e) { return op_solve(e, UHM_LU_PIV, 2, 1, 0, 0, 1); }

  bool op_solve_lu_piv_1_r_without_merge     (Element e) { return op_solve(e, UHM_LU_PIV, 1, 0, 1, 0, 0); }
  bool op_solve_lu_piv_1_r_with_merge        (Element e) { return op_solve(e, UHM_LU_PIV, 1, 0, 1, 1, 0); }
  bool op_solve_lu_piv_2_r_without_branch    (Element e) { return op_solve(e, UHM_LU_PIV, 2, 0, 1, 0, 0); }
  bool op_solve_lu_piv_2_r_with_branch       (Element e) { return op_solve(e, UHM_LU_PIV, 2, 0, 1, 0, 1); }
  // --------------------------------------------------------------
  bool op_check_lu_piv_1                     (Element e) { return op_check(e, UHM_LU_PIV, 1); }
  bool op_check_lu_piv_2                     (Element e) { return op_check(e, UHM_LU_PIV, 2); } 
  // --------------------------------------------------------------
  bool op_qr_without_merge                   (Element e) { return op_decompose(e, UHM_QR, 0, 0); }
  bool op_qr_with_merge_and_free             (Element e) { return op_decompose(e, UHM_QR, 1, 1); }
  bool op_qr_with_merge_and_no_free          (Element e) { return op_decompose(e, UHM_QR, 1, 0); }
  bool op_qr_with_merge_and_ooc              (Element e) { return op_decompose(e, UHM_QR, 1, 2); }
  // --------------------------------------------------------------
  bool op_solve_qr_1_x_without_merge         (Element e) { return op_solve(e, UHM_QR, 1, 1, 0, 0, 0); }
  bool op_solve_qr_1_x_with_merge            (Element e) { return op_solve(e, UHM_QR, 1, 1, 0, 1, 0); }
  bool op_solve_qr_2_x_without_branch        (Element e) { return op_solve(e, UHM_QR, 2, 1, 0, 0, 0); }
  bool op_solve_qr_2_x_with_branch           (Element e) { return op_solve(e, UHM_QR, 2, 1, 0, 0, 1); }

  bool op_solve_qr_1_r_without_merge         (Element e) { return op_solve(e, UHM_QR, 1, 0, 1, 0, 0); }
  bool op_solve_qr_1_r_with_merge            (Element e) { return op_solve(e, UHM_QR, 1, 0, 1, 1, 0); }
  bool op_solve_qr_2_r_without_branch        (Element e) { return op_solve(e, UHM_QR, 2, 0, 1, 0, 0); }
  bool op_solve_qr_2_r_with_branch           (Element e) { return op_solve(e, UHM_QR, 2, 0, 1, 0, 1); }
  // --------------------------------------------------------------
  bool op_check_qr_1                         (Element e) { return op_check(e, UHM_QR, 1); }
  bool op_check_qr_2                         (Element e) { return op_check(e, UHM_QR, 2); } 
  // --------------------------------------------------------------

  // --------------------------------------------------------------
  bool op_solve_lu_piv_1_x_without_merge_ooc (Element e) { 
    assert(op_create_buffer_ooc(e));
    assert(op_read_ooc(e));
    assert(op_solve(e, UHM_LU_PIV, 1, 1, 0, 0, 0));
    assert(op_free_buffer_ooc(e));
    return true;
  }
  bool op_solve_lu_piv_1_x_with_merge_ooc    (Element e) { 
    assert(op_create_buffer_ooc(e));
    assert(op_read_ooc(e));
    assert(op_solve(e, UHM_LU_PIV, 1, 1, 0, 1, 0)); 
    assert(op_free_buffer_ooc(e));
    return true;
  }
  bool op_solve_lu_piv_2_x_without_branch_ooc(Element e) { 
    assert(op_create_buffer_ooc(e));
    assert(op_read_ooc(e));
    assert(op_solve(e, UHM_LU_PIV, 2, 1, 0, 0, 0)); 
    assert(op_free_buffer_ooc(e));
    return true;
  }
  bool op_solve_lu_piv_2_x_with_branch_ooc   (Element e) { 
    assert(op_create_buffer_ooc(e));
    assert(op_read_ooc(e));
    assert(op_solve(e, UHM_LU_PIV, 2, 1, 0, 0, 1)); 
    assert(op_free_buffer_ooc(e));
    return true;
  }
  
  bool op_solve_lu_piv_1_r_without_merge_ooc (Element e) { 
    assert(op_create_buffer_ooc(e));
    assert(op_read_ooc(e));
    assert(op_solve(e, UHM_LU_PIV, 1, 0, 1, 0, 0)); 
    assert(op_free_buffer_ooc(e));
    return true;
  }
  bool op_solve_lu_piv_1_r_with_merge_ooc    (Element e) { 
    assert(op_create_buffer_ooc(e));
    assert(op_read_ooc(e));
    assert( op_solve(e, UHM_LU_PIV, 1, 0, 1, 1, 0)); 
    assert(op_free_buffer_ooc(e));
    return true;
  }
  bool op_solve_lu_piv_2_r_without_branch_ooc(Element e) { 
    assert(op_create_buffer_ooc(e));
    assert(op_read_ooc(e));
    assert(op_solve(e, UHM_LU_PIV, 2, 0, 1, 0, 0)); 
    assert(op_free_buffer_ooc(e));
    return true;
  }
  bool op_solve_lu_piv_2_r_with_branch_ooc   (Element e) { 
    assert(op_create_buffer_ooc(e));
    assert(op_read_ooc(e));
    assert(op_solve(e, UHM_LU_PIV, 2, 0, 1, 0, 1)); 
    assert(op_free_buffer_ooc(e));
    return true;
  }
  // --------------------------------------------------------------
  bool op_check_lu_piv_1_ooc                 (Element e) { 
    assert(op_create_buffer_ooc(e));
    assert(op_read_ooc(e));
    assert(op_check(e, UHM_LU_PIV, 1)); 
    assert(op_free_buffer_ooc(e));
    return true;
  }
  bool op_check_lu_piv_2_ooc                 (Element e) { 
    assert(op_create_buffer_ooc(e));
    assert(op_read_ooc(e));
    assert(op_check(e, UHM_LU_PIV, 2)); 
    assert(op_free_buffer_ooc(e));
    return true;
  } 

  // --------------------------------------------------------------
  bool op_check_solution(Element e) {
    assert(element_valid(e) && e->is_matrix_created());
    e->get_matrix()->check_solution();
    return true;
  }
  bool op_improve_solution(Element e) {
    assert(element_valid(e) && e->is_matrix_created());
    e->get_matrix()->improve_solution();
    return true;
  }
  bool op_test(Element e) {
    assert(element_valid(e));
    std::pair< int,int > n_dof = e->get_n_dof();
    printf("id : %d, (%d %d)\n", 
	   e->get_id(), n_dof.first, n_dof.second);
    return true;
  }
}
