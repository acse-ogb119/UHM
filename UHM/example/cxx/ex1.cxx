#include "uhm.hxx"

#define UHM_ERROR_TOL 1.0e-5

int main (int argc, char **argv)
{

  // ---------------------------------------------------------
  // ** libFLAME initialization

  FLA_Init();

  int n_threads = 4, blocksize = 10;

  // ** create a mesh object 
  uhm::Mesh m = new uhm::Mesh_;

  // ** set environments
  uhm::set_num_threads     (n_threads);
  uhm::set_hier_block_size (blocksize);


  // ---------------------------------------------------------
  // ** mesh connectivity
  //           5    11   6
  //            +-------+ 
  //           / \  2  /
  //        8 /  9\   / 10
  //         /  1  \ /
  //        +-------+
  //       3    7    4

  int elt[2][7]     = { { 3, 4, 5, /**/   7, 9, 8, /**/  1 },
                        { 4, 6, 5, /**/  10,11, 9, /**/  2 } };
  int elt_n_dof[7]  =   { 1, 1, 1, /**/   3, 3, 3, /**/  3 };

  // ** mirroring the nodal connectivity
  uhm::Element e[2];

  for (int j=0;j<2;++j) {
    
    // ** add new element
    e[j] = m->add_element();

    for (int i=0;i<7;++i) {

      int nod       = elt[j][i];
      int n_dof_nod = elt_n_dof[i];

      // ** add new node and assign the nodes to element
      e[j]->add_node( m->add_node(nod, n_dof_nod) );

    }
  }

  // ---------------------------------------------------------
  // ** analysis

  uhm::build_tree(m);
  printf( "New mesh info :: < n_nodes %d, n_elements %d >\n", 
	  m->get_n_nodes(), m->get_n_elements() );

  // ** prevent further mesh modification
  m->lock();

  unsigned int n_dof_mesh = m->get_n_dof();

  int datatype  = UHM_REAL;
  int n_dof_elt = 15;
  int n_rhs     = 1;
  
  // ** workspace for shcur complements will be created dynamically
  int is_schur  = false;   

  // ** create symbolic matrices then allocate buffer
  m->create_matrix_without_buffer ( datatype, n_rhs );
  m->create_matrix_buffer         ( is_schur );

  // ---------------------------------------------------------
  // ** interface to unassembled matrices : copy_in

  linal::Flat_ A, B;
  A.create( datatype, n_dof_elt, n_dof_elt );
  B.create( datatype, n_dof_elt, n_rhs );

  FLA_Random_matrix( ~A );
  FLA_Random_matrix( ~B );
  
  for (int j=0;j<2;++j) {
    m->copy_in(e[j], datatype, n_dof_elt, n_dof_elt,
               &(elt[j][0]), UHM_LHS, (void*)A.get_buffer() );
    m->copy_in(e[j], datatype, n_dof_elt, n_rhs,
               &(elt[j][0]), UHM_RHS, (void*)B.get_buffer() );

    char title[32];
    sprintf(title, "- elt %d : unassembled b -", j);
    B.disp(title);
  }

  // ** make RHS ready
  m->set_rhs();

  // ---------------------------------------------------------
  // ** Factorization, Solution and Check

  m->lu_piv_with_free();
  m->solve_lu_piv();
  m->check_lu_piv();

  double residual = m->get_residual();

  // ---------------------------------------------------------
  // ** Interface to unassembled matrices : copy_out
  for (int j=0;j<2;++j) {
    m->copy_out(e[j], datatype, n_dof_elt, n_rhs,
                &(elt[j][0]), UHM_RHS, (void*)B.get_buffer() );
    char title[32];
    sprintf(title, "- elt %d : solution -", j);
    B.disp(title);
  }

  A.free();
  B.free();

  // ** unlock the mesh
  m->unlock();

  // ----------------------------------------------------------------
  // ** report

  printf("==== Report =====\n");
  printf("Number of RHS          = %d\n", n_rhs);
  printf("Number of threads      = %d\n", uhm::get_num_threads());
  printf("NDOF                   = %d\n", n_dof_mesh);
  printf("--------------------------\n");

  if (uhm::is_multithreading_enable())
    printf("Openmp is used\n");
  else
    printf("Openmp is NOT used\n");

  if (uhm::is_hier_matrix_enable()) 
    printf("Hier-Matrix is used ( blocksize = %d )\n", 
	   uhm::get_hier_block_size());
  else
    printf("Hier-Matrix is NOT used \n");

  printf("Residual              = %E\n", residual);
  if ( residual < UHM_ERROR_TOL ) 
    printf("TESTING UHM : **** PASS **** \n");
  else 
    printf("TESTING UHM : **** FAIL **** \n");

  // ** delete a mesh object
  delete m;

  // ** libFLAME finalization
  FLA_Finalize();

  return 0;
}


