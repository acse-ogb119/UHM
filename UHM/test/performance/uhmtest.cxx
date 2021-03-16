#include "uhm.hxx"

#define UHM_ERROR_TOL 1.0e-5

int main (int argc, char **argv)
{
  FLA_Init();

  uhm::Mesh m;

  // input check
  if (argc != 5) {
    printf("Try : uhm [n_thread][decomposition][blocksize][input_file]\n");
    return 0;
  }

  int n_threads, decomposition, blocksize, svd_cutoff;
  double rel_thres;
  char *filename;
  n_threads     = atoi( (argv[1]) );
  decomposition = atoi( (argv[2]) );
  blocksize     = atoi( (argv[3]) );
  filename      = argv[4];

  double t_base, t_tmp, t_build_tree, t_decompose, t_solve;
  double f_decompose, f_solve, m_estimate, m_used, m_max_used;

  printf( "BEGIN : Import mesh from file : %s\n", filename );
  m = new uhm::Mesh_;
  m->import_file( filename );

  m->get_n_dof();
  printf( "END   : Import \n" );

  uhm::set_num_threads(n_threads);

  printf( "BEGIN : Build tree \n" );
  t_base       = uhm::timer();
  uhm::build_tree(m);
  t_build_tree = uhm::timer() - t_base;
  printf( "END   : Build tree\n" );
  printf( "New mesh info < n_nodes %d, n_elements %d >\n", 
	  m->get_n_nodes(), m->get_n_elements() );

  m->lock();

  int datatype = UHM_REAL; 
  int n_rhs    = 1;
  int is_schur = false;
  uhm::set_hier_block_size(blocksize);

  printf( "BEGIN : Create matrix without buffer \n" );
  m->create_matrix_without_buffer( datatype, n_rhs );
  printf( "END   : Create matrix without buffer \n" );

  printf( "BEGIN : Create buffer \n" );
  m->create_matrix_buffer(is_schur);
  printf( "END   : Create buffer \n" );

  switch (decomposition) {
  case UHM_CHOL: 
    m->random_spd_matrix(); 
    m->triangularize();
    break;
  default: 
    m->random_matrix(); 
    break;
  }

  m->set_rhs();

  // ----------------------------------------------------------------
  t_decompose = 1.0e9;
  t_solve     = 1.0e9;

  switch (decomposition) {
  case UHM_CHOL:
    printf("BEGIN : CHOL Decomposition\n");
    t_base      = uhm::timer();
    m->chol_with_free();
    t_tmp = uhm::timer() - t_base;
    t_decompose = min(t_tmp, t_decompose);
    printf("END   : CHOL Decomposition\n");
    break;
  case UHM_LU_NOPIV:
    printf("BEGIN : LU_NOPIV Decomposition\n");
    t_base      = uhm::timer();
    m->lu_nopiv_with_free();
    t_tmp = uhm::timer() - t_base;
    t_decompose = min(t_tmp, t_decompose);
    printf("END   : LU_NOPIV Decomposition\n");
    break;
  case UHM_LU_INCPIV:
    printf("BEGIN : LU_INCPIV Decomposition\n");
    t_base      = uhm::timer();
    m->lu_incpiv_with_free();
    t_tmp = uhm::timer() - t_base;
    t_decompose = min(t_tmp, t_decompose);
    printf("END   : LU_INCPIV Decomposition\n");
    break;
  case UHM_LU_PIV:
    printf("BEGIN : LU_PIV Decomposition\n");
    t_base      = uhm::timer();
    m->lu_piv_with_free();
    t_tmp = uhm::timer() - t_base;
    t_decompose = min(t_tmp, t_decompose);
    printf("END   : LU_PIV Decomposition\n");
    break;
  case UHM_QR:
    printf("BEGIN : QR Decomposition\n");
    t_base      = uhm::timer();
    m->qr_with_free();
    t_tmp = uhm::timer() - t_base;
    t_decompose = min(t_tmp, t_decompose);
    printf("END   : QR Decomposition\n");
    break;
  }

  m_used      = uhm::matrix_buffer_used();
  m_max_used  = uhm::matrix_max_buffer_used();
  
  switch (decomposition) {
  case UHM_CHOL:
    printf("BEGIN : CHOL Solve\n");
    t_base  = uhm::timer();
    m->solve_chol();
    t_tmp = uhm::timer() - t_base;
    t_solve = min(t_tmp, t_solve);
    printf("END   : CHOL Solve\n");
    break;
  case UHM_LU_NOPIV:
    printf("BEGIN : LU_NOPIV Solve\n");
    t_base  = uhm::timer();
    m->solve_lu_nopiv();
    t_tmp = uhm::timer() - t_base;
    t_solve = min(t_tmp, t_solve);
    printf("END   : LU_NOPIV Solve\n");
    break;
  case UHM_LU_INCPIV:
  case UHM_LU_PIV:
    printf("BEGIN : LU_PIV Solve\n");
    t_base  = uhm::timer();
    m->solve_lu_piv();
    t_tmp = uhm::timer() - t_base;
    t_solve = min(t_tmp, t_solve);
    printf("END   : LU_PIV Solve\n");
    break;
  case UHM_QR:
    printf("BEGIN : QR Solve\n");
    t_base  = uhm::timer();
    m->solve_qr();
    t_tmp = uhm::timer() - t_base;
    t_solve = min(t_tmp, t_solve);
    printf("END   : QR Solve\n");
    break;
  }
  
  switch (decomposition) {
  case UHM_CHOL:
    printf("BEGIN : CHOL Check\n");
    m->check_chol();
    printf("END   : CHOL Check\n");
    break;
  case UHM_LU_NOPIV:
    printf("BEGIN : LU_NOPIV Check\n");
    m->check_lu_nopiv();
    printf("END   : LU_NOPIV Check\n");
    break;
  case UHM_LU_INCPIV:
  case UHM_LU_PIV:
    printf("BEGIN : LU_PIV Check\n");
    m->check_lu_piv();
    printf("END   : LU_PIV Check\n");
    break;
  case UHM_QR:
    printf("BEGIN : QR Check\n");
    m->check_qr();
    printf("END   : QR Check\n");
    break;
  }

  // ----------------------------------------------------------------
  unsigned int n_dof     = m->get_n_dof();
  unsigned int n_nonzero_factor;
  
  m->estimate_cost( decomposition, datatype, n_rhs,
                    f_decompose, f_solve, n_nonzero_factor, m_estimate );

  printf("==== Report =====\n");
  printf("Number of RHS          = %d\n", n_rhs);
  printf("Number of threads      = %d\n", uhm::get_num_threads());
  printf("Decomposition          = %d\n", decomposition);
  printf("CHOL(1), LU_NOPIV(2), LU_PIV(3), LU_INCPIV(4), QR(5)\n");
  printf("Blocksize              = %d\n", uhm::get_hier_block_size());
  printf("NDOF                   = %d\n", n_dof);
  printf("Non-zero entries       = %d\n", n_nonzero_factor);
  printf("Sparsity               = %E\n", 
	 (double)n_nonzero_factor/(double)(n_dof*n_dof));
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

  printf("--------------------------\n");
  printf("Time build tree (s)   = %E\n", t_build_tree);
  printf("Time tree / Time decom= %E\n", t_build_tree/t_decompose);
  printf("--------------------------\n");
  printf("Buffer estimate (MB)  = %6.0lf\n", m_estimate/1.0e6);
  printf("Buffer used     (MB)  = %6.0lf\n", m_used/1.0e6);
  printf("Max buffer used (MB)  = %6.0lf\n", m_max_used/1.0e6);
  printf("--------------------------\n");
  printf("Time decom (s)        = %E\n", t_decompose);
  printf("Time solve (s)        = %E\n", t_solve);
  printf("--------------------------\n");
  printf("FLOP decom (GFLOP)    = %6.3lf\n", f_decompose/1.0e9);
  printf("FLOP solve (GFLOP)    = %6.3lf\n", f_solve/1.0e9);
  printf("--------------------------\n");
  printf("FLOP decom (GFLOP/s)  = %6.3lf\n", f_decompose/t_decompose/1.0e9);
  printf("FLOP solve (GFLOP/s)  = %6.3lf\n", f_solve/t_solve/1.0e9);
  printf("--------------------------\n");
  switch (decomposition) {
  case UHM_CHOL:
  case UHM_LU_NOPIV:
  case UHM_LU_PIV:
  case UHM_LU_INCPIV:
  case UHM_QR:

    double residual = m->get_residual();

    printf("Residual              = %E\n", residual);
    if ( residual < UHM_ERROR_TOL ) 
      printf("TESTING UHM : **** PASS **** \n");
    else 
      printf("TESTING UHM : **** FAIL **** \n");
  }

  m->unlock();

  delete m;

  FLA_Finalize();
  return 0;
}


