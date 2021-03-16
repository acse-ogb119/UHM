#include "uhm.hxx"
#include "uhm/interf/pardiso.hxx"

#define UHM_ERROR_TOL 1.0e-5

int main (int argc, char **argv)
{
  FLA_Init();

  uhm::Mesh m;
  uhm::Scheduler s;

  // input check
  if (argc != 5) {
    printf("Try : uhm [n_thread][decomposition][is_sym][input_file]\n");
    return 0;
  }

  int n_threads, decomposition, is_sym;
  char *filename;

  n_threads     = atoi( (argv[1]) );
  decomposition = atoi( (argv[2]) );
  is_sym        = atoi( (argv[3]) );
  filename      = argv[4];

  printf( "BEGIN : Import mesh from file : %s\n", filename );
  m = new uhm::Mesh_;
  m->import_file( filename );
  printf( "END   : Import \n" );

  uhm::set_num_threads(n_threads);

  printf( "BEGIN : Build tree \n" );
  uhm::build_tree(m);
  printf( "END   : Build tree\n" );
  printf( "New mesh info < n_nodes %d, n_elements %d >\n", 
	  m->get_n_nodes(), m->get_n_elements() );

  m->lock();

  int datatype = UHM_REAL; 
  int n_rhs    = 1;

  uhm::set_hier_block_size(256);

  printf( "BEGIN : Create matrix without buffer \n" );
  m->create_matrix_without_buffer( datatype, n_rhs );
  printf( "END   : Create matrix without buffer \n" );

  printf( "BEGIN : Create buffer \n" );
  m->create_leaf_matrix_buffer();
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
  uhm::Pardiso pardiso = new uhm::Pardiso_(UHM_REAL);

  double 
    t_base, 
    t_export_matrix, 
    t_init, 
    t_analyze, 
    t_decompose, 
    t_solve, 
    t_rhs_import, 
    t_finalize;

  // Step 0 : export matrix
  t_base       = uhm::timer();
  m->export_matrix(pardiso, n_rhs, is_sym);
  t_export_matrix = uhm::timer() - t_base;

  // Step 1 : setting
  for (int i=1;i<65;++i)
    pardiso->set_iparm(i, 0);

  int iparm_8[2] = { 13, 8 };
  int iparm_11[2] = { 1, 0};
    
  pardiso->set_iparm(1, 1); // 0 - all default, not 0 - supply all
  pardiso->set_iparm(2, 2); // 0 - MMD, 2 - ND METIS
  pardiso->set_iparm(3, uhm::get_num_threads()); // num of threads
  pardiso->set_iparm(7,0); // num of iterative refinement
  pardiso->set_iparm(10, iparm_8[is_sym]); // 13 - non sym, 8 sym indefinite
  pardiso->set_iparm(11, iparm_11[is_sym]); // 1 - non sym, 0 - sym
  pardiso->set_iparm(13, 1); // 1 - normal matching, 2 - advanced
  pardiso->set_iparm(21, 1); // 0 - 1x1 pivot, 1 - 2x2 bunch-kaufman pivot
  pardiso->set_iparm(24, 1); // 0 - one level parallel, 1 - two level parallel
  pardiso->set_iparm(25, 1); // 0 - seq solve, 1 - par solve
  pardiso->set_iparm(52, 1); // num of dist solver : multi thread only should be 1

  pardiso->set_show_n_rhs(10);
    
  // Step 2 : solver initialize
  t_base       = uhm::timer();
  pardiso->init();
  t_init = uhm::timer() - t_base;

  t_base       = uhm::timer();
  pardiso->analyze();
  t_analyze = uhm::timer() - t_base;

  t_base       = uhm::timer();
  pardiso->decompose();
  t_decompose = uhm::timer() - t_base;

  t_base       = uhm::timer();
  pardiso->solve();
  t_solve = uhm::timer() - t_base;

  t_base       = uhm::timer();
  pardiso->export_matrix(m);
  t_rhs_import = uhm::timer() - t_base;

  t_base       = uhm::timer();
  pardiso->finalize();
  t_finalize = uhm::timer() - t_base;

  printf("==== Report =====\n");
  printf("Number of RHS          = %d\n", n_rhs);
  printf("Number of threads      = %d\n", uhm::get_num_threads());
  printf("Time export (s)        = %E\n", t_export_matrix);
  printf("Time init (s)          = %E\n", t_init);
  printf("Time analyze (s)       = %E\n", t_analyze);
  printf("Time decom (s)         = %E\n", t_decompose);
  printf("Time solve (s)         = %E\n", t_solve);
  printf("Time rhs import (s)    = %E\n", t_rhs_import);
  printf("Time finalize (s)      = %E\n", t_finalize);
  printf("--------------------------\n");


  delete pardiso;

  m->unlock();

  delete m;

  FLA_Finalize();
  return 0;
}


