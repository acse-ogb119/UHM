
#include "dmumps_c.h"
#include "zmumps_c.h"

#include "uhm.hxx"
#include "uhm/interf/mumps.hxx"

#define UHM_ERROR_TOL 1.0e-5

int main (int argc, char **argv)
{
  FLA_Init();

  uhm::Mesh m;
  uhm::Scheduler s;

  // input check
  if (argc != 5) {
    printf("Try : mumpstest [n_thread][decomposition][is_sym][input_file]\n");
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

  double 
    t_base, 
    t_export_matrix, 
    t_init, 
    t_analyze, 
    t_decompose, 
    t_solve, 
    t_rhs_import, 
    t_finalize;

  uhm::Mumps mumps = new uhm::Mumps_(UHM_REAL);
  int mumps_sym;
  if (is_sym) mumps_sym = 2;
  else mumps_sym = 0;

  // Step -1 : initialization
  t_base       = uhm::timer();
  mumps->set_par(1); // host working
  mumps->set_sym(mumps_sym); // 0-unsym, 1-spd, 2-general sym
  mumps->set_comm(-987654); // mpi comm world

  mumps->init();
  t_init = uhm::timer() - t_base;

  mumps->set_show_n_rhs(10);

  // Step 0 : export matrix
  t_base       = uhm::timer();
  m->export_matrix(mumps, n_rhs, is_sym);
  t_export_matrix = uhm::timer() - t_base;

  // Step 1 : setting

  // icntl setting
  mumps->set_icntl(1, 6); // stream for error message : default 6 
  mumps->set_icntl(2, 0); // stream for warning message : default 0
  mumps->set_icntl(3, 6); // stream for output message : default 0
  mumps->set_icntl(4, 3); // error level - only message printed
  mumps->set_icntl(5, 0); // 0-assembled, 1-elemental
  mumps->set_icntl(6, 7); // permutation and scaling : default  7 automatic
  mumps->set_icntl(7, 7); 
  // ordering : 0-AMD, 1-user, 2-AMF, 3-SCOTCH, 4-PORD, 5-METIS, 6-QAMD, 7-auto

  mumps->set_icntl(8, 77); // default scaling
  mumps->set_icntl(9, 1);  // solve Ax = b, otherwise A^T x = b


  mumps->set_icntl(13,0); // ScaLAPACK used in root
  mumps->set_icntl(18,0); // input matrix centralized on host

  
   // cntl setting 
  mumps->set_cntl(1, 0.01); // threshold pivot default 0.01


  // Step 2 : solver initialize
  t_base       = uhm::timer();
  mumps->analyze();
  t_analyze = uhm::timer() - t_base;

  t_base       = uhm::timer();
  mumps->decompose();
  t_decompose = uhm::timer() - t_base;

  t_base       = uhm::timer();
  mumps->solve();
  t_solve = uhm::timer() - t_base;

  t_base       = uhm::timer();
  mumps->export_matrix(m);
  t_rhs_import = uhm::timer() - t_base;

  t_base       = uhm::timer();
  mumps->finalize();
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


  delete mumps;

  m->unlock();

  delete m;

  FLA_Finalize();
  return 0;
}


