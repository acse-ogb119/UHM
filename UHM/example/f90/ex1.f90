
program main_ex1
  implicit none
  include "uhm/wrapper/fort.f90"

  integer, parameter :: &
       n_dof_elt = 15, &           ! n_dof for each element
       n_rhs = 1, &                ! n_rhs for the linear system
       datatype = UHM_REAL, &      ! datatype
       is_schur = 0                ! do not create workspace initially

  ! ** nodal connectivity
  integer, dimension( 7, 2 ), parameter :: &
       elt = reshape( (/3,4,5,  7, 9,8,  1, &
                        4,6,5, 10,11,9,  2/), (/7,2/) )
  ! ** element dofs
  integer, dimension( 7 ), parameter :: &
       elt_n_dof = (/1,1,1,  3,3,3,  3/)

  ! ** n_threads : for linux, 'cat /proc/cpuinfo'
  ! ** blocksize : for most modern computer architecture, 192 or 256
  integer :: n_threads = 4, blocksize = 10

  ! ** bind c_pointer
  integer( kind = 8 ) :: mesh, e( 2 ), n

  ! ** temporary variables
  integer :: &
       n_dof_mesh, n_dof_nod,  &
       nod, i, j, k, dummy, iflag

  ! ** elemental matrices
  real*8 :: A(n_dof_elt, n_dof_elt), B(n_dof_elt, n_rhs), res

  !---------------------------------------------------------
  ! ** libFLAME initialization

  call uhm_initialize_fla

  ! ** create a mesh object 
  call uhm_mesh_create         (mesh);

  ! ** set environments
  call uhm_set_num_threads     (n_threads);
  call uhm_set_hier_block_size (blocksize);

  ! ---------------------------------------------------------
  ! ** mesh connectivity
  !           5    11   6
  !            +-------+ 
  !           / \  2  /
  !        8 /  9\   / 10
  !         /  1  \ /
  !        +-------+
  !       3    7    4

  ! ** mirroring the nodal connectivity
  do j=1, 2

     ! ** add new element
     call uhm_add_element( mesh, e( j ) )

     do i=1, 7

        nod       = elt( i, j )
        n_dof_nod = elt_n_dof( i )
        dummy     = 0

        ! ** add new node and assign the nodes to element
        call uhm_add_node( mesh, nod, dummy, n_dof_nod, n_dof_nod, n ) 
        call uhm_element_add_node( e( j ), n );

     end do

  end do

  ! ---------------------------------------------------------
  ! ** analysis

  call uhm_build_tree( mesh );

  ! ** prevent further mesh modification
  call uhm_lock( mesh )
  call uhm_get_n_dof( mesh, n_dof_mesh )

  ! ** create symbolic matrices then allocate buffer
  call uhm_create_matrix_without_buffer ( mesh, datatype, n_rhs );
  call uhm_create_element_matrix_buffer ( mesh, is_schur );

  ! ---------------------------------------------------------
  ! ** interface to unassembled matrices : copy_in

  call srand( 100 )

  do j=1, n_dof_elt
     do i=1, n_dof_elt
        A( i , j ) = rand()
     end do
  end do

  do j=1, n_rhs
     do i=1, n_dof_elt
        B( i , j ) = rand()
     end do
  end do

  do j=1, 2
     call uhm_copy_in( &
          mesh, e( j ), datatype, n_dof_elt, n_dof_elt, &
          UHM_PHYSICS_SINGLE, elt( : , j ), UHM_LHS, A );
     call uhm_copy_in( &
          mesh, e( j ), datatype, n_dof_elt, n_rhs, &
          UHM_PHYSICS_SINGLE, elt( : , j ), UHM_RHS, B );

     write(*,*) ' '
     write(*,*) ' element = ', j
     do k=1, n_dof_elt
        write(*,*) 'b = ', B( k , : )
     end do
  end do

  ! ** make RHS ready
  call uhm_set_rhs(mesh);

  ! ---------------------------------------------------------
  ! ** Factorization, Solution and Check

  call uhm_lu_piv_with_free (mesh);
  call uhm_solve_lu_piv     (mesh);
  call uhm_check_lu_piv     (mesh);

  call uhm_get_residual     (mesh, res);

  ! ---------------------------------------------------------
  ! ** Interface to unassembled matrices : copy_out
  do j=1, 2
     call uhm_copy_out( &
          mesh, e( j ), datatype, n_dof_elt, n_rhs, &
          UHM_PHYSICS_SINGLE, elt( : , j ), UHM_RHS, B );

     write(*,*) ' '
     write(*,*) ' element = ', j
     do k=1, n_dof_elt
        write(*,*) 'solution = ', B( k, : )
     end do

  end do

  ! ** unlock the mesh
  call uhm_unlock( mesh );

  ! ----------------------------------------------------------------
  ! ** report

  write(*,*) '==== Report ====='
  write(*,*) 'Number of RHS          = ', n_rhs
  write(*,*) 'Number of threads      = ', n_threads
  write(*,*) 'NDOF                   = ', n_dof_mesh
  write(*,*) '--------------------------'

  call uhm_is_multithreading_enable(iflag)
  
  if (iflag.eq.1) then
     write(*,*) 'Openmp is used'
  else
     write(*,*) 'Openmp is NOT used'
  end if

  call uhm_is_hier_matrix_enable( iflag )
  call uhm_get_hier_block_size  ( blocksize )
  if (iflag.eq.1) then
     write(*,*) 'Hier-Matrix is used ( blocksize =  )', blocksize 
  else
     write(*,*) 'Hier-Matrix is NOT used ';
  end if

  write(*,*) 'Residual              = ', res

  ! ** delete a mesh object
  call uhm_mesh_delete(mesh)

  ! ** libFLAME finalization
  call uhm_finalize_fla

end program main_ex1
