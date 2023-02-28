 program create_and_diagonalise

    use hdf5_read_m
    use blacs_arrays
    use global_variables
    use hdf5
    implicit none

    complex(kind=8), dimension(:), allocatable :: mf_v
    complex(kind=8), dimension(:), allocatable :: mf_c

    integer :: nv,nc
    integer :: n1,n2,n3,n4,n5,n6,p,q
    double precision  :: wcoul0,vol


!!!!!!
     !!!!! load global group sizes !!!!! 
!!!!! 
       !reads nc and nv from bsemat.h5 file
       !reads E_qp of all the valence and conduction bands
       ! calculates n1 , n2 (1+2) n3 (1+2+3) , n4, n5 ,n6( 1+2+3+4+5+6)
       ! 1 - tt-(h ne hp)(l ne lp) 2-ss (h ne hp)(l ne lp) 3 hh-ll 4- hhp -ll 5 hh - llp 6 h-l (single)

       call load_global_mat_size()





!!!!!!
    !!!!!!global matrix Size !!!!!!
!!!!!       
        pzheevx_vars%m = 30!sys_var%n6  !n6       ! No of rows
        pzheevx_vars%n = 30!sys_var%n6 !n6       ! No of columns
        pzheevx_vars%mb = 4!25             ! No of block row
        pzheevx_vars%nb = 4!25            ! No of block columns


!!!!!!!
        !!!!!! initiate_mpi !!!!!
!!!!!!
        !initialize  mpi
        !get mpi%rank and mpi%size_

        call  initiate_mpi()

!!!!!!!
        !!!!!! setup blacs variables and grid !!!!!
!!!!!!
        !initialie blacs
        !create processor grid with grid%nprow vs grid%ncol
        !create row and proc id to each of the processors grid%mypcol grid%myprow


        call  setup_blacs_grid()

!!!!!!!!
      !!!!!!! compute required local storage and allocate local arrays !!!!!
!!!!!!!!
       !calculate  hamiltonian%size_ = hamiltonian%locq*hamiltonian%lld
       !intializes the array descriptor desca
       !allocates the local hamiltonian%mat

       call allocate_Hamiltonian()

!!!!!!!
      !!!!! compute and distribute the matrix in a block cyclic fashion !!!! 
!!!!!!!
       call distribute_matrix()


!!!!!
    !!!!! diagonalize the matrix 
!!!!
      call diagonalize_hamiltonian()

end program create_and_diagonalise

