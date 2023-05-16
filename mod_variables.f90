module global_variables



    type blacs_info
        integer  :: context
        integer  :: nprow, npcol, myprow, mypcol, rank, size_
    end type blacs_info

    type(blacs_info) :: grid
   type real_arrays
        double precision, allocatable, dimension(:) :: mat
        integer , dimension(9) :: desca
        integer :: size_, lld, locq
    end type real_arrays


    type(real_arrays) :: hamiltonian,evec


    type scalapack_variables
        integer  :: mb, nb, il, iu, comp_num_eval, comp_num_evec,m,n
        double precision :: vl, vu, abstol, orfac
        character(1) :: range_, comp_evec
    end type scalapack_variables

    type global_parameters 
            integer :: n1,n2,n3,n4,n5,n6,p,q,nc,nv
            double precision :: vol, wcoul0
            complex(kind=8), dimension(:), allocatable :: mf_v
            complex(kind=8), dimension(:), allocatable :: mf_c

    end type global_parameters

    type(global_parameters) :: sys_var

    type(scalapack_variables) :: pzheevx_vars

    type mpi_var
        integer  :: context
        integer  :: rank, size_,group_size,row_group
    end type mpi_var

    type(mpi_var) :: mpi



    integer , parameter :: BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,     &
                          CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,    &
                          RSRC_ = 7, CSRC_ = 8, LLD_ = 9

    double precision, allocatable, dimension(:) :: eval

















end module global_variables
