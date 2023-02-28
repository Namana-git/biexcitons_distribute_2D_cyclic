module blacs_arrays

     use global_variables
     use create_H
     use hdf5_read_m
     include 'mpif.h'

     
     public :: setup_blacs_grid, allocate_Hamiltonian,initiate_mpi,distribute_matrix,&
               load_global_mat_size
 

     contains

     subroutine setup_blacs_grid()

         call blacs_pinfo(grid%rank, grid%size_)

         grid%nprow=int(sqrt(real(grid%size_)))
         grid%npcol = grid%nprow


         call blacs_get(-1,0,grid%context)
         call blacs_gridinit(grid%context,'Col-major',grid%nprow,grid%npcol)
         call blacs_gridinfo(grid%context,grid%nprow,grid%npcol,grid%myprow,grid%mypcol)
         
         print*,grid%nprow,grid%npcol,grid%myprow,grid%mypcol



     end subroutine setup_blacs_grid

     subroutine allocate_Hamiltonian()
         integer :: rsrc, csrc, info, AllocateStatus
         integer, external :: numroc

         rsrc = 0
         csrc = 0


         hamiltonian%locq = numroc(pzheevx_vars%m, pzheevx_vars%mb, grid%mypcol, csrc, grid%npcol)
         hamiltonian%locq = max(hamiltonian%locq,1)
         hamiltonian%lld = numroc(pzheevx_vars%n, pzheevx_vars%nb, grid%myprow, rsrc, grid%nprow)
         hamiltonian%lld = max(hamiltonian%lld,1)

         hamiltonian%size_ = hamiltonian%locq*hamiltonian%lld

         print*,"hamiltonian%locq,hamiltonian%lld ,hamiltonian%size_",hamiltonian%locq,hamiltonian%lld,hamiltonian%size_

         call descinit(hamiltonian%desca, pzheevx_vars%m, pzheevx_vars%n, pzheevx_vars%mb, &
                  pzheevx_vars%nb, rsrc, csrc, grid%context, hamiltonian%lld, info)             ! intializes the array descriptor desca

      
        allocate(hamiltonian%mat(Hamiltonian%size_), STAT = AllocateStatus)
        if (AllocateStatus.ne.0) STOP "Insufficient memory for A"

        evec%locq = numroc(pzheevx_vars%m, pzheevx_vars%mb, grid%mypcol, csrc, grid%npcol)
        evec%locq = max(evec%locq,1)
        evec%lld = numroc(pzheevx_vars%n, pzheevx_vars%nb, grid%myprow, rsrc, grid%nprow)
        evec%lld = max(evec%lld,1)

        evec%size_ = evec%locq*evec%lld
        call descinit(evec%desca, pzheevx_vars%m, pzheevx_vars%n, pzheevx_vars%mb, &
                  pzheevx_vars%nb, rsrc, csrc, grid%context, evec%lld, info)

       
        allocate(evec%mat(evec%size_))
        allocate(eval(pzheevx_vars%m)) 


     end subroutine allocate_Hamiltonian

     subroutine initiate_mpi()
         integer :: ierror 

         call MPI_INIT(ierror)
         call MPI_COMM_RANK(MPI_COMM_WORLD,mpi%rank,ierror)
         call MPI_COMM_SIZE(MPI_COMM_WORLD,mpi%size_,ierror)




     end subroutine initiate_mpi

     subroutine load_global_mat_size()
         
         integer :: nv,nc
         integer :: n1,n2,n3,n4,n5,n6,p,q

         call read_inp(sys_var%nv,sys_var%nc,sys_var%vol)
         allocate(sys_var%mf_v(sys_var%nv))
         allocate(sys_var%mf_c(sys_var%nc))
         call load_mf(sys_var%nv,sys_var%nc,sys_var%mf_v,sys_var%mf_c)
         call load_weff(sys_var%wcoul0)
         nv = sys_var%nv
         nc = sys_var%nc
         n1= nc*(nc-1)*nv*(nv-1)/4     !size of 1
         n2 = 2*n1  !size of 1 + 2
         p = nc*(nc-1)*nv/2
         q = nv*(nv-1)*nc/2
         n3 = (2*n1) + (nv*nc)        !size of 1+2+3
         n4 = (2*n1) + (nv*nc) + q    !size of 1+2+3+4
         n5 = (2*n1) + (nv*nc) + q + p !size of 1+2+3+4+5
         n6 = (2*n1) + (2*(nv*nc)) + q + p !size of 1+2+3+4+5+6

         sys_var%n1 = n1
         sys_var%n2 = n2 
         sys_var%n3 = n3
         sys_var%n4 = n4
         sys_var%n5 = n5
         sys_var%n6 = n6
         sys_var%p = p
         sys_var%q = q

         !print*,"sys_var%nv,sys_var%nc,sys_var%vol",sys_var%nv,sys_var%nc,sys_var%vol
        
         !print*,sys_var%n1,sys_var%n2,sys_var%n3,sys_var%n4,sys_var%n5,sys_var%n6





                 
     end subroutine load_global_mat_size

     subroutine distribute_matrix()
        complex(kind=8) , dimension(13) :: H1
        integer                         :: j

        !do i = 1,13 
         !  H1(i) = cmplx(real(i+7+mpi%rank),0.0)
           !print*,"i,H1(i)",i,H1(i)

        !end do
       ! print*,"i,H1(i)",i,H1(i)
 
 

        !call distrubute_rows_2D_cyclo(H1,13,1,8)

               call create_H1()
               !call gather_rows_test()
               !print*,grid%myprow,grid%mypcol,mpi%rank,"grid%myprow,grid%mypcol,mpi%rank"

               !do j = 1,hamiltonian%size_
                  !if (grid%myprow == 0 .and. grid%mypcol == 0)then
                      !print*,"j,hamiltonian%mat(j)",j,hamiltonian%mat(j)

                  !end if

               !end do

     end subroutine distribute_matrix




subroutine diagonalize_hamiltonian()

  
   

    implicit none

    integer :: ia, ja
    double precision, allocatable, dimension(:) :: gap, rwork
    double complex, allocatable, dimension(:) :: work
    integer, allocatable, dimension(:) :: iwork, ifail, iclustr
    integer :: lwork, liwork, lrwork, lgap, lifail, liclustr, info
    integer :: i
    

    ia = 1
    ja = 1
    lgap = grid%nprow*grid%npcol
    lifail = pzheevx_vars%m
    liclustr = 2*lgap

    allocate(gap(lgap))
    allocate(ifail(lifail))
    allocate(iclustr(liclustr))
    lwork = -1
    lrwork = -1
    liwork = -1
    allocate(work(1))
    allocate(rwork(1))
    allocate(iwork(1))


    call pzheevx('V', 'I', 'U', pzheevx_vars%m,  &
    hamiltonian%mat, ia, ja, hamiltonian%desca, pzheevx_vars%vl,pzheevx_vars%vu, &
    1, 10, pzheevx_vars%abstol, pzheevx_vars%comp_num_eval, &
    pzheevx_vars%comp_num_evec, eval, pzheevx_vars%orfac, evec%mat, 1, 1, evec%desca, work, &
    lwork, rwork, lrwork, iwork, liwork, ifail, iclustr, gap, info)

    lwork  = int(abs(work(1)))+1
    lrwork = int(abs(rwork(1)))+1
    liwork = int(abs(iwork(1)))+1

    deallocate(work)
    deallocate(rwork)
    deallocate(iwork)

    allocate(work(lwork))
    allocate(rwork(lrwork))
    allocate(iwork(liwork))

    call pzheevx('V', 'I', 'U', pzheevx_vars%m,   &
    hamiltonian%mat, ia, ja, hamiltonian%desca, pzheevx_vars%vl, pzheevx_vars%vu, &
    1, 10, pzheevx_vars%abstol, pzheevx_vars%comp_num_eval, &
    pzheevx_vars%comp_num_evec,eval, pzheevx_vars%orfac, evec%mat, 1, 1, evec%desca, work, &
    lwork, rwork, lrwork, iwork, liwork, ifail, iclustr, gap, info)





    deallocate(work)
    deallocate(iwork)
    deallocate(rwork)
    deallocate(gap)
    deallocate(ifail)
    deallocate(iclustr)

    ! print eigenvalues
    if (mpi%rank == 0 ) then
       do i=1,10
          print*, eval(i)*13.605698059
       enddo
    end if
    
    deallocate(eval)



end subroutine diagonalize_hamiltonian




end module blacs_arrays

