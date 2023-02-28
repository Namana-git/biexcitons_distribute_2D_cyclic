module create_H

  use hdf5
  use global_variables
  use hdf5_read_m
  include 'mpif.h'
  


  private

  public :: create_H1,create_H2,create_H12,create_H3,create_H13,&
            create_H23,create_H4,create_H14,create_H24,create_H34,&
            create_H5,create_H15,create_H25,create_H35,create_H45,&
            create_H6,create_H16,create_H26,create_H36,create_H46,create_H56
  contains


  subroutine create_H1()
      !The ij pairs are distributed amonh mpi%size_ processors.
      ! if proc_1 gets ij_1 no of ijpairs. The no of rows that proc_1 calculates is ij_1*(nv)*(nv -1)/2

      double precision :: fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16


      integer :: nb,i,j,alpha,beta,gama,eta,m,n,imatrix,keyword,c1,c2,c,x,A,y,max_ijp,q,r      
      double precision            :: element
      complex(kind=8),dimension(:),allocatable                            :: H1
      complex(kind=8),dimension(:,:),allocatable                            :: H_b

      complex(kind=8), allocatable :: bse_mat_e(:,:,:,:)
      complex(kind=8), allocatable :: bse_mat_h(:,:,:,:)
      complex(kind=8), allocatable :: bse_mat(:,:,:,:)
      integer :: ierror  
      integer ::  chunk_size,itr_length
      integer , dimension(:,:) , allocatable :: ij_mat
      integer, dimension(:),allocatable :: row_ij
      integer :: cij,cab,cijp,c1p,p_contrib
      
      



       keyword = 1
       fac_hd=0
       fac_hx=0
       fac_ed=0
       fac_ex=0
      fac_1=0
      fac_2=0
      fac_3=0
      fac_4=0
      fac_5=0
      fac_6=0
      fac_7=0
      fac_8=0
      fac_9=0
      fac_10=0
      fac_11=0
      fac_12=0
      fac_13=0
      fac_14=0
      fac_15=0
      fac_16=0

      itr_length =sys_var%nc*(sys_var%nc-1)/2  !total no of ij pairs 
      call calculate_chunk_size(chunk_size) !no of ij pairs each processor holds 
      q = int(itr_length/mpi%size_)
      r = mod(itr_length,mpi%size_)
      if (r == 0)then
         max_ij=q
         print*,"tne number of proc divide the exactly the no ijpairs" 
      else
         max_ij=q+1
         print*,"tne number of proc donot divide the exactly the no ijpairs"
      end if
      !print*,"mpi%rank,max_ij",mpi%rank,max_ij
     
    
      allocate(bse_mat_h(sys_var%nv,sys_var%nv,sys_var%nv,sys_var%nv))
      allocate(bse_mat_e(sys_var%nc,1,1,sys_var%nc))
      allocate(bse_mat(sys_var%nv,sys_var%nv,2,sys_var%nc))

      allocate(row_ij(chunk_size))
      A = sys_var%nv*(sys_var%nv-1)/2
      allocate(ij_mat(2,chunk_size))
      allocate(H1(sys_var%n1))
      H1 = cmplx(0.0,0.0)
      !print*,"myprow,mypcol,mpi%rank",grid%myprow,grid%mypcol,mpi%rank
      call distribute_ij(mpi%rank,sys_var%nc,A,chunk_size,ij_mat,row_ij)
      c1=0
      c2=0
      cij=0
      cijp=0
      cab=0
      c1p = 0
      p_contrib=0
      !print*,"c1,c2",c1,c2


         
      do x = 1,max_ij ! ij pairs are distributed among processors
        
         !if (x < chunk_size+1) then   
            !p_contrib =1
            !print*,"mpi%rank,x,chunk_size,p_contrib",mpi%rank,x,chunk_size,p_contrib
         !end if
         cab = 0 

         do beta = 1,sys_var%nv-1
             do alpha = beta+1,sys_var%nv
               
             
               !print*,"c1,c2",c1,c2
               H1 = cmplx(0.0,0.0)
               if(x < chunk_size+1) then 
                 
                   !print*,"mpi%rank,x,chunk_size,p_contrib,",mpi%rank,x,chunk_size,p_contrib
                  i = ij_mat(1,x)
                  j = ij_mat(2,x)
                  cij = row_ij(x)
                  cab = cab+1
                  c1 = cij+cab
                   !print*,"hey,i am at the beggining",c1,x
                  !print*,"mpi%rank,x,chunk_size,alpha,beta,i,j,c1,p_contrib",mpi%rank,x,chunk_size,alpha,beta,i,j,c1,p_contrib 
 
                  do imatrix = 3,4
                    call fac_values(sys_var%wcoul0,sys_var%vol,keyword,imatrix,fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16)
                    if (imatrix .eq.3) then
                       call load_Hhh(sys_var%nv,bse_mat_h)
                    end if
                    if (imatrix .eq. 3) then
                       call load_Hee(i,j,sys_var%nc,sys_var%nv,bse_mat_e)
                       call load_Heh(i,j,sys_var%nc,sys_var%nv,imatrix,bse_mat)
                    else if (imatrix .eq. 4) then
                       call load_Heh(i,j,sys_var%nc,sys_var%nv,imatrix,bse_mat)
                    end if
                    c2 =0
                    do n = 1,sys_var%nc-1
                      do m = n+1,sys_var%nc
                          do eta = 1,sys_var%nv-1
                             do gama = eta+1,sys_var%nv
                                 c2 = c2+1
                                 !print*,grid%myprow,c1,c2,i,j
                                 if (imatrix==3 .and. c1==c2) then
                                      H1(c2)=H1(c2)+ (sys_var%mf_c(j) + sys_var%mf_c(i) - sys_var%mf_v(beta) -sys_var%mf_v(alpha))
                                 end if

                                 if(alpha==gama .and. beta==eta) then
                                    if (imatrix .ne. 4) then
                                     !print*,"2",H(c1,c2),bse_mat(n,2,i,m)
                                       H1(c2)=H1(c2) + (fac_ed *bse_mat_e(n,1,1,m)) &
                                          + (fac_ex * bse_mat_e(m,1,1,n))
                                    end if
                                 end if
                                 if(i==m .and. j==n) then
                                   if (imatrix .ne. 4) then
                                      H1(c2) =H1(c2) + (fac_hd * bse_mat_h(eta,beta,alpha,gama)) &
                                       + (fac_hx * bse_mat_h(gama,beta,alpha,eta))   !<
                                       !print*, alpha,beta,gama,eta,bse_mat_h(eta,beta,alpha,gama)
                                   end if
                                 end if

                                 if(j==n .and. beta==gama) then
                                   if (imatrix .ne. 4) then
                                      H1(c2) =H1(c2) + fac_1 * bse_mat(alpha,eta,1,m)
                                   end if
                                   if (imatrix == 4) then
                                      H1(c2)=H1(c2) + fac_1 * bse_mat(alpha,eta,1,m)
                                      !print*,alpha,eta,i,m, bse_mat(alpha,eta,1,m)
                                   end if
                                 end if
                                 if(j==m .and. beta==gama) then
                                   if (imatrix .ne. 4) then
                                     H1(c2)=H1(c2) + fac_2 * bse_mat(alpha,eta,1,n)
                                   end if
                                   if (imatrix == 4) then
                                     H1(c2)=H1(c2) + fac_2 * bse_mat(alpha,eta,1,n)
                                   end if
                                end if


                                 if(i==n .and. beta==gama) then
                                   if (imatrix .ne. 4) then
                                      H1(c2)=H1(c2) + fac_3 * bse_mat(alpha,eta,2,m)
                                      !print*,alpha,eta,j,m, bse_mat(alpha,eta,2,m) 
                                      !print*,"5", i,j,alpha,beta,m,n,gama,eta
                                      !print*, H(c1,c2),c1,c2
                                   end if
                                   if (imatrix == 4) then
                                     H1(c2)=H1(c2) + fac_3 * bse_mat(alpha,eta,2,m)
                                     !print*,alpha,eta,j,m, bse_mat(alpha,eta,2,m)
                                   end if
                                 end if
                                 if(i==m .and. beta==gama) then
                                   if (imatrix .ne. 4) then
                                     H1(c2)=H1(c2) + fac_4 * bse_mat(alpha,eta,2,n)
                                    !print*,"6", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                                   end if
                                   if (imatrix == 4) then
                                     H1(c2)=H1(c2) + fac_4 * bse_mat(alpha,eta,2,n)
                                   end if
                                 end if
                                 if(j==n .and. alpha==eta) then
                                   if (imatrix .ne. 4) then
                                     H1(c2)=H1(c2) + fac_5 * bse_mat(beta,gama,1,m)
                                    !print*, H(c1,c2),c1,c2
                                   end if
                                   if (imatrix == 4) then
                                     H1(c2)=H1(c2) + fac_5 * bse_mat(beta,gama,1,m)
                                   end if
                                 end if
                                 if(j==m .and. alpha==eta) then
                                   if (imatrix .ne. 4) then
                                     H1(c2)=H1(c2) + fac_6 * bse_mat(beta,gama,1,n)
                                    !print*,"8", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                                   end if
                                   if (imatrix == 4) then
                                     H1(c2)=H1(c2) + fac_6 * bse_mat(beta,gama,1,n)
                                   end if
                                 end if
                                 if(i==n .and. alpha==eta) then
                                   if (imatrix .ne. 4) then
                                     H1(c2)=H1(c2) + fac_7 * bse_mat(beta,gama,2,m)
                                    !print*,"9", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                                   end if
                                   if (imatrix == 4) then
                                     H1(c2)=H1(c2) + fac_7 * bse_mat(beta,gama,2,m)
                                   end if
                                 end if

                                 if(i==m .and. alpha==eta) then
                                   if (imatrix .ne. 4) then
                                    H1(c2)=H1(c2) + fac_8 * bse_mat(beta,gama,2,n)
                                    !print*,"10", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                                   end if
                                   if (imatrix == 4) then
                                     H1(c2)=H1(c2) + fac_8 * bse_mat(beta,gama,2,n)
                                   end if
                                 end if
                                 if(j==n .and. beta==eta) then
                                   if (imatrix .ne. 4) then
                                    !print*,"11",H(c1,c2),bse_mat(alpha,gama,i,m)
                                     H1(c2)=H1(c2) +fac_9 * bse_mat(alpha,gama,1,m)
                                    !print*,"11", i,j,alpha,beta,m,n,gama,eta
                                    !print*,c1,c2,H(c1,c2)
                                   end if
                                   if (imatrix == 4) then
                                     H1(c2)=H1(c2) + fac_9 * bse_mat(alpha,gama,1,m)
                                   end if
                                  !print*,c1,c2,H(1,8),
                                 end if
                                 if(j==m .and. beta==eta) then
                                   if (imatrix .ne. 4) then
                                      H1(c2)=H1(c2) + fac_10 * bse_mat(alpha,gama,1,n)
                                     !print*,"12", i,j,alpha,beta,m,n,gama,eta
                                     !print*,c1,c2,H(c1,c2) 
                                   end if
                                   if (imatrix == 4) then
                                     H1(c2)=H1(c2) + fac_10 * bse_mat(alpha,gama,1,n)
                                   end if
                                 end if
                                  if(i==n .and. beta==eta) then
                                   if (imatrix .ne. 4) then
                                     H1(c2)=H1(c2) + fac_11 * bse_mat(alpha,gama,2,m)
                                    !print*,"13", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                                  end if
                                  if (imatrix == 4) then
                                     H1(c2) =H1(c2) + fac_11 * bse_mat(alpha,gama,2,m)
                                  end if
                                 end if
                                 if(i==m .and. beta==eta) then
                                   if (imatrix .ne. 4) then
                                    !print*,"14",H(c1,c2),bse_mat(alpha,gama,2,n)
                                    H1(c2) =H1(c2) + fac_12 * bse_mat(alpha,gama,2,n)
                                    !print*,"14", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                                   end if
                                   if (imatrix == 4) then
                                     H1(c2)=H1(c2) + fac_12 * bse_mat(alpha,gama,2,n)
                                   end if
                                  !print*,c1,c2,H(10,2),H(2,10)
                                 end if
                                 if(j==n .and. alpha==gama) then
                                   if (imatrix .ne. 4) then
                                     H1(c2)=H1(c2) + fac_13 * bse_mat(beta,eta,1,m)
                                    !print*,"15", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                                   end if
                                   if (imatrix == 4) then
                                     H1(c2)=H1(c2) + fac_13 * bse_mat(beta,eta,1,m)
                                   end if
                                 end if
                                 if(j==m .and. alpha==gama) then
                                   if (imatrix .ne. 4) then
                                     H1(c2)=H1(c2) + fac_14 * bse_mat(beta,eta,1,n)
                                    !print*,"16", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                                   end if
                                   if (imatrix == 4) then
                                     H1(c2)=H1(c2) + fac_14 * bse_mat(beta,eta,1,n)
                                   end if
                                 end if
                                 if(i==n .and. alpha==gama) then
                                   if (imatrix .ne. 4) then
                                     H1(c2)=H1(c2) + fac_15 * bse_mat(beta,eta,2,m)
                                    !print*,"17", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                                   end if
                                   if (imatrix == 4) then
                                     H1(c2)=H1(c2) + fac_15 * bse_mat(beta,eta,2,m)
                                   end if
                                 end if
                                 if(i==m .and. alpha==gama) then
                                   if (imatrix .ne. 4) then
                                    H1(c2)=H1(c2) + fac_16 * bse_mat(beta,eta,2,n)
                                    !print*,"18", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                                   end if
                                   if (imatrix == 4) then
                                     H1(c2)=H1(c2) + fac_16 * bse_mat(beta,eta,2,n)
                                  end if
                                 end if
                                 if (imatrix == 4) then
                                    !if (mpi%rank ==0 )then  
                                     ! print*,c1,c2,H1(c2) 
                                    !end if
                                     

                                 end if
                             end do
                          end do
                       end do
                    end do
                  end do
                 ! if (mpi%rank == 1)then
                  !print*,"mpi%rank,c1",mpi%rank,c1
                 ! end if
                end if
                 call distrubute_rows_2D_cyclo(H1,sys_var%n1,c1,1,x,chunk_size)
                 call MPI_Barrier(MPI_COMM_WORLD, error)


              end do
         end do
      end do
      

      call MPI_Barrier(MPI_COMM_WORLD, error)





    deallocate(H1)

  end subroutine create_H1  
 subroutine create_H2(nv,nc,mf_v,mf_c,wcoul0,s,vol,chunk)

      double precision :: fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16

      integer, intent(in) :: nv,nc,s,chunk

      integer :: nb,i,j,alpha,beta,gama,eta,m,n,imatrix,keyword,c1,c2,c,x,A,y
      double precision            :: element
      double precision ::wcoul0,vol,a1,b1,T1,T2,T3
      complex(kind=8),dimension(:,:),allocatable                            :: H2
      double precision, dimension(:,:),allocatable                          :: data_out
      complex(kind=8), allocatable :: bse_mat_e(:,:,:,:)
      complex(kind=8), allocatable :: bse_mat_h(:,:,:,:)
      complex(kind=8), allocatable :: bse_mat(:,:,:,:)
      CHARACTER(LEN=11), PARAMETER :: filename1 = "singlet.h5" ! File name
      CHARACTER(LEN=5), PARAMETER :: dsetname1 = "CISD2"    ! Dataset name
      CHARACTER(LEN=11), PARAMETER :: filename2 = "singletX.h5" ! File name
      CHARACTER(LEN=5), PARAMETER :: dsetname2 = "HCISD"    ! Dataset name

      INTEGER, PARAMETER   :: DRANK = 3 ! Dataset rank 1 dimension for real or complex 2 is row-index 3 is col-index
      INTEGER(SIZE_T), PARAMETER :: NUMP = 2 ! Number of points selected
      complex(kind=8) ,DIMENSION(nv) ,intent(in) :: mf_v
      complex(kind=8),DiMENSION(nc),intent(in) :: mf_c
      INTEGER(HID_T) :: file1_id       ! File1 identifier
      INTEGER(HID_T) :: dset1_id       ! Dataset1 identifier
      INTEGER(HID_T) :: file2_id       ! File2 identifier
      INTEGER(HID_T) :: dset2_id       ! Dataset2 identifier
      INTEGER(HID_T) :: dataspace1     ! Dataspace identifier
          ! memspace identifier
              
      INTEGER(HSIZE_T), DIMENSION(2) :: dimsm
                                                   ! Memory dataspace dimensions
      INTEGER(HSIZE_T), DIMENSION(2) :: dimsf
                                                   ! File dataspace dimensions
      INTEGER(HSIZE_T), DIMENSION(DRANK,NUMP) :: coord ! Elements coordinates
                                                      ! in the file
      DOUBLE PRECISION, DIMENSION(2) :: val =(/3.14,1.414/)  ! Values to write


      INTEGER :: memrank = 1  ! Rank of the dataset in memory
      DOUBLE PRECISION,DIMENSION(:,:,:),allocatable :: bufnew

      INTEGER :: error,info  ! Error flag
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims

      integer(HID_T) :: filespace     ! Dataspace identifier in file
      integer(HID_T) :: memspace      ! Dataspace identifier in mem
      integer(HID_T) :: dset_id       ! Dataset identifier
      integer(HSIZE_T) :: count(2), offset(2)
      INTEGER(HID_T) :: plist_id,fapl_id      ! Property list identifier

      !!!!!!!!!!!mpivariables!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer :: nproc, chunk_size,itr_length
      integer , dimension(:,:) , allocatable :: ij_mat
      integer, dimension(:),allocatable :: row_ij
      integer :: comm, rank, numproc, ierror
      integer :: cij,cab,cijp,c1p
      !call MPI_INIT(ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)



       keyword =2
       fac_hd=0
       fac_hx=0
       fac_ed=0
       fac_ex=0
      fac_1=0
      fac_2=0
      fac_3=0
      fac_4=0
      fac_5=0
      fac_6=0
      fac_7=0
      fac_8=0
      fac_9=0
      fac_10=0
      fac_11=0
      fac_12=0
      fac_13=0
      fac_14=0
      fac_15=0
      fac_16=0

      itr_length = nc*(nc-1)/2
      chunk_size = chunk
      print*,chunk_size,nproc
      !chunk_size = nc*(nc-1)/2


      allocate(bse_mat_h(nv,nv,nv,nv))
      allocate(bse_mat_e(nc,1,1,nc))
      allocate(bse_mat(nv,nv,2,nc))

      allocate(row_ij(chunk_size))
      A = nv*(nv-1)/2
      nrows_p= chunk_size*A
      allocate(ij_mat(2,chunk_size))
      allocate(H2(nrows_p,s))
      c = rank+1
      H2 = cmplx(0.0,0.0)
      call cpu_time(t1)
      !CALL h5screate_simple_f(1, dimsm, memspace, error)
       call distribute_ij(rank,nc,A,chunk_size,ij_mat,row_ij)
      c1=0
      c2=0
      cij=0
      cijp=0
      cab=0
      do imatrix = 3,4
       call fac_values(wcoul0,vol,keyword,imatrix,fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16)
         if (imatrix==3) then
            call load_Hhh(nv,bse_mat_h)
         end if
         c1 = 0

         do x = 1,chunk_size

      
            i = ij_mat(1,x)
            j = ij_mat(2,x)
            cijp = row_ij(x)
            cij =(x-1)*A
            cab =0

            if (imatrix .eq. 3) then
                call load_Hee(i,j,nc,nv,bse_mat_e)
                call load_Heh(i,j,nc,nv,imatrix,bse_mat)
            else if (imatrix .eq. 4) then
                call load_Heh(i,j,nc,nv,imatrix,bse_mat)
            end if
            do beta = 1,nv-1
               do alpha = beta+1,nv
                  cab= cab+1
                  c1 = cij+cab
                  c1p = cijp+cab
                  c2 = 0
                  print*,rank,c1,i,j,alpha,beta
                   do n = 1,nc-1
                      do m = n+1,nc
                          do eta = 1,nv-1
                             do gama = eta+1,nv
                                 c2 = c2+1
                                 !print*,rank,c1,c2
                                 if (imatrix==3 .and. c1p==c2) then
                                      H2(c1,c2)=H2(c1,c2)+ (mf_c(j) +mf_c(i) -mf_v(beta) -mf_v(alpha))
                                  end if

                                 if(alpha==gama .and. beta==eta) then
                                    if (imatrix .ne. 4) then
                                     !print*,"2",H(c1,c2),bse_mat(n,2,i,m)
                                       H2(c1,c2)=H2(c1,c2) + (fac_ed *bse_mat_e(n,1,1,m)) &
                                          + (fac_ex * bse_mat_e(m,1,1,n))
                                    end if
                                 end if
                                 if(i==m .and. j==n) then
                                   if (imatrix .ne. 4) then
                                      H2(c1,c2) =H2(c1,c2) + (fac_hd * bse_mat_h(eta,beta,alpha,gama)) &
                                       + (fac_hx * bse_mat_h(gama,beta,alpha,eta))   !<
                                       !print*, alpha,beta,gama,eta,bse_mat_h(eta,beta,alpha,gama)
                                   end if
                                 end if

                                 if(j==n .and. beta==gama) then
                                   if (imatrix .ne. 4) then
                                      H2(c1,c2) =H2(c1,c2) + fac_1 * bse_mat(alpha,eta,1,m)
                                   end if
                                   if (imatrix == 4) then
                                      H2(c1,c2)=H2(c1,c2) + fac_1 * bse_mat(alpha,eta,1,m)
                                      !print*,alpha,eta,i,m, bse_mat(alpha,eta,1,m)
                                   end if
                                 end if
                                 if(j==m .and. beta==gama) then
                                   if (imatrix .ne. 4) then
                                     H2(c1,c2)=H2(c1,c2) + fac_2 * bse_mat(alpha,eta,1,n)
                                   end if
                                   if (imatrix == 4) then
                                     H2(c1,c2)=H2(c1,c2) + fac_2 * bse_mat(alpha,eta,1,n)
                                   end if
                                end if
                                 if(i==n .and. beta==gama) then
                                   if (imatrix .ne. 4) then
                                      H2(c1,c2)=H2(c1,c2) + fac_3 * bse_mat(alpha,eta,2,m)
                                      !print*,alpha,eta,j,m, bse_mat(alpha,eta,2,m) 
                                      !print*,"5", i,j,alpha,beta,m,n,gama,eta
                                      !print*, H(c1,c2),c1,c2
                                   end if
                                   if (imatrix == 4) then
                                     H2(c1,c2)=H2(c1,c2) + fac_3 * bse_mat(alpha,eta,2,m)
                                     !print*,alpha,eta,j,m, bse_mat(alpha,eta,2,m)
                                   end if
                                 end if
                                 if(i==m .and. beta==gama) then
                                   if (imatrix .ne. 4) then
                                     H2(c1,c2)=H2(c1,c2) + fac_4 * bse_mat(alpha,eta,2,n)
                                    !print*,"6", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                                   end if
                                   if (imatrix == 4) then
                                     H2(c1,c2)=H2(c1,c2) + fac_4 * bse_mat(alpha,eta,2,n)
                                   end if
                                 end if
                                 if(j==n .and. alpha==eta) then
                                   if (imatrix .ne. 4) then
                                     H2(c1,c2)=H2(c1,c2) + fac_5 * bse_mat(beta,gama,1,m)
                                    !print*, H(c1,c2),c1,c2
                                   end if
                                   if (imatrix == 4) then
                                     H2(c1,c2)=H2(c1,c2) + fac_5 * bse_mat(beta,gama,1,m)
                                   end if
                                 end if
                                 if(j==m .and. alpha==eta) then
                                   if (imatrix .ne. 4) then
                                     H2(c1,c2)=H2(c1,c2) + fac_6 * bse_mat(beta,gama,1,n)
                                    !print*,"8", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                                   end if
                                   if (imatrix == 4) then
                                     H2(c1,c2)=H2(c1,c2) + fac_6 * bse_mat(beta,gama,1,n)
                                   end if
                                 end if
                                 if(i==n .and. alpha==eta) then
                                   if (imatrix .ne. 4) then
                                     H2(c1,c2)=H2(c1,c2) + fac_7 * bse_mat(beta,gama,2,m)
                                    !print*,"9", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                                   end if
                                   if (imatrix == 4) then
                                      H2(c1,c2)=H2(c1,c2) + fac_7 * bse_mat(beta,gama,2,m)
                                   end if
                                 end if

                                 if(i==m .and. alpha==eta) then
                                   if (imatrix .ne. 4) then
                                    H2(c1,c2)=H2(c1,c2) + fac_8 * bse_mat(beta,gama,2,n)
                                    !print*,"10", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                                   end if
                                   if (imatrix == 4) then
                                     H2(c1,c2)=H2(c1,c2) + fac_8 * bse_mat(beta,gama,2,n)
                                   end if
                                 end if
                                 if(j==n .and. beta==eta) then
                                   if (imatrix .ne. 4) then
                                    !print*,"11",H(c1,c2),bse_mat(alpha,gama,i,m)
                                     H2(c1,c2)=H2(c1,c2) +fac_9 * bse_mat(alpha,gama,1,m)
                                    !print*,"11", i,j,alpha,beta,m,n,gama,eta
                                    !print*,c1,c2,H(c1,c2)
                                   end if
                                   if (imatrix == 4) then
                                     H2(c1,c2)=H2(c1,c2) + fac_9 * bse_mat(alpha,gama,1,m)
                                   end if
                                  !print*,c1,c2,H(1,8),
                                 end if
                                 if(j==m .and. beta==eta) then
                                   if (imatrix .ne. 4) then
                                      H2(c1,c2)=H2(c1,c2) + fac_10 * bse_mat(alpha,gama,1,n)
                                     !print*,"12", i,j,alpha,beta,m,n,gama,eta
                                     !print*,c1,c2,H(c1,c2) 
                                   end if
                                   if (imatrix == 4) then
                                     H2(c1,c2)=H2(c1,c2) + fac_10 * bse_mat(alpha,gama,1,n)
                                   end if
                                 end if
                                  if(i==n .and. beta==eta) then
                                   if (imatrix .ne. 4) then
                                     H2(c1,c2)=H2(c1,c2) + fac_11 * bse_mat(alpha,gama,2,m)
                                    !print*,"13", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                                  end if
                                  if (imatrix == 4) then
                                     H2(c1,c2) =H2(c1,c2) + fac_11 * bse_mat(alpha,gama,2,m)
                                  end if
                                 end if
                                 if(i==m .and. beta==eta) then
                                   if (imatrix .ne. 4) then
                                    !print*,"14",H(c1,c2),bse_mat(alpha,gama,2,n)
                                    H2(c1,c2) =H2(c1,c2) + fac_12 * bse_mat(alpha,gama,2,n)
                                    !print*,"14", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                                   end if
                                   if (imatrix == 4) then
                                     H2(c1,c2)=H2(c1,c2) + fac_12 * bse_mat(alpha,gama,2,n)
                                   end if
                                  !print*,c1,c2,H(10,2),H(2,10)
                                 end if
                                 if(j==n .and. alpha==gama) then
                                   if (imatrix .ne. 4) then
                                     H2(c1,c2)=H2(c1,c2) + fac_13 * bse_mat(beta,eta,1,m)
                                    !print*,"15", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                                   end if
                                   if (imatrix == 4) then
                                     H2(c1,c2)=H2(c1,c2) + fac_13 * bse_mat(beta,eta,1,m)
                                   end if
                                 end if
                                 if(j==m .and. alpha==gama) then
                                   if (imatrix .ne. 4) then
                                     H2(c1,c2)=H2(c1,c2) + fac_14 * bse_mat(beta,eta,1,n)
                                    !print*,"16", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                                   end if
                                   if (imatrix == 4) then
                                     H2(c1,c2)=H2(c1,c2) + fac_14 * bse_mat(beta,eta,1,n)
                                   end if
                                 end if
                                 if(i==n .and. alpha==gama) then
                                   if (imatrix .ne. 4) then
                                     H2(c1,c2)=H2(c1,c2) + fac_15 * bse_mat(beta,eta,2,m)
                                    !print*,"17", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                                   end if
                                   if (imatrix == 4) then
                                     H2(c1,c2)=H2(c1,c2) + fac_15 * bse_mat(beta,eta,2,m)
                                   end if
                                 end if
                                 if(i==m .and. alpha==gama) then
                                   if (imatrix .ne. 4) then
                                    H2(c1,c2)=H2(c1,c2) + fac_16 * bse_mat(beta,eta,2,n)
                                    !print*,"18", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                                   end if
                                   if (imatrix == 4) then
                                     H2(c1,c2)=H2(c1,c2) + fac_16 * bse_mat(beta,eta,2,n)
                                  end if
                                 end if
                                 !print*,H2(c1,c2),c1,c2,imatrix

                             end do
                          end do
                      end do
                   end do
                end do
            end do
         end do
      end do
      call cpu_time(t2)
      CALL h5open_f(error)
      !CALL h5fcreate_f(filename1,  H5F_ACC_TRUNC_F, file1_id, error)


      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)

      !call h5fcreate_f(filename1, H5F_ACC_TRUNC_F, file1_id, error, &
      !            access_prp = plist_id)
      CALL h5fopen_f (filename1, H5F_ACC_RDWR_F, file1_id, error,access_prp = plist_id)


       !CALL h5fopen_f (filename2, H5F_ACC_RDWR_F, file2_id, error,access_prp = plist_id)
       CALL h5pclose_f(plist_id, error)
       dimsf(1) = s
       dimsf(2) = s
       CALL h5screate_simple_f(2, dimsf, dataspace1, error)

       CALL h5dcreate_f(file1_id, dsetname1,H5T_NATIVE_DOUBLE,dataspace1,dset1_id, error)
       CALL h5sclose_f(dataspace1, error)
     !
       !CALL h5dopen_f(file2_id, dsetname2, dset2_id, error)
     !
     ! Get  dataspace identifier.
     !
       count(1) = nrows_p
       count(2) = s
       dimsm(1)=nrows_p
       dimsm(2)=s
       dimsf(1) = s
       dimsf(2) = s
       !call h5screate_simple_f(2, dimsm, memspace, error)
       !allocate(data_out(count(1),count(2)))
       !data_out = real(H2)
       offset(1) = rank*nrows_p
       offset(2) = 0
      ! print*,"rank,count(1),count(2),offset(1),offset(2),s",rank,count(1),count(2),offset(1),offset(2),s
       !if ( rank == 1) then
        ! do x =1,6
         !   do y=1,30
          !     print*,H2(x,y),x,y
           ! end do
         !end do
       !end if 
       !CALL h5dget_space_f(dset1_id, dataspace1, error)

       !CALL h5dget_space_f(dset1_id, dataspace1, error)
       call h5screate_simple_f(2, count, memspace, error)
      ! CALL h5sselect_hyperslab_f (memspace, H5S_SELECT_SET_F, offset,&
       !        count, error)
       CALL h5dget_space_f(dset1_id, dataspace1, error)
       CALL h5sselect_hyperslab_f (dataspace1, H5S_SELECT_SET_F, offset,&
               count, error)


       CALL h5pcreate_f(H5P_DATASET_XFER_F, fapl_id, error)
       CALL h5pset_dxpl_mpio_f(fapl_id, H5FD_MPIO_COLLECTIVE_F, error)
       call h5screate_simple_f(2, dimsm, memspace, error)
       CALL h5dwrite_f(dset1_id, H5T_NATIVE_DOUBLE, real(H2), dimsf , error, &
                    mem_space_id=memspace, file_space_id=dataspace1)
       call h5pclose_f(fapl_id, error)

       CALL h5sclose_f(memspace, error)
       !call h5pclose_f(fapl_id, error)

       CALL h5dclose_f(dset1_id, error)

       CALL h5fclose_f(file1_id, error)

       CALL h5close_f(error)
       call cpu_time(t3)
                           
       print*,t2-t1,t3-t2


  end subroutine create_H2
    subroutine create_H12(nv,nc,mf_v,mf_c,wcoul0,s,vol,chunk)

      double precision :: fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16

      integer, intent(in) :: nv,nc,s,chunk

      integer :: nb,i,j,alpha,beta,gama,eta,m,n,imatrix,keyword,c1,c2,c,x,A,y
      double precision            :: element
      double precision ::wcoul0,vol,a1,b1,T1,T2,T3
      complex(kind=8),dimension(:,:),allocatable                            :: H12
      double precision, dimension(:,:),allocatable                          :: data_out
      complex(kind=8), allocatable :: bse_mat_e(:,:,:,:)
      complex(kind=8), allocatable :: bse_mat_h(:,:,:,:)
      complex(kind=8), allocatable :: bse_mat(:,:,:,:)
      CHARACTER(LEN=11), PARAMETER :: filename1 = "singlet.h5" ! File name
      CHARACTER(LEN=5), PARAMETER :: dsetname1 = "CIS12"    ! Dataset name
      CHARACTER(LEN=11), PARAMETER :: filename2 = "singletX.h5" ! File name
      CHARACTER(LEN=5), PARAMETER :: dsetname2 = "HCISD"    ! Dataset name

      INTEGER, PARAMETER   :: DRANK = 3 ! Dataset rank 1 dimension for real or complex 2 is row-index 3 is col-index
      INTEGER(SIZE_T), PARAMETER :: NUMP = 2 ! Number of points selected
      complex(kind=8) ,DIMENSION(nv) ,intent(in) :: mf_v
      complex(kind=8),DiMENSION(nc),intent(in) :: mf_c
      INTEGER(HID_T) :: file1_id       ! File1 identifier
      INTEGER(HID_T) :: dset1_id       ! Dataset1 identifier
      INTEGER(HID_T) :: file2_id       ! File2 identifier
      INTEGER(HID_T) :: dset2_id       ! Dataset2 identifier
      INTEGER(HID_T) :: dataspace1     ! Dataspace identifier
          ! memspace identifier

      INTEGER(HSIZE_T), DIMENSION(2) :: dimsm
                                                   ! Memory dataspace dimensions
      INTEGER(HSIZE_T), DIMENSION(2) :: dimsf
                                                   ! File dataspace dimensions
      INTEGER(HSIZE_T), DIMENSION(DRANK,NUMP) :: coord ! Elements coordinates
                                                      ! in the file
      DOUBLE PRECISION, DIMENSION(2) :: val =(/3.14,1.414/)  ! Values to write


      INTEGER :: memrank = 1  ! Rank of the dataset in memory
      DOUBLE PRECISION,DIMENSION(:,:,:),allocatable :: bufnew

      INTEGER :: error,info  ! Error flag
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims

      integer(HID_T) :: filespace     ! Dataspace identifier in file
      integer(HID_T) :: memspace      ! Dataspace identifier in mem
      integer(HID_T) :: dset_id       ! Dataset identifier
      integer(HSIZE_T) :: count(2), offset(2)
      INTEGER(HID_T) :: plist_id,fapl_id      ! Property list identifier

      !!!!!!!!!!!mpivariables!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer :: nproc, chunk_size,itr_length
      integer , dimension(:,:) , allocatable :: ij_mat
      integer, dimension(:),allocatable :: row_ij
      integer :: comm, rank, numproc, ierror
      integer :: cij,cab,cijp,c1p
      !call MPI_INIT(ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)



       keyword =3
       fac_hd=0
       fac_hx=0
       fac_ed=0
       fac_ex=0
      fac_1=0
      fac_2=0
      fac_3=0
      fac_4=0
      fac_5=0
      fac_6=0
      fac_7=0
      fac_8=0
      fac_9=0
      fac_10=0
      fac_11=0
      fac_12=0
      fac_13=0
      fac_14=0
      fac_15=0
      fac_16=0

      itr_length = nc*(nc-1)/2
      chunk_size = chunk
      print*,chunk_size,nproc
      !chunk_size = nc*(nc-1)/2


      allocate(bse_mat_h(nv,nv,nv,nv))
      allocate(bse_mat_e(nc,1,1,nc))
      allocate(bse_mat(nv,nv,2,nc))

      allocate(row_ij(chunk_size))
      A = nv*(nv-1)/2
      nrows_p= chunk_size*A
      allocate(ij_mat(2,chunk_size))
      allocate(H12(nrows_p,s))
      c = rank+1
      H12 = cmplx(0.0,0.0)
  
      call cpu_time(t1)
      !CALL h5screate_simple_f(1, dimsm, memspace, error)
       call distribute_ij(rank,nc,A,chunk_size,ij_mat,row_ij)
      c1=0
      c2=0
      cij=0
      cijp=0
      cab=0
      c1p = 0

      do imatrix = 3,4
       call fac_values(wcoul0,vol,keyword,imatrix,fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16)
         if (imatrix==3) then
            call load_Hhh(nv,bse_mat_h)
         end if
         c1 = 0

         do x = 1,chunk_size

            !print*,x,ij_mat(1,x),ij_mat(2,x)

            i = ij_mat(1,x)
            j = ij_mat(2,x)
            cijp = row_ij(x)
            cij =(x-1)*A
            cab =0

            if (imatrix .eq. 3) then
                call load_Hee(i,j,nc,nv,bse_mat_e)
                call load_Heh(i,j,nc,nv,imatrix,bse_mat)
            else if (imatrix .eq. 4) then
                call load_Heh(i,j,nc,nv,imatrix,bse_mat)
            end if
            do beta = 1,nv-1
               do alpha = beta+1,nv
                  cab= cab+1
                  c1 = cij+cab
                  c1p = cijp+cab
                  c2 = 0
                  print*,rank,c1,i,j,alpha,beta
                   do n = 1,nc-1
                      do m = n+1,nc
                          do eta = 1,nv-1
                             do gama = eta+1,nv
                                 c2 = c2+1
                                 !print*,rank,c1,c2

                                 if(alpha==gama .and. beta==eta) then
                                    if (imatrix .ne. 4) then
                                     !print*,"2",H(c1,c2),bse_mat(n,2,i,m)
                                       H12(c1,c2)=H12(c1,c2) + (fac_ed *bse_mat_e(n,1,1,m)) &
                                          + (fac_ex * bse_mat_e(m,1,1,n))
                                    end if
                                 end if
                                 if(i==m .and. j==n) then
                                   if (imatrix .ne. 4) then
                                      H12(c1,c2) =H12(c1,c2) + (fac_hd * bse_mat_h(eta,beta,alpha,gama)) &
                                       + (fac_hx * bse_mat_h(gama,beta,alpha,eta))   !<
                                       !print*, alpha,beta,gama,eta,bse_mat_h(eta,beta,alpha,gama)
                                   end if
                                 end if

                                 if(j==n .and. beta==gama) then
                                   if (imatrix .ne. 4) then
                                      H12(c1,c2) =H12(c1,c2) + fac_1 * bse_mat(alpha,eta,1,m)
                                   end if
                                   if (imatrix == 4) then
                                      H12(c1,c2)=H12(c1,c2) + fac_1 * bse_mat(alpha,eta,1,m)
                                      !print*,alpha,eta,i,m, bse_mat(alpha,eta,1,m)
                                   end if
                                 end if
                                 if(j==m .and. beta==gama) then
                                   if (imatrix .ne. 4) then
                                     H12(c1,c2)=H12(c1,c2) + fac_2 * bse_mat(alpha,eta,1,n)
                                   end if
                                   if (imatrix == 4) then
                                     H12(c1,c2)=H12(c1,c2) + fac_2 * bse_mat(alpha,eta,1,n)
                                   end if
                                end if
                                 if(i==n .and. beta==gama) then
                                   if (imatrix .ne. 4) then
                                      H12(c1,c2)=H12(c1,c2) + fac_3 * bse_mat(alpha,eta,2,m)
                                      !print*,alpha,eta,j,m, bse_mat(alpha,eta,2,m)
                                      !print*,"5", i,j,alpha,beta,m,n,gama,eta
                                      !print*, H(c1,c2),c1,c2
                                   end if
                                   if (imatrix == 4) then
                                     H12(c1,c2)=H12(c1,c2) + fac_3 * bse_mat(alpha,eta,2,m)
                                     !print*,alpha,eta,j,m, bse_mat(alpha,eta,2,m)
                                   end if
                                 end if
                                 if(i==m .and. beta==gama) then
                                   if (imatrix .ne. 4) then
                                     H12(c1,c2)=H12(c1,c2) + fac_4 * bse_mat(alpha,eta,2,n)
                                    !print*,"6", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                                   end if
                                   if (imatrix == 4) then
                                     H12(c1,c2)=H12(c1,c2) + fac_4 * bse_mat(alpha,eta,2,n)
                                   end if
                                 end if
                                 if(j==n .and. alpha==eta) then
                                   if (imatrix .ne. 4) then
                                     H12(c1,c2)=H12(c1,c2) + fac_5 * bse_mat(beta,gama,1,m)
                                    !print*, H(c1,c2),c1,c2
                                   end if
                                   if (imatrix == 4) then
                                     H12(c1,c2)=H12(c1,c2) + fac_5 * bse_mat(beta,gama,1,m)
                                   end if
                                 end if
                                 if(j==m .and. alpha==eta) then
                                   if (imatrix .ne. 4) then
                                     H12(c1,c2)=H12(c1,c2) + fac_6 * bse_mat(beta,gama,1,n)
                                    !print*,"8", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                                   end if
                                   if (imatrix == 4) then
                                     H12(c1,c2)=H12(c1,c2) + fac_6 * bse_mat(beta,gama,1,n)
                                   end if
                                 end if
                                 if(i==n .and. alpha==eta) then
                                   if (imatrix .ne. 4) then
                                     H12(c1,c2)=H12(c1,c2) + fac_7 * bse_mat(beta,gama,2,m)
                                    !print*,"9", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                                   end if
                                   if (imatrix == 4) then
                                      H12(c1,c2)=H12(c1,c2) + fac_7 * bse_mat(beta,gama,2,m)
                                   end if
                                 end if

                                 if(i==m .and. alpha==eta) then
                                   if (imatrix .ne. 4) then
                                    H12(c1,c2)=H12(c1,c2) + fac_8 * bse_mat(beta,gama,2,n)
                                    !print*,"10", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                                   end if
                                   if (imatrix == 4) then
                                     H12(c1,c2)=H12(c1,c2) + fac_8 * bse_mat(beta,gama,2,n)
                                   end if
                                 end if
                                 if(j==n .and. beta==eta) then
                                   if (imatrix .ne. 4) then
                                    !print*,"11",H(c1,c2),bse_mat(alpha,gama,i,m)
                                     H12(c1,c2)=H12(c1,c2) +fac_9 * bse_mat(alpha,gama,1,m)
                                    !print*,"11", i,j,alpha,beta,m,n,gama,eta
                                    !print*,c1,c2,H(c1,c2)
                                   end if
                                   if (imatrix == 4) then
                                     H12(c1,c2)=H12(c1,c2) + fac_9 * bse_mat(alpha,gama,1,m)
                                   end if
                                  !print*,c1,c2,H(1,8),
                                 end if
                                 if(j==m .and. beta==eta) then
                                   if (imatrix .ne. 4) then
                                      H12(c1,c2)=H12(c1,c2) + fac_10 * bse_mat(alpha,gama,1,n)
                                     !print*,"12", i,j,alpha,beta,m,n,gama,eta
                                     !print*,c1,c2,H(c1,c2)
                                   end if
                                   if (imatrix == 4) then
                                     H12(c1,c2)=H12(c1,c2) + fac_10 * bse_mat(alpha,gama,1,n)
                                   end if
                                 end if
                                  if(i==n .and. beta==eta) then
                                   if (imatrix .ne. 4) then
                                     H12(c1,c2)=H12(c1,c2) + fac_11 * bse_mat(alpha,gama,2,m)
                                    !print*,"13", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                                  end if
                                  if (imatrix == 4) then
                                     H12(c1,c2) =H12(c1,c2) + fac_11 * bse_mat(alpha,gama,2,m)
                                  end if
                                 end if
                                 if(i==m .and. beta==eta) then
                                   if (imatrix .ne. 4) then
                                    !print*,"14",H(c1,c2),bse_mat(alpha,gama,2,n)
                                    H12(c1,c2) =H12(c1,c2) + fac_12 * bse_mat(alpha,gama,2,n)
                                    !print*,"14", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                                   end if
                                   if (imatrix == 4) then
                                     H12(c1,c2)=H12(c1,c2) + fac_12 * bse_mat(alpha,gama,2,n)
                                   end if
                                  !print*,c1,c2,H(10,2),H(2,10)
                                 end if
                                 if(j==n .and. alpha==gama) then
                                   if (imatrix .ne. 4) then
                                     H12(c1,c2)=H12(c1,c2) + fac_13 * bse_mat(beta,eta,1,m)
                                    !print*,"15", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                                   end if
                                   if (imatrix == 4) then
                                     H12(c1,c2)=H12(c1,c2) + fac_13 * bse_mat(beta,eta,1,m)
                                   end if
                                 end if
                                 if(j==m .and. alpha==gama) then
                                   if (imatrix .ne. 4) then
                                     H12(c1,c2)=H12(c1,c2) + fac_14 * bse_mat(beta,eta,1,n)
                                    !print*,"16", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                                   end if
                                   if (imatrix == 4) then
                                     H12(c1,c2)=H12(c1,c2) + fac_14 * bse_mat(beta,eta,1,n)
                                   end if
                                 end if
                                 if(i==n .and. alpha==gama) then
                                   if (imatrix .ne. 4) then
                                     H12(c1,c2)=H12(c1,c2) + fac_15 * bse_mat(beta,eta,2,m)
                                    !print*,"17", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                                   end if
                                   if (imatrix == 4) then
                                     H12(c1,c2)=H12(c1,c2) + fac_15 * bse_mat(beta,eta,2,m)
                                   end if
                                 end if
                                 if(i==m .and. alpha==gama) then
                                   if (imatrix .ne. 4) then
                                    H12(c1,c2)=H12(c1,c2) + fac_16 * bse_mat(beta,eta,2,n)
                                    !print*,"18", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                                   end if
                                   if (imatrix == 4) then
                                     H12(c1,c2)=H12(c1,c2) + fac_16 * bse_mat(beta,eta,2,n)
                                  end if
                                 end if
                                 !print*,H12(c1,c2),c1,c2,imatrix

                             end do
                          end do
                      end do
                   end do
                end do
            end do
         end do
      end do
      call cpu_time(t2)
      CALL h5open_f(error)
      !CALL h5fcreate_f(filename1,  H5F_ACC_TRUNC_F, file1_id, error)


      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)

      !call h5fcreate_f(filename1, H5F_ACC_TRUNC_F, file1_id, error, &
       !           access_prp = plist_id)
      CALL h5fopen_f (filename1, H5F_ACC_RDWR_F, file1_id, error,access_prp = plist_id)


       !CALL h5fopen_f (filename2, H5F_ACC_RDWR_F, file2_id, error,access_prp = plist_id)
       CALL h5pclose_f(plist_id, error)
       dimsf(1) = s
       dimsf(2) = s
       CALL h5screate_simple_f(2, dimsf, dataspace1, error)

       CALL h5dcreate_f(file1_id, dsetname1,H5T_NATIVE_DOUBLE,dataspace1,dset1_id, error)
       CALL h5sclose_f(dataspace1, error)
     !
       !CALL h5dopen_f(file2_id, dsetname2, dset2_id, error)
     !
     ! Get  dataspace identifier.
     !
      count(1) = nrows_p
       count(2) = s
       dimsm(1)=nrows_p
       dimsm(2)=s
       dimsf(1) = s
       dimsf(2) = s
       !call h5screate_simple_f(2, dimsm, memspace, error)
       !allocate(data_out(count(1),count(2)))
       !data_out = real(H12)
       offset(1) = rank*nrows_p
       offset(2) = 0
       call h5screate_simple_f(2, count, memspace, error)
       CALL h5dget_space_f(dset1_id, dataspace1, error)
       CALL h5sselect_hyperslab_f (dataspace1, H5S_SELECT_SET_F, offset,&
               count, error)


       CALL h5pcreate_f(H5P_DATASET_XFER_F, fapl_id, error)
       CALL h5pset_dxpl_mpio_f(fapl_id, H5FD_MPIO_COLLECTIVE_F, error)
       call h5screate_simple_f(2, dimsm, memspace, error)
       CALL h5dwrite_f(dset1_id, H5T_NATIVE_DOUBLE, real(H12), dimsf , error, &
                    mem_space_id=memspace, file_space_id=dataspace1)
       call h5pclose_f(fapl_id, error)

       CALL h5sclose_f(memspace, error)
       !call h5pclose_f(fapl_id, error)

       CALL h5dclose_f(dset1_id, error)

       CALL h5fclose_f(file1_id, error)

       CALL h5close_f(error)
       call cpu_time(t3)

       print*,t2-t1,t3-t2


  end subroutine create_H12
      subroutine create_H3(nv,nc,mf_v,mf_c,wcoul0,s,vol)
      double precision :: fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16

      integer, intent(in) :: nv,nc,s
      integer :: nb,i,j,alpha,beta,gama,eta,m,n,imatrix,keyword,c1,c2,c,x,A
      double precision            :: element
      double precision ::wcoul0,vol,a1,b1
      complex(kind=8),dimension(:,:),allocatable                            :: H3
      complex(kind=8), allocatable :: bse_mat_e(:,:,:,:)
      complex(kind=8), allocatable :: bse_mat_h(:,:,:,:)
      complex(kind=8), allocatable :: bse_mat(:,:,:,:)
      CHARACTER(LEN=11), PARAMETER :: filename1 = "singlet.h5" ! File name
      CHARACTER(LEN=5), PARAMETER :: dsetname1 = "CISD3"    ! Dataset name
      CHARACTER(LEN=11), PARAMETER :: filename2 = "singletX.h5" ! File name
      CHARACTER(LEN=5), PARAMETER :: dsetname2 = "HCISD"    ! Dataset name
      INTEGER, PARAMETER   :: DRANK = 3 ! Dataset rank 1 dimension for real or complex 2 is row-index 3 is col-index
      INTEGER(SIZE_T), PARAMETER :: NUMP = 2 ! Number of points selected
      complex(kind=8) ,DIMENSION(nv) ,intent(in) :: mf_v
      complex(kind=8),DiMENSION(nc),intent(in) :: mf_c
      INTEGER(HID_T) :: file1_id       ! File1 identifier
      INTEGER(HID_T) :: dset1_id       ! Dataset1 identifier
      INTEGER(HID_T) :: file2_id       ! File2 identifier
      INTEGER(HID_T) :: dset2_id       ! Dataset2 identifier
      INTEGER(HID_T) :: dataspace1     ! Dataspace identifier
          ! memspace identifier
      INTEGER(HID_T) :: dataspace2     ! Dataspace identifier
            ! memspace identifier
      INTEGER(HID_T) :: fapl_id

      INTEGER(HSIZE_T), DIMENSION(2) :: dimsm
                                                   ! Memory dataspace dimensions
      INTEGER(HSIZE_T), DIMENSION(2) :: dimsf
      integer(HSIZE_T) :: count(2), offset(2)
 
     

      INTEGER :: error,info  ! Error flag
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims

      integer(HID_T) :: filespace     ! Dataspace identifier in file
      integer(HID_T) :: memspace      ! Dataspace identifier in mem
      integer(HID_T) :: dset_id       ! Dataset identifier
      INTEGER(HID_T) :: plist_id     ! Property list identifier
      !!!!!!!!!!!mpivariables!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer :: nproc, chunk_size,itr_length
      integer , dimension(:,:) , allocatable :: ij_mat
      integer, dimension(:),allocatable :: row_ij
      integer :: comm, rank, numproc, ierror
      integer :: cij,cab,cijp,c1p
      !call MPI_INIT(ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)

     
       keyword = 4

       fac_hd=0
       fac_hx=0
       fac_ed=0
       fac_ex=0
      fac_1=0
      fac_2=0
      fac_3=0
      fac_4=0
      fac_5=0
      fac_6=0
      fac_7=0
      fac_8=0
      fac_9=0
      fac_10=0
      fac_11=0
      fac_12=0
      fac_13=0
      fac_14=0
      fac_15=0
      fac_16=0

       itr_length = nc*(nc-1)/2
      if (mod(nc,2)==0) then
         nproc = nc-1
         chunk_size = nc/2
      else
         nproc = nc
         chunk_size = (nc-1)/2
      end if
      print*,chunk_size,nproc

      allocate(bse_mat_h(nv,nv,nv,nv))
      allocate(bse_mat_e(nc,1,1,nc))
      allocate(bse_mat(nv,nv,2,nc))
      !allocate(row_ij(chunk_size))
      A = nv
      !allocate(ij_mat(2,chunk_size))
      c = rank+1
      nrows_p= A
      allocate(H3(nrows_p,nc*nv))

      H3 = cmplx(0.0,0.0)
      

      c1=0
      c2=0
      cij=0
      cijp=0
      cijp=0
      c1p = 0
      cab=0
      do imatrix = 3,4
         call fac_values(wcoul0,vol,keyword,imatrix,fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16)

         i = c
         cij = 0
         cijp = (c-1)*A
         cab = 0
         if (imatrix==3) then
            call load_Hhh(nv,bse_mat_h)
         end if

         if (imatrix .eq. 3) then
                call load_Hee(i,i,nc,nv,bse_mat_e)
                call load_Heh(i,i,nc,nv,imatrix,bse_mat)
         else if (imatrix .eq. 4) then
                call load_Heh(i,i,nc,nv,imatrix,bse_mat)
         end if

         do alpha = 1,nv
            cab = cab + 1
            c1 = cij + cab
            c1p = cijp + cab
            c2 = 0
            do m = 1,nc
               do gama = 1,nv
                  c2 = c2 + 1
                  !print*,i,alpha,m,gama,c1,c2,c1p
                  if (imatrix==3 .and. c1p==c2) then
                      H3(c1,c2)=H3(c1,c2)+ (2*mf_c(i) - 2*mf_v(alpha))
                  end if

                  if(i==m) then
                     if (imatrix .ne. 4) then
                                     !print*,"bb",c1,c2, H(5,11),H(71,1),imatrix
                                      !print*,"1",H(c1,c2),bse_mat(eta,beta,alpha,gama), bse_mat(gama,beta,alpha,eta)
                         H3(c1,c2) = H3(c1,c2) + (fac_hd * bse_mat_h(gama,alpha,alpha,gama)) &
                                       + (fac_hx * bse_mat_h(gama,alpha,alpha,gama))   !<
                                     !print*,"1", i,j,alpha,beta,m,n,gama,eta
                                     !print*,c1,c2,H(c1,c2) 
                     end if
                  end if
                  if(alpha==gama) then
                     if (imatrix .ne. 4) then
                                     !print*,"2",H(c1,c2),bse_mat(n,j+nv,1,m)
                        H3(c1,c2)=H3(c1,c2) + (fac_ed *bse_mat_e(m,1,1,m)) &
                                       + (fac_ex * bse_mat_e(m,1,1,m))
                                    !print*,"2", i,j,alpha,beta,m,n,gama,eta
                                    !print*,c1,c2,H(c1,c2)
                     end if
                  end if

                  if(i==m .and. alpha==gama) then
                     if (imatrix .ne. 4) then
                         H3(c1,c2) = H3(c1,c2) + fac_1 * bse_mat(alpha,gama,1,m)
                         H3(c1,c2) = H3(c1,c2) + fac_2 * bse_mat(alpha,gama,1,m)
                         H3(c1,c2) = H3(c1,c2) + fac_3 * bse_mat(alpha,gama,1,m)
                         H3(c1,c2) = H3(c1,c2) + fac_4 * bse_mat(alpha,gama,1,m)
                         H3(c1,c2) = H3(c1,c2) +  fac_5 * bse_mat(alpha,gama,1,m)
                         H3(c1,c2) = H3(c1,c2) + fac_6 * bse_mat(alpha,gama,1,m)
                         H3(c1,c2) = H3(c1,c2) + fac_7 * bse_mat(alpha,gama,1,m)
                         H3(c1,c2) = H3(c1,c2) + fac_8 * bse_mat(alpha,gama,1,m)
                         H3(c1,c2) = H3(c1,c2) + fac_9 * bse_mat(alpha,gama,1,m)
                         H3(c1,c2) = H3(c1,c2) + fac_10 * bse_mat(alpha,gama,1,m)
                         H3(c1,c2) = H3(c1,c2) + fac_11 * bse_mat(alpha,gama,1,m)
                         H3(c1,c2) = H3(c1,c2) + fac_12* bse_mat(alpha,gama,1,m)
                         H3(c1,c2) = H3(c1,c2) + fac_13 * bse_mat(alpha,gama,1,m)
                         H3(c1,c2) = H3(c1,c2) + fac_14 * bse_mat(alpha,gama,1,m)
                         H3(c1,c2) = H3(c1,c2) + fac_15 * bse_mat(alpha,gama,1,m)
                         H3(c1,c2) = H3(c1,c2) + fac_16 * bse_mat(alpha,gama,1,m)

                       end if
                       if (imatrix == 4) then
                         H3(c1,c2) = H3(c1,c2) + fac_1 * bse_mat(alpha,gama,1,m)
                         H3(c1,c2) = H3(c1,c2) + fac_2 * bse_mat(alpha,gama,1,m)
                         H3(c1,c2) = H3(c1,c2) + fac_3 * bse_mat(alpha,gama,1,m)
                         H3(c1,c2) = H3(c1,c2) + fac_4 * bse_mat(alpha,gama,1,m)
                         H3(c1,c2) = H3(c1,c2) + fac_5 * bse_mat(alpha,gama,1,m)
                         H3(c1,c2) = H3(c1,c2) + fac_6 * bse_mat(alpha,gama,1,m)
                         H3(c1,c2) = H3(c1,c2) + fac_7 * bse_mat(alpha,gama,1,m)
                         H3(c1,c2) = H3(c1,c2) + fac_8 * bse_mat(alpha,gama,1,m)
                         H3(c1,c2) = H3(c1,c2) + fac_9 * bse_mat(alpha,gama,1,m)
                         H3(c1,c2) = H3(c1,c2) + fac_10 * bse_mat(alpha,gama,1,m)
                         H3(c1,c2) = H3(c1,c2) + fac_11 * bse_mat(alpha,gama,1,m)
                         H3(c1,c2) = H3(c1,c2) + fac_12 * bse_mat(alpha,gama,1,m)
                         H3(c1,c2) = H3(c1,c2) + fac_13 * bse_mat(alpha,gama,1,m)
                         H3(c1,c2) = H3(c1,c2) + fac_14 * bse_mat(alpha,gama,1,m)
                         H3(c1,c2) = H3(c1,c2) + fac_15 * bse_mat(alpha,gama,1,m)
                         H3(c1,c2) = H3(c1,c2) + fac_16 * bse_mat(alpha,gama,1,m)
                        
                     end if
                  end if
                  !print*,H3(c1,c2),c1p,c2
               end do
            end do
          end do
       end do
       call cpu_time(t2)
      CALL h5open_f(error)
      !CALL h5fcreate_f(filename1,  H5F_ACC_TRUNC_F, file1_id, error)


      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)

      !call h5fcreate_f(filename1, H5F_ACC_TRUNC_F, file1_id, error, &
       !           access_prp = plist_id)
      CALL h5fopen_f (filename1, H5F_ACC_RDWR_F, file1_id, error,access_prp = plist_id)


       !CALL h5fopen_f (filename2, H5F_ACC_RDWR_F, file2_id, error,access_prp = plist_id)
       CALL h5pclose_f(plist_id, error)
       dimsf(1) = nc*nv
       dimsf(2) = nc*nv
       CALL h5screate_simple_f(2, dimsf, dataspace1, error)

       CALL h5dcreate_f(file1_id, dsetname1,H5T_NATIVE_DOUBLE,dataspace1,dset1_id, error)
       CALL h5sclose_f(dataspace1, error)
     !
       !CALL h5dopen_f(file2_id, dsetname2, dset2_id, error)
     !
     ! Get  dataspace identifier.
     !
      count(1) = nrows_p
       count(2) = nc*nv
       dimsm(1)=nrows_p
       dimsm(2)=nc*nv
       dimsf(1) = nc*nv
       dimsf(2) = nc*nv
       !print*,"rank,count(1),count(2),offset(1),offset(2),nc*nv",rank,count(1),count(2),offset(1),offset(2),nc*nv

       !call h5screate_simple_f(2, dimsm, memspace, error)
       !allocate(data_out(count(1),count(2)))
       !data_out = real(H12)
       offset(1) = rank*nrows_p
       offset(2) = 0
       print*,"rank,count(1),count(2),offset(1),offset(2),nc*nv",rank,count(1),count(2),offset(1),offset(2),nc*nv

       call h5screate_simple_f(2, count, memspace, error)
       CALL h5dget_space_f(dset1_id, dataspace1, error)
       CALL h5sselect_hyperslab_f (dataspace1, H5S_SELECT_SET_F, offset,&
               count, error)


       CALL h5pcreate_f(H5P_DATASET_XFER_F, fapl_id, error)
       CALL h5pset_dxpl_mpio_f(fapl_id, H5FD_MPIO_COLLECTIVE_F, error)
       call h5screate_simple_f(2, dimsm, memspace, error)
       CALL h5dwrite_f(dset1_id, H5T_NATIVE_DOUBLE, real(H3), dimsf , error, &
                    mem_space_id=memspace, file_space_id=dataspace1)
       call h5pclose_f(fapl_id, error)

       CALL h5sclose_f(memspace, error)
       CALL h5dclose_f(dset1_id, error)

       CALL h5fclose_f(file1_id, error)

       CALL h5close_f(error)
       call cpu_time(t3)

       print*,t2-t1,t3-t2

                                    
       
  end subroutine create_H3
    subroutine create_H13(nv,nc,mf_v,mf_c,wcoul0,s,vol,chunk)
      double precision :: fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16

      integer, intent(in) :: nv,nc,s,chunk
      integer :: nb,i,j,alpha,beta,gama,eta,m,n,imatrix,keyword,c1,c2,c,x,A
      double precision            :: element
      double precision ::wcoul0,vol,a1,b1
      complex(kind=8),dimension(:,:),allocatable                           :: H13
      complex(kind=8), allocatable :: bse_mat_e(:,:,:,:)
      complex(kind=8), allocatable :: bse_mat_h(:,:,:,:)
      complex(kind=8), allocatable :: bse_mat(:,:,:,:)
      CHARACTER(LEN=11), PARAMETER :: filename1 = "singlet.h5" ! File name
      CHARACTER(LEN=5), PARAMETER :: dsetname1 = "CIS13"    ! Dataset name
      CHARACTER(LEN=11), PARAMETER :: filename2 = "singletX.h5" ! File name
      CHARACTER(LEN=5), PARAMETER :: dsetname2 = "HCISD"    ! Dataset name

      INTEGER, PARAMETER   :: DRANK = 3 ! Dataset rank 1 dimension for real or complex 2 is row-index 3 is col-index
      INTEGER(SIZE_T), PARAMETER :: NUMP = 2 ! Number of points selected
      complex(kind=8) ,DIMENSION(nv) ,intent(in) :: mf_v
      complex(kind=8),DiMENSION(nc),intent(in) :: mf_c
      INTEGER(HID_T) :: file1_id       ! File1 identifier
      INTEGER(HID_T) :: dset1_id       ! Dataset1 identifier
      INTEGER(HID_T) :: file2_id       ! File2 identifier
      INTEGER(HID_T) :: dset2_id       ! Dataset2 identifier
      INTEGER(HID_T) :: dataspace1     ! Dataspace identifier
          ! memspace identifier
      INTEGER(HID_T) :: dataspace2     ! Dataspace identifier
            ! memspace identifier
      INTEGER(HID_T) :: fapl_id

      INTEGER(HSIZE_T), DIMENSION(2) :: dimsm
                                                   ! Memory dataspace dimensions
      INTEGER(HSIZE_T), DIMENSION(2) :: dimsf
                                                   ! File dataspace dimensions
      

      INTEGER :: error,info  ! Error flag
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims
      integer(HSIZE_T) :: count(2), offset(2)

      integer(HID_T) :: filespace     ! Dataspace identifier in file
      integer(HID_T) :: memspace      ! Dataspace identifier in mem
      integer(HID_T) :: dset_id       ! Dataset identifier
      INTEGER(HID_T) :: plist_id      ! Property list identifier
      !!!!!!!!!!!mpivariables!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer :: nproc, chunk_size,itr_length
      integer , dimension(:,:) , allocatable :: ij_mat
      integer, dimension(:),allocatable :: row_ij
      integer :: comm, rank, numproc, ierror
      integer :: cij,cab,cijp,c1p
      !call MPI_INIT(ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)

     
       keyword = 5

       fac_hd=0
       fac_hx=0
       fac_ed=0
       fac_ex=0
      fac_1=0
      fac_2=0
      fac_3=0
      fac_4=0
      fac_5=0
      fac_6=0
      fac_7=0
      fac_8=0
      fac_9=0
      fac_10=0
      fac_11=0
      fac_12=0
      fac_13=0
      fac_14=0
      fac_15=0
      fac_16=0

       itr_length = nc*(nc-1)/2
       chunk_size = chunk
      print*,chunk_size,nproc

      allocate(bse_mat_h(nv,nv,nv,nv))
      allocate(bse_mat_e(nc,1,1,nc))
      allocate(bse_mat(nv,nv,2,nc))

      allocate(row_ij(chunk_size))
      A = (nv*(nv-1))/2
      allocate(ij_mat(2,chunk_size))
      c = rank+1
      nrows_p= chunk_size*A
      allocate(H13(nrows_p,nc*nv))
      H13 = cmplx(0.0,0.0)
      


      !print*,"proc_id+1",c
      call distribute_ij(rank,nc,A,chunk_size,ij_mat,row_ij)
      c1=0
      c2=0
      cij=0
      cab=0
      do imatrix = 3,4
       call fac_values(wcoul0,vol,keyword,imatrix,fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16)
         if (imatrix==3) then
            call load_Hhh(nv,bse_mat_h)
         end if
         do x = 1,chunk_size
            print*,x,ij_mat(1,x),ij_mat(2,x)

            i = ij_mat(1,x)
            j = ij_mat(2,x)
            cij = (x-1)*A

            cab =0

            if (imatrix .eq. 3) then
                call load_Hee(i,j,nc,nv,bse_mat_e)
                call load_Heh(i,j,nc,nv,imatrix,bse_mat)
            else if (imatrix .eq. 4) then
                call load_Heh(i,j,nc,nv,imatrix,bse_mat)
            end if
            do beta = 1,nv-1
               do alpha = beta+1,nv
                  cab = cab +1
                  c1 = cij+cab
                  c2 = 0
                  do m = 1,nc
                     do gama = 1,nv
                        c2 = c2 +1
                        if(j==m .and. beta==gama) then
                           if (imatrix .ne. 4) then
                               H13(c1,c2) = H13(c1,c2) + fac_1 * bse_mat(alpha,gama,1,m)
                               H13(c1,c2) = H13(c1,c2) + fac_2 * bse_mat(alpha,gama,1,m)
                               H13(c1,c2) = H13(c1,c2) + fac_9 * bse_mat(alpha,gama,1,m)
                               H13(c1,c2) = H13(c1,c2) + fac_10 * bse_mat(alpha,gama,1,m)
                               !print*,"3", i,j,alpha,beta,m,n,gama,c1,c2,H13(c1,c2)
                                    !print*, H(c1,c2),c1,c2
                           end if
                           if (imatrix == 4) then
                               H13(c1,c2) = H13(c1,c2) + fac_1 * bse_mat(alpha,gama,1,m)
                               H13(c1,c2) = H13(c1,c2) + fac_2 * bse_mat(alpha,gama,1,m)
                               H13(c1,c2) = H13(c1,c2) + fac_9 * bse_mat(alpha,gama,1,m)
                               H13(c1,c2) = H13(c1,c2) + fac_10 * bse_mat(alpha,gama,1,m)
                               !print*,"3", i,j,alpha,beta,m,n,gama,c1,c2,H13(c1,c2)
                           end if
                        end if

                        if(i==m .and. beta==gama) then
                           if (imatrix .ne. 4) then
                               H13(c1,c2)=H13(c1,c2) + fac_3 * bse_mat(alpha,gama,2,m)
                               H13(c1,c2)=H13(c1,c2) + fac_4 * bse_mat(alpha,gama,2,m)
                               H13(c1,c2)=H13(c1,c2) + fac_11 * bse_mat(alpha,gama,2,m)
                               H13(c1,c2)=H13(c1,c2) + fac_12 * bse_mat(alpha,gama,2,m)
                               !print*,"3", i,j,alpha,beta,m,n,gama,c1,c2,H13(c1,c2)
                                      !print*,"5", i,j,alpha,beta,m,n,gama,eta
                                      !print*, H(c1,c2),c1,c2
                           end if
                           if (imatrix == 4) then
                               H13(c1,c2)=H13(c1,c2) + fac_3 * bse_mat(alpha,gama,2,m)
                               H13(c1,c2)=H13(c1,c2) + fac_4 * bse_mat(alpha,gama,2,m)
                               H13(c1,c2)=H13(c1,c2) + fac_11 * bse_mat(alpha,gama,2,m)
                               H13(c1,c2)=H13(c1,c2) + fac_12 * bse_mat(alpha,gama,2,m)
                               !print*,"3","hey", i,j,alpha,beta,m,n,gama,c1,c2,H13(c1,c2)
                           end if
                        end if

                        if(j==m .and. alpha==gama) then
                           if (imatrix .ne. 4) then
                               !print*,"3", i,j,alpha,beta,m,n,gama,c1,c2,H13(c1,c2),bse_mat(beta,gama,1,m),&
                                                        !fac_5 * bse_mat(beta,gama,1,m),fac_6 * bse_mat(beta,gama,1,m)
                               H13(c1,c2)=H13(c1,c2) + fac_5 * bse_mat(beta,gama,1,m)
                               H13(c1,c2)=H13(c1,c2) + fac_6 * bse_mat(beta,gama,1,m)
                               H13(c1,c2)=H13(c1,c2) + fac_13 * bse_mat(beta,gama,1,m)
                               H13(c1,c2)=H13(c1,c2) + fac_14 * bse_mat(beta,gama,1,m)
                                !print*,"3", i,j,alpha,beta,m,n,gama,c1,c2,H13(c1,c2),bse_mat(beta,gama,1,m),&
                                                        !fac_5 * bse_mat(beta,gama,1,m),fac_6 * bse_mat(beta,gama,1,m)
                               !print*, fac_5,fac_6
                           end if
                           if (imatrix == 4) then
                               H13(c1,c2)=H13(c1,c2) + fac_5 * bse_mat(beta,gama,1,m)
                               H13(c1,c2)=H13(c1,c2) + fac_6 * bse_mat(beta,gama,1,m)
                               H13(c1,c2)=H13(c1,c2) + fac_13 * bse_mat(beta,gama,1,m)
                               H13(c1,c2)=H13(c1,c2) + fac_14 * bse_mat(beta,gama,1,m)
                               !print*,"3", i,j,alpha,beta,m,n,gama,c1,c2,H13(c1,c2)
                           end if
                        end if

                        if(i==m .and. alpha==gama) then
                           if (imatrix .ne. 4) then
                               H13(c1,c2)=H13(c1,c2) + fac_7 * bse_mat(beta,gama,2,m)
                               H13(c1,c2)=H13(c1,c2) + fac_8 * bse_mat(beta,gama,2,m)
                               H13(c1,c2)=H13(c1,c2) + fac_15 * bse_mat(beta,gama,2,m)
                               H13(c1,c2)=H13(c1,c2) + fac_16 * bse_mat(beta,gama,2,m)
                                    !print*,"9", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                               !print*,"3", i,j,alpha,beta,m,n,gama,c1,c2,H13(c1,c2),bse_mat(beta,gama,1,m)

                            end if
                            if (imatrix == 4) then
                                H13(c1,c2)=H13(c1,c2) + fac_7 * bse_mat(beta,gama,2,m)
                                H13(c1,c2)=H13(c1,c2) + fac_8 * bse_mat(beta,gama,2,m)
                                H13(c1,c2)=H13(c1,c2) + fac_15 * bse_mat(beta,gama,2,m)
                                H13(c1,c2)=H13(c1,c2) + fac_16 * bse_mat(beta,gama,2,m)
                                !print*,"3", i,j,alpha,beta,m,n,gama,c1,c2,H13(c1,c2),bse_mat(beta,gama,1,m)
                            end if
                        end if
                       
                     end do
                  end do
               end do
            end do
         end do
      end do
       call cpu_time(t2)
      CALL h5open_f(error)
      !CALL h5fcreate_f(filename1,  H5F_ACC_TRUNC_F, file1_id, error)


      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)

      !call h5fcreate_f(filename1, H5F_ACC_TRUNC_F, file1_id, error, &
       !           access_prp = plist_id)
      CALL h5fopen_f (filename1, H5F_ACC_RDWR_F, file1_id, error,access_prp = plist_id)


       !CALL h5fopen_f (filename2, H5F_ACC_RDWR_F, file2_id, error,access_prp = plist_id)
       CALL h5pclose_f(plist_id, error)
       dimsf(1) = s
       dimsf(2) = nc*nv
       CALL h5screate_simple_f(2, dimsf, dataspace1, error)

       CALL h5dcreate_f(file1_id, dsetname1,H5T_NATIVE_DOUBLE,dataspace1,dset1_id, error)
       CALL h5sclose_f(dataspace1, error)
     !
       !CALL h5dopen_f(file2_id, dsetname2, dset2_id, error)
     !
     ! Get  dataspace identifier.
     !
      count(1) = nrows_p
       count(2) = nc*nv
       dimsm(1)=nrows_p
       dimsm(2)=nc*nv
       dimsf(1) = s
       dimsf(2) = nc*nv
       !call h5screate_simple_f(2, dimsm, memspace, error)
       !allocate(data_out(count(1),count(2)))
       !data_out = real(H12)
       offset(1) = rank*nrows_p
       offset(2) = 0
       call h5screate_simple_f(2, count, memspace, error)
       CALL h5dget_space_f(dset1_id, dataspace1, error)
       CALL h5sselect_hyperslab_f (dataspace1, H5S_SELECT_SET_F, offset,&
               count, error)


       CALL h5pcreate_f(H5P_DATASET_XFER_F, fapl_id, error)
       CALL h5pset_dxpl_mpio_f(fapl_id, H5FD_MPIO_COLLECTIVE_F, error)
       call h5screate_simple_f(2, dimsm, memspace, error)
       CALL h5dwrite_f(dset1_id, H5T_NATIVE_DOUBLE, real(H13), dimsf , error, &
                    mem_space_id=memspace, file_space_id=dataspace1)
       call h5pclose_f(fapl_id, error)

       CALL h5sclose_f(memspace, error)
       CALL h5dclose_f(dset1_id, error)

       CALL h5fclose_f(file1_id, error)

       CALL h5close_f(error)
       call cpu_time(t3)

       print*,t2-t1,t3-t2


  end subroutine create_H13
    subroutine create_H23(nv,nc,mf_v,mf_c,wcoul0,s,vol,chunk)
      double precision :: fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16

      integer, intent(in) :: nv,nc,s,chunk
      integer :: nb,i,j,alpha,beta,gama,eta,m,n,imatrix,keyword,c1,c2,c,x,A
      double precision            :: element
      double precision ::wcoul0,vol,a1,b1
      complex(kind=8),dimension(:,:),allocatable    :: H23
      complex(kind=8), allocatable :: bse_mat_e(:,:,:,:)
      complex(kind=8), allocatable :: bse_mat_h(:,:,:,:)
      complex(kind=8), allocatable :: bse_mat(:,:,:,:)
      CHARACTER(LEN=11), PARAMETER :: filename1 = "singlet.h5" ! File name
      CHARACTER(LEN=5), PARAMETER :: dsetname1 = "CIS23"    ! Dataset name
      CHARACTER(LEN=11), PARAMETER :: filename2 = "singletX.h5" ! File name
      CHARACTER(LEN=5), PARAMETER :: dsetname2 = "HCISD"    ! Dataset name

      INTEGER, PARAMETER   :: DRANK = 3 ! Dataset rank 1 dimension for real or complex 2 is row-index 3 is col-index
      INTEGER(SIZE_T), PARAMETER :: NUMP = 2 ! Number of points selected
      complex(kind=8) ,DIMENSION(nv) ,intent(in) :: mf_v
      complex(kind=8),DiMENSION(nc),intent(in) :: mf_c
      INTEGER(HID_T) :: file1_id       ! File1 identifier
      INTEGER(HID_T) :: dset1_id       ! Dataset1 identifier
      INTEGER(HID_T) :: file2_id       ! File2 identifier
      INTEGER(HID_T) :: dset2_id       ! Dataset2 identifier
      INTEGER(HID_T) :: dataspace1     ! Dataspace identifier
          ! memspace identifier
      INTEGER(HID_T) :: dataspace2     ! Dataspace identifier
            ! memspace identifier
      INTEGER(HID_T) :: fapl_id

      INTEGER(HSIZE_T), DIMENSION(2) :: dimsm
                                                   ! Memory dataspace dimensions
      INTEGER(HSIZE_T), DIMENSION(2) :: dimsf
                                                   ! File dataspace dimensions
      

      INTEGER :: error,info  ! Error flag
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims
      integer(HSIZE_T) :: count(2), offset(2)

      integer(HID_T) :: filespace     ! Dataspace identifier in file
      integer(HID_T) :: memspace      ! Dataspace identifier in mem
      integer(HID_T) :: dset_id       ! Dataset identifier
      INTEGER(HID_T) :: plist_id      ! Property list identifier
      !!!!!!!!!!!mpivariables!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer :: nproc, chunk_size,itr_length
      integer , dimension(:,:) , allocatable :: ij_mat
      integer, dimension(:),allocatable :: row_ij
      integer :: comm, rank, numproc, ierror
      integer :: cij,cab,cijp,c1p
      !call MPI_INIT(ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)

     
       keyword = 6

       fac_hd=0
       fac_hx=0
       fac_ed=0
       fac_ex=0
      fac_1=0
      fac_2=0
      fac_3=0
      fac_4=0
      fac_5=0
      fac_6=0
      fac_7=0
      fac_8=0
      fac_9=0
      fac_10=0
      fac_11=0
      fac_12=0
      fac_13=0
      fac_14=0
      fac_15=0
      fac_16=0

       itr_length = nc*(nc-1)/2
       chunk_size = chunk
      print*,chunk_size,nproc

      allocate(bse_mat_h(nv,nv,nv,nv))
      allocate(bse_mat_e(nc,1,1,nc))
      allocate(bse_mat(nv,nv,2,nc))

      allocate(row_ij(chunk_size))
      A = (nv*(nv-1))/2
      allocate(ij_mat(2,chunk_size))
      c = rank+1
      nrows_p= chunk_size*A
      allocate(H23(nrows_p,nc*nv))
      H23 = cmplx(0.0,0.0)
      


      !print*,"proc_id+1",c
      call distribute_ij(rank,nc,A,chunk_size,ij_mat,row_ij)
      c1=0
      c2=0
      cij=0
      cab=0
      do imatrix = 3,4
       call fac_values(wcoul0,vol,keyword,imatrix,fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16)
         if (imatrix==3) then
            call load_Hhh(nv,bse_mat_h)
         end if
         do x = 1,chunk_size
            print*,x,ij_mat(1,x),ij_mat(2,x)

            i = ij_mat(1,x)
            j = ij_mat(2,x)
            cij = (x-1)*A

            cab =0

            if (imatrix .eq. 3) then
                call load_Hee(i,j,nc,nv,bse_mat_e)
                call load_Heh(i,j,nc,nv,imatrix,bse_mat)
            else if (imatrix .eq. 4) then
                call load_Heh(i,j,nc,nv,imatrix,bse_mat)
            end if
            do beta = 1,nv-1
               do alpha = beta+1,nv
                  cab = cab +1
                  c1 = cij+cab
                  c2 = 0
                  do m = 1,nc
                     do gama = 1,nv
                        c2 = c2 +1
                        if(j==m .and. beta==gama) then
                           if (imatrix .ne. 4) then
                               H23(c1,c2) = H23(c1,c2) + fac_1 * bse_mat(alpha,gama,1,m)
                               H23(c1,c2) = H23(c1,c2) + fac_2 * bse_mat(alpha,gama,1,m)
                               H23(c1,c2) = H23(c1,c2) + fac_9 * bse_mat(alpha,gama,1,m)
                               H23(c1,c2) = H23(c1,c2) + fac_10 * bse_mat(alpha,gama,1,m)
                               !print*,"3", i,j,alpha,beta,m,n,gama,c1,c2,H23(c1,c2)
                                    !print*, H(c1,c2),c1,c2
                           end if
                           if (imatrix == 4) then
                               H23(c1,c2) = H23(c1,c2) + fac_1 * bse_mat(alpha,gama,1,m)
                               H23(c1,c2) = H23(c1,c2) + fac_2 * bse_mat(alpha,gama,1,m)
                               H23(c1,c2) = H23(c1,c2) + fac_9 * bse_mat(alpha,gama,1,m)
                               H23(c1,c2) = H23(c1,c2) + fac_10 * bse_mat(alpha,gama,1,m)
                               !print*,"3", i,j,alpha,beta,m,n,gama,c1,c2,H23(c1,c2)
                           end if
                        end if

                        if(i==m .and. beta==gama) then
                           if (imatrix .ne. 4) then
                               H23(c1,c2)=H23(c1,c2) + fac_3 * bse_mat(alpha,gama,2,m)
                               H23(c1,c2)=H23(c1,c2) + fac_4 * bse_mat(alpha,gama,2,m)
                               H23(c1,c2)=H23(c1,c2) + fac_11 * bse_mat(alpha,gama,2,m)
                               H23(c1,c2)=H23(c1,c2) + fac_12 * bse_mat(alpha,gama,2,m)
                               !print*,"3", i,j,alpha,beta,m,n,gama,c1,c2,H23(c1,c2)
                                      !print*,"5", i,j,alpha,beta,m,n,gama,eta
                                      !print*, H(c1,c2),c1,c2
                           end if
                           if (imatrix == 4) then
                               H23(c1,c2)=H23(c1,c2) + fac_3 * bse_mat(alpha,gama,2,m)
                               H23(c1,c2)=H23(c1,c2) + fac_4 * bse_mat(alpha,gama,2,m)
                               H23(c1,c2)=H23(c1,c2) + fac_11 * bse_mat(alpha,gama,2,m)
                               H23(c1,c2)=H23(c1,c2) + fac_12 * bse_mat(alpha,gama,2,m)
                               !print*,"3","hey", i,j,alpha,beta,m,n,gama,c1,c2,H23(c1,c2)
                           end if
                        end if

                        if(j==m .and. alpha==gama) then
                           if (imatrix .ne. 4) then
                               !print*,"3", i,j,alpha,beta,m,n,gama,c1,c2,H23(c1,c2),bse_mat(beta,gama,1,m),&
                                                        !fac_5 * bse_mat(beta,gama,1,m),fac_6 * bse_mat(beta,gama,1,m)
                               H23(c1,c2)=H23(c1,c2) + fac_5 * bse_mat(beta,gama,1,m)
                               H23(c1,c2)=H23(c1,c2) + fac_6 * bse_mat(beta,gama,1,m)
                               H23(c1,c2)=H23(c1,c2) + fac_13 * bse_mat(beta,gama,1,m)
                               H23(c1,c2)=H23(c1,c2) + fac_14 * bse_mat(beta,gama,1,m)
                                !print*,"3", i,j,alpha,beta,m,n,gama,c1,c2,H23(c1,c2),bse_mat(beta,gama,1,m),&
                                                        !fac_5 * bse_mat(beta,gama,1,m),fac_6 * bse_mat(beta,gama,1,m)
                               !print*, fac_5,fac_6
                           end if
                           if (imatrix == 4) then
                               H23(c1,c2)=H23(c1,c2) + fac_5 * bse_mat(beta,gama,1,m)
                               H23(c1,c2)=H23(c1,c2) + fac_6 * bse_mat(beta,gama,1,m)
                               H23(c1,c2)=H23(c1,c2) + fac_13 * bse_mat(beta,gama,1,m)
                               H23(c1,c2)=H23(c1,c2) + fac_14 * bse_mat(beta,gama,1,m)
                               !print*,"3", i,j,alpha,beta,m,n,gama,c1,c2,H23(c1,c2)
                           end if
                        end if

                        if(i==m .and. alpha==gama) then
                           if (imatrix .ne. 4) then
                               H23(c1,c2)=H23(c1,c2) + fac_7 * bse_mat(beta,gama,2,m)
                               H23(c1,c2)=H23(c1,c2) + fac_8 * bse_mat(beta,gama,2,m)
                               H23(c1,c2)=H23(c1,c2) + fac_15 * bse_mat(beta,gama,2,m)
                               H23(c1,c2)=H23(c1,c2) + fac_16 * bse_mat(beta,gama,2,m)
                                    !print*,"9", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                               !print*,"3", i,j,alpha,beta,m,n,gama,c1,c2,H23(c1,c2),bse_mat(beta,gama,1,m)

                            end if
                            if (imatrix == 4) then
                                H23(c1,c2)=H23(c1,c2) + fac_7 * bse_mat(beta,gama,2,m)
                                H23(c1,c2)=H23(c1,c2) + fac_8 * bse_mat(beta,gama,2,m)
                                H23(c1,c2)=H23(c1,c2) + fac_15 * bse_mat(beta,gama,2,m)
                                H23(c1,c2)=H23(c1,c2) + fac_16 * bse_mat(beta,gama,2,m)
                                !print*,"3", i,j,alpha,beta,m,n,gama,c1,c2,H23(c1,c2),bse_mat(beta,gama,1,m)
                            end if
                        end if
                       
                     end do
                  end do
               end do
            end do
         end do
      end do
       call cpu_time(t2)
      CALL h5open_f(error)
      !CALL h5fcreate_f(filename1,  H5F_ACC_TRUNC_F, file1_id, error)


      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)

      !call h5fcreate_f(filename1, H5F_ACC_TRUNC_F, file1_id, error, &
       !           access_prp = plist_id)
      CALL h5fopen_f (filename1, H5F_ACC_RDWR_F, file1_id, error,access_prp = plist_id)


       !CALL h5fopen_f (filename2, H5F_ACC_RDWR_F, file2_id, error,access_prp = plist_id)
       CALL h5pclose_f(plist_id, error)
       dimsf(1) = s
       dimsf(2) = nc*nv
       CALL h5screate_simple_f(2, dimsf, dataspace1, error)

       CALL h5dcreate_f(file1_id, dsetname1,H5T_NATIVE_DOUBLE,dataspace1,dset1_id, error)
       CALL h5sclose_f(dataspace1, error)
     !
       !CALL h5dopen_f(file2_id, dsetname2, dset2_id, error)
     !
     ! Get  dataspace identifier.
     !
      count(1) = nrows_p
       count(2) = nc*nv
       dimsm(1)=nrows_p
       dimsm(2)=nc*nv
       dimsf(1) = s
       dimsf(2) = nc*nv
       !call h5screate_simple_f(2, dimsm, memspace, error)
       !allocate(data_out(count(1),count(2)))
       !data_out = real(H12)
       offset(1) = rank*nrows_p
       offset(2) = 0
       call h5screate_simple_f(2, count, memspace, error)
       CALL h5dget_space_f(dset1_id, dataspace1, error)
       CALL h5sselect_hyperslab_f (dataspace1, H5S_SELECT_SET_F, offset,&
               count, error)


       CALL h5pcreate_f(H5P_DATASET_XFER_F, fapl_id, error)
       CALL h5pset_dxpl_mpio_f(fapl_id, H5FD_MPIO_COLLECTIVE_F, error)
       call h5screate_simple_f(2, dimsm, memspace, error)
       CALL h5dwrite_f(dset1_id, H5T_NATIVE_DOUBLE, real(H23), dimsf , error, &
                    mem_space_id=memspace, file_space_id=dataspace1)
       call h5pclose_f(fapl_id, error)

       CALL h5sclose_f(memspace, error)
       CALL h5dclose_f(dset1_id, error)

       CALL h5fclose_f(file1_id, error)

       CALL h5close_f(error)
       call cpu_time(t3)

       print*,t2-t1,t3-t2


  end subroutine create_H23
  subroutine create_H4(nv,nc,mf_v,mf_c,wcoul0,s,vol,q,n4)
      double precision :: fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16

      integer, intent(in) :: nv,nc,s,q,n4
      integer :: nb,i,j,alpha,beta,gama,eta,m,n,imatrix,keyword,c1,c2,c,x,A
      double precision            :: element
      double precision ::wcoul0,vol,a1,b1
      complex(kind=8),dimension(:,:),allocatable     :: H4
      complex(kind=8), allocatable :: bse_mat_e(:,:,:,:)
      complex(kind=8), allocatable :: bse_mat_h(:,:,:,:)
      complex(kind=8), allocatable :: bse_mat(:,:,:,:)
      CHARACTER(LEN=11), PARAMETER :: filename1 = "singlet.h5" ! File name
      CHARACTER(LEN=5), PARAMETER :: dsetname1 = "CISD4"    ! Dataset name
      CHARACTER(LEN=11), PARAMETER :: filename2 = "singletX.h5" ! File name
      CHARACTER(LEN=5), PARAMETER :: dsetname2 = "HCISD"    ! Dataset name
      INTEGER, PARAMETER   :: DRANK = 3 ! Dataset rank 1 dimension for real or complex 2 is row-index 3 is col-index
      INTEGER(SIZE_T), PARAMETER :: NUMP = 2 ! Number of points selected
      complex(kind=8) ,DIMENSION(nv) ,intent(in) :: mf_v
      complex(kind=8),DiMENSION(nc),intent(in) :: mf_c
      INTEGER(HID_T) :: file1_id       ! File1 identifier
      INTEGER(HID_T) :: dset1_id       ! Dataset1 identifier
      INTEGER(HID_T) :: file2_id       ! File2 identifier
      INTEGER(HID_T) :: dset2_id       ! Dataset2 identifier
      INTEGER(HID_T) :: dataspace1     ! Dataspace identifier
          ! memspace identifier
      INTEGER(HID_T) :: dataspace2     ! Dataspace identifier
            ! memspace identifier
      INTEGER(HID_T) :: fapl_id

      INTEGER(HSIZE_T), DIMENSION(2) :: dimsm
                                                   ! Memory dataspace dimensions
      INTEGER(HSIZE_T), DIMENSION(2) :: dimsf
                                                   ! File dataspace dimensions
      
      INTEGER :: memrank = 1  ! Rank of the dataset in memory
      DOUBLE PRECISION,DIMENSION(:,:,:),allocatable :: bufnew

      INTEGER :: error,info  ! Error flag
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims
      integer(HSIZE_T) :: count(2), offset(2)

      integer(HID_T) :: filespace     ! Dataspace identifier in file
      integer(HID_T) :: memspace      ! Dataspace identifier in mem
      integer(HID_T) :: dset_id       ! Dataset identifier
      INTEGER(HID_T) :: plist_id      ! Property list identifier
      !!!!!!!!!!!mpivariables!!!!!!!!!!!!!!!!!!!!!!!!!!!
       integer :: nproc, chunk_size,itr_length
      integer , dimension(:,:) , allocatable :: ij_mat
      integer, dimension(:),allocatable :: row_ij
      integer :: comm, rank, numproc, ierror
      integer :: cij,cab,cijp,c1p
      !call MPI_INIT(ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)

       keyword = 7
         fac_hd=0
       fac_hx=0
       fac_ed=0
       fac_ex=0
      fac_1=0
      fac_2=0
      fac_3=0
      fac_4=0
      fac_5=0
      fac_6=0
      fac_7=0
      fac_8=0
      fac_9=0
      fac_10=0
      fac_11=0
      fac_12=0
      fac_13=0
      fac_14=0
      fac_15=0
      fac_16=0

       itr_length = nc*(nc-1)/2
      if (mod(nc,2)==0) then
         nproc = nc-1
         chunk_size = nc/2
      else
         nproc = nc
         chunk_size = (nc-1)/2
      end if
      print*,chunk_size,nproc

      allocate(bse_mat_h(nv,nv,nv,nv))
      allocate(bse_mat_e(nc,1,1,nc))
      allocate(bse_mat(nv,nv,2,nc))
      !allocate(row_ij(chunk_size))
      A = (nv*(nv-1))/2
      !allocate(ij_mat(2,chunk_size))
      c = rank+1
      nrows_p= A
      allocate(H4(nrows_p,q))
      H4 = cmplx(0.0,0.0)
     

      c1=0
      c2=0
      cij=0
      c1p =0
      cijp = 0
      cab=0
      do imatrix = 3,4
         call fac_values(wcoul0,vol,keyword,imatrix,fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16)

         i = c
         cij = 0
         cijp = (c-1)*A
         
         cab = 0
         if (imatrix==3) then
            call load_Hhh(nv,bse_mat_h)
         end if

         if (imatrix .eq. 3) then
                call load_Hee(i,i,nc,nv,bse_mat_e)
                call load_Heh(i,i,nc,nv,imatrix,bse_mat)
         else if (imatrix .eq. 4) then
                call load_Heh(i,i,nc,nv,imatrix,bse_mat)
         end if
         do beta = 1,nv-1
            do alpha = beta+1,nv
               cab = cab + 1
               c1 = cij + cab
               c1p = cijp + cab
               c2=0
               do m = 1,nc
                 do eta = 1,nv-1
                    do gama = eta+1,nv
                       c2 = c2 + 1
                       if (imatrix==3 .and. c1p==c2) then
                                   H4(c1,c2)=H4(c1,c2) + (2*mf_c(i)) -mf_v(beta) -mf_v(alpha)
                       end if
                       if(i==m) then
                           if (imatrix .ne. 4) then
                                     !print*,"bb",c1,c2, H(5,11),H(71,1),imatrix
                                      !print*,"1",H(c1,c2),bse_mat(eta,beta,alpha,gama), bse_mat(gama,beta,alpha,eta)
                               H4(c1,c2) = H4(c1,c2) + (fac_hd * bse_mat_h(eta,beta,alpha,gama)) &
                                        + (fac_hx * bse_mat_h(gama,beta,alpha,eta))   !<
                                     !print*,"1", i,j,alpha,beta,m,n,gama,eta
                                     !print*,c1,c2,H(c1,c2) 
                           end if
                        end if
                        if(alpha==gama .and. beta==eta) then
                           if (imatrix .ne. 4) then
                                     !print*,"2",H(c1,c2),bse_mat(n,j+nv,1,m)
                              H4(c1,c2)=H4(c1,c2) + (fac_ed *bse_mat_e(m,1,1,m)) &
                                       + (fac_ex * bse_mat_e(m,1,1,m))
                                    !print*,"2", i,j,alpha,beta,m,n,gama,eta
                                    !print*,c1,c2,H(c1,c2)
                           end if
                        end if
                        if(i==m .and. beta==gama) then
                           if (imatrix .ne. 4) then
                               H4(c1,c2)= H4(c1,c2) + fac_1 * bse_mat(alpha,eta,1,m)
                               H4(c1,c2)=H4(c1,c2) + fac_2 * bse_mat(alpha,eta,1,m)
                               H4(c1,c2)=H4(c1,c2) + fac_3 * bse_mat(alpha,eta,1,m)
                               H4(c1,c2) = H4(c1,c2) + fac_4 * bse_mat(alpha,eta,1,m)
                                    !print*,"3", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                           end if
                           if (imatrix == 4) then
                               H4(c1,c2)= H4(c1,c2) + fac_1 * bse_mat(alpha,eta,1,m)
                               H4(c1,c2)=H4(c1,c2) + fac_2 * bse_mat(alpha,eta,1,m)
                               H4(c1,c2)=H4(c1,c2) + fac_3 * bse_mat(alpha,eta,1,m)
                               H4(c1,c2)=H4(c1,c2) + fac_4 * bse_mat(alpha,eta,1,m)
                           end if
                        end if
                        if(i==m .and. alpha==eta) then
                          if (imatrix .ne. 4) then
                              H4(c1,c2)=H4(c1,c2) + fac_5 * bse_mat(beta,gama,1,m)
                              H4(c1,c2)=H4(c1,c2) + fac_6 * bse_mat(beta,gama,1,m)
                              H4(c1,c2)=H4(c1,c2) + fac_7 * bse_mat(beta,gama,1,m)
                              H4(c1,c2)=H4(c1,c2) + fac_8 * bse_mat(beta,gama,1,m)

                                    !print*, H(c1,c2),c1,c2
                          end if
                           if (imatrix == 4) then
                             H4(c1,c2)=H4(c1,c2) + fac_5 * bse_mat(beta,gama,1,m)
                             H4(c1,c2)=H4(c1,c2) + fac_6 * bse_mat(beta,gama,1,m)
                             H4(c1,c2)=H4(c1,c2) + fac_7 * bse_mat(beta,gama,1,m)
                             H4(c1,c2)=H4(c1,c2) + fac_8 * bse_mat(beta,gama,1,m)

                          end if
                        end if
                        if(i==m .and. beta==eta) then
                           if (imatrix .ne. 4) then

                              H4(c1,c2)=H4(c1,c2) + fac_9 * bse_mat(alpha,gama,1,m)
                              H4(c1,c2)=H4(c1,c2) + fac_10 * bse_mat(alpha,gama,1,m)
                              H4(c1,c2)=H4(c1,c2) + fac_11 * bse_mat(alpha,gama,1,m)
                              H4(c1,c2)=H4(c1,c2) + fac_12 * bse_mat(alpha,gama,1,m)
                                    !print*,"11", i,j,alpha,beta,m,n,gama,eta
                                    !print*,c1,c2,H(c1,c2)
                           end if
                           if (imatrix == 4) then
                              H4(c1,c2)=H4(c1,c2) + fac_9 * bse_mat(alpha,gama,1,m)
                              H4(c1,c2)=H4(c1,c2) + fac_10 * bse_mat(alpha,gama,1,m)
                              H4(c1,c2)=H4(c1,c2) + fac_11 * bse_mat(alpha,gama,1,m)
                              H4(c1,c2)=H4(c1,c2) + fac_12 * bse_mat(alpha,gama,1,m)

                           end if
                                  !print*,c1,c2,H(1,8),
                        end if
                        if(i==m .and. alpha==gama) then
                           if (imatrix .ne. 4) then
                               H4(c1,c2)=H4(c1,c2) + fac_13 * bse_mat(beta,eta,1,m)
                               H4(c1,c2)=H4(c1,c2) + fac_14 * bse_mat(beta,eta,1,m)
                               H4(c1,c2)=H4(c1,c2) + fac_15 * bse_mat(beta,eta,1,m)
                               H4(c1,c2)=H4(c1,c2) + fac_16 * bse_mat(beta,eta,1,m)
                                    !print*,"15", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                           end if
                           if (imatrix == 4) then
                              H4(c1,c2)=H4(c1,c2) + fac_13 * bse_mat(beta,eta,1,m)
                              H4(c1,c2)=H4(c1,c2) + fac_14 * bse_mat(beta,eta,1,m)
                              H4(c1,c2)=H4(c1,c2) + fac_15 * bse_mat(beta,eta,1,m)
                              H4(c1,c2)=H4(c1,c2) + fac_16 * bse_mat(beta,eta,1,m)
                           end if

                        end if
                     
                    end do 
                  end do
                end do
              end do
           end do
        end do
      call cpu_time(t2)
      CALL h5open_f(error)
      !CALL h5fcreate_f(filename1,  H5F_ACC_TRUNC_F, file1_id, error)


      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)

      !call h5fcreate_f(filename1, H5F_ACC_TRUNC_F, file1_id, error, &
       !           access_prp = plist_id)
      CALL h5fopen_f (filename1, H5F_ACC_RDWR_F, file1_id, error,access_prp = plist_id)


       !CALL h5fopen_f (filename2, H5F_ACC_RDWR_F, file2_id, error,access_prp = plist_id)
       CALL h5pclose_f(plist_id, error)
       dimsf(1) = q
       dimsf(2) = q
       CALL h5screate_simple_f(2, dimsf, dataspace1, error)

       CALL h5dcreate_f(file1_id, dsetname1,H5T_NATIVE_DOUBLE,dataspace1,dset1_id, error)
       CALL h5sclose_f(dataspace1, error)
     !
       !CALL h5dopen_f(file2_id, dsetname2, dset2_id, error)
     !
     ! Get  dataspace identifier.
     !
      count(1) = nrows_p
       count(2) = q
       dimsm(1)=nrows_p
       dimsm(2)=q
       dimsf(1) = q
       dimsf(2) = q
       !print*,"rank,count(1),count(2),offset(1),offset(2),nc*nv",rank,count(1),count(2),offset(1),offset(2),nc*nv

       !call h5screate_simple_f(2, dimsm, memspace, error)
       !allocate(data_out(count(1),count(2)))
       !data_out = real(H12)
       offset(1) = rank*nrows_p
       offset(2) = 0
       print*,"rank,count(1),count(2),offset(1),offset(2),nc*nv",rank,count(1),count(2),offset(1),offset(2),nc*nv

       call h5screate_simple_f(2, count, memspace, error)
       CALL h5dget_space_f(dset1_id, dataspace1, error)
       CALL h5sselect_hyperslab_f (dataspace1, H5S_SELECT_SET_F, offset,&
               count, error)


       CALL h5pcreate_f(H5P_DATASET_XFER_F, fapl_id, error)
       CALL h5pset_dxpl_mpio_f(fapl_id, H5FD_MPIO_COLLECTIVE_F, error)
       call h5screate_simple_f(2, dimsm, memspace, error)
       CALL h5dwrite_f(dset1_id, H5T_NATIVE_DOUBLE, real(H4), dimsf , error, &
                    mem_space_id=memspace, file_space_id=dataspace1)
       call h5pclose_f(fapl_id, error)
       
       CALL h5sclose_f(memspace, error)
       CALL h5dclose_f(dset1_id, error)

       CALL h5fclose_f(file1_id, error)

       CALL h5close_f(error)
       call cpu_time(t3)

       print*,t2-t1,t3-t2


  end subroutine create_H4
    subroutine create_H14(nv,nc,mf_v,mf_c,wcoul0,s,vol,q,n4,chunk)
      double precision :: fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16

      integer, intent(in) :: nv,nc,s,q,n4,chunk
      integer :: nb,i,j,alpha,beta,gama,eta,m,n,imatrix,keyword,c1,c2,c,x,A
      double precision            :: element
      double precision ::wcoul0,vol,a1,b1
      complex(kind=8) ,dimension(:,:),allocatable                          :: H14
      complex(kind=8), allocatable :: bse_mat_e(:,:,:,:)
      complex(kind=8), allocatable :: bse_mat_h(:,:,:,:)
      complex(kind=8), allocatable :: bse_mat(:,:,:,:)
      CHARACTER(LEN=11), PARAMETER :: filename1 = "singlet.h5" ! File name
      CHARACTER(LEN=5), PARAMETER :: dsetname1 = "CIS14"    ! Dataset name
      CHARACTER(LEN=11), PARAMETER :: filename2 = "singletX.h5" ! File name
      CHARACTER(LEN=5), PARAMETER :: dsetname2 = "HCISD"    ! Dataset name

      INTEGER, PARAMETER   :: DRANK = 3 ! Dataset rank 1 dimension for real or complex 2 is row-index 3 is col-index
      INTEGER(SIZE_T), PARAMETER :: NUMP = 2 ! Number of points selected
      complex(kind=8) ,DIMENSION(nv) ,intent(in) :: mf_v
      complex(kind=8),DiMENSION(nc),intent(in) :: mf_c
      INTEGER(HID_T) :: file1_id       ! File1 identifier
      INTEGER(HID_T) :: dset1_id       ! Dataset1 identifier
      INTEGER(HID_T) :: file2_id       ! File2 identifier
      INTEGER(HID_T) :: dset2_id       ! Dataset2 identifier
      INTEGER(HID_T) :: dataspace1     ! Dataspace identifier
          ! memspace identifier
      INTEGER(HID_T) :: dataspace2     ! Dataspace identifier
            ! memspace identifier
      INTEGER(HID_T) :: fapl_id

      INTEGER(HSIZE_T), DIMENSION(2) :: dimsm
                                                   ! Memory dataspace dimensions
      INTEGER(HSIZE_T), DIMENSION(2) :: dimsf
                                                   ! File dataspace dimensions
     
      INTEGER :: memrank = 1  ! Rank of the dataset in memory
      DOUBLE PRECISION,DIMENSION(:,:,:),allocatable :: bufnew
      integer(HSIZE_T) :: count(2), offset(2)


      INTEGER :: error,info  ! Error flag
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims

      integer(HID_T) :: filespace     ! Dataspace identifier in file
      integer(HID_T) :: memspace      ! Dataspace identifier in mem
      integer(HID_T) :: dset_id       ! Dataset identifier
      INTEGER(HID_T) :: plist_id      ! Property list identifier
       !!!!!!!!!!!mpivariables!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer :: nproc, chunk_size,itr_length
      integer , dimension(:,:) , allocatable :: ij_mat
      integer, dimension(:),allocatable :: row_ij
      integer :: comm, rank, numproc, ierror
      integer :: cij,cab
      !call MPI_INIT(ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)

       keyword = 8
       fac_hd=0
       fac_hx=0
       fac_ed=0
       fac_ex=0
      fac_1=0
      fac_2=0
      fac_3=0
      fac_4=0
      fac_5=0
      fac_6=0
      fac_7=0
      fac_8=0
      fac_9=0
      fac_10=0
      fac_11=0
      fac_12=0
      fac_13=0
      fac_14=0
      fac_15=0
      fac_16=0
       itr_length = nc*(nc-1)/2
       chunk_size = chunk
      print*,chunk_size,nproc

      allocate(bse_mat_h(nv,nv,nv,nv))
      allocate(bse_mat_e(nc,1,1,nc))
      allocate(bse_mat(nv,nv,2,nc))

      allocate(row_ij(chunk_size))
      A = (nv*(nv-1))/2
      allocate(ij_mat(2,chunk_size))
      c = rank+1
      nrows_p= chunk_size*A
      allocate(H14(nrows_p,q))

      H14 = cmplx(0.0,0.0)
      

      !print*,"proc_id+1",c
      call distribute_ij(rank,nc,A,chunk_size,ij_mat,row_ij)
      c1=0
      c2=0
      cij=0
      cab=0
      do imatrix = 3,4
       call fac_values(wcoul0,vol,keyword,imatrix,fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16)
         if (imatrix==3) then
            call load_Hhh(nv,bse_mat_h)
         end if
         do x = 1,chunk_size
            print*,x,ij_mat(1,x),ij_mat(2,x)

            i = ij_mat(1,x)
            j = ij_mat(2,x)
            cij = (x-1)*A
            cab =0

            if (imatrix .eq. 3) then
                call load_Hee(i,j,nc,nv,bse_mat_e)
                call load_Heh(i,j,nc,nv,imatrix,bse_mat)
            else if (imatrix .eq. 4) then
                call load_Heh(i,j,nc,nv,imatrix,bse_mat)
            end if
            do beta = 1,nv-1
               do alpha = beta+1,nv
                  cab = cab +1
                  c1 = cij+cab
                  c2 = 0
                   do m = 1,nc
                      do eta = 1,nv-1
                         do gama = eta+1,nv 
                            c2 = c2 + 1
                            if(alpha==gama .and. beta==eta) then
                               if (imatrix .ne. 4) then
                                     !print*,"2",H(c1,c2),bse_mat(n,2,1,m)
                                   H14(c1,c2)=H14(c1,c2) + (fac_ed *bse_mat_e(m,1,1,m)) &
                                       + (fac_ex * bse_mat_e(m,1,1,m))
                                    !print*,"2", i,j,alpha,beta,m,n,gama,eta
                                    !print*,c1,c2,H(c1,c2)
                               end if
                            end if 
                            if(j==m .and. beta==gama) then
                               if (imatrix .ne. 4) then
                                   H14(c1,c2)= H14(c1,c2) + fac_1 * bse_mat(alpha,eta,1,m)
                                   H14(c1,c2)=H14(c1,c2) + fac_2 * bse_mat(alpha,eta,1,m)
                                    !print*,"3", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                               end if
                               if (imatrix == 4) then 
                                   H14(c1,c2)= H14(c1,c2) + fac_1 * bse_mat(alpha,eta,1,m)
                                   H14(c1,c2)=H14(c1,c2) + fac_2 * bse_mat(alpha,eta,1,m)
                               end if
                            end if
                            if(i==m .and. beta==gama) then
                               if (imatrix .ne. 4) then
                                   H14(c1,c2)=H14(c1,c2) + fac_3 * bse_mat(alpha,eta,2,m)
                                   H14(c1,c2)=H14(c1,c2) + fac_4 * bse_mat(alpha,eta,2,m)
                                      !print*,"5", i,j,alpha,beta,m,n,gama,eta
                                      !print*, H(c1,c2),c1,c2
                               end if 
                               if (imatrix == 4) then 
                                   H14(c1,c2)=H14(c1,c2) + fac_3 * bse_mat(alpha,eta,2,m)
                                   H14(c1,c2)=H14(c1,c2) + fac_4 * bse_mat(alpha,eta,2,m)                               
                               end if
                            end if 
                            if(j==m .and. alpha==eta) then
                               if (imatrix .ne. 4) then
                                   H14(c1,c2)=H14(c1,c2) + fac_5 * bse_mat(beta,gama,1,m)
                                   H14(c1,c2)=H14(c1,c2) + fac_6 * bse_mat(beta,gama,1,m)
                                    !print*, H(c1,c2),c1,c2
                               end if
                               if (imatrix == 4) then 
                                   H14(c1,c2)=H14(c1,c2) + fac_5 * bse_mat(beta,gama,1,m)
                                   H14(c1,c2)=H14(c1,c2) + fac_6 * bse_mat(beta,gama,1,m)                                
                               end if 
                            end if 
                            if(i==m .and. alpha==eta) then
                               if (imatrix .ne. 4) then
                                   H14(c1,c2)=H14(c1,c2) + fac_7 * bse_mat(beta,gama,2,m)
                                   H14(c1,c2)=H14(c1,c2) + fac_8 * bse_mat(beta,gama,2,m)
                                    !print*,"9", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                               end if
                               if (imatrix == 4) then 
                                  H14(c1,c2)=H14(c1,c2) + fac_7 * bse_mat(beta,gama,2,m)
                                  H14(c1,c2)=H14(c1,c2) + fac_8 * bse_mat(beta,gama,2,m)                                 
                               end if  
                            end if 
                            if(j==m .and. beta==eta) then
                               if (imatrix .ne. 4) then
                                    !print*,"11",H(c1,c2),bse_mat(alpha,gama,1,m)
                                  H14(c1,c2)=H14(c1,c2) +fac_9 * bse_mat(alpha,gama,1,m)
                                  H14(c1,c2)=H14(c1,c2) + fac_10 * bse_mat(alpha,gama,1,m)
                                    !print*,"11", i,j,alpha,beta,m,n,gama,eta
                                    !print*,c1,c2,H(c1,c2)
                               end if 
                               if (imatrix == 4) then 
                                  H14(c1,c2)=H14(c1,c2) +fac_9 * bse_mat(alpha,gama,1,m)
                                  H14(c1,c2)=H14(c1,c2) + fac_10 * bse_mat(alpha,gama,1,m)                               
                               end if  
                                  !print*,c1,c2,H(1,8),
                            end if 
                            if(i==m .and. beta==eta) then
                               if (imatrix .ne. 4) then
                                  H14(c1,c2)=H14(c1,c2) + fac_11 * bse_mat(alpha,gama,2,m)
                                  H14(c1,c2)=H14(c1,c2) + fac_12 * bse_mat(alpha,gama,2,m)
                                    !print*,"13", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                               end if 
                               if (imatrix == 4) then 
                                   H14(c1,c2)=H14(c1,c2) + fac_11 * bse_mat(alpha,gama,2,m)
                                   H14(c1,c2)=H14(c1,c2) + fac_12 * bse_mat(alpha,gama,2,m)                               
                               end if 
                            end if 
                            if(j==m .and. alpha==gama) then
                               if (imatrix .ne. 4) then
                                   H14(c1,c2)=H14(c1,c2) + fac_13 * bse_mat(beta,eta,1,m)
                                   H14(c1,c2)=H14(c1,c2) + fac_14 * bse_mat(beta,eta,1,m)
                                    !print*,"15", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                               end if 
                               if (imatrix == 4) then 
                                   H14(c1,c2)=H14(c1,c2) + fac_13 * bse_mat(beta,eta,1,m)
                                   H14(c1,c2)=H14(c1,c2) + fac_14 * bse_mat(beta,eta,1,m)                                  
                               end if 
                             end if 
                            if(i==m.and. alpha==gama) then
                               if (imatrix .ne. 4) then
                                   H14(c1,c2)=H14(c1,c2) + fac_15 * bse_mat(beta,eta,2,m)
                                   H14(c1,c2)=H14(c1,c2) + fac_16 * bse_mat(beta,eta,2,m)
                                    !print*,"17", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                               end if
                               if (imatrix == 4) then 
                                   H14(c1,c2)=H14(c1,c2) + fac_15 * bse_mat(beta,eta,2,m)
                                   H14(c1,c2)=H14(c1,c2) + fac_16 * bse_mat(beta,eta,2,m)
                                                                  
                               end if  
                            end if 
                            
                         end do
                      end do
                   end do
               end do
            end do
         end do
      end do
      CALL h5open_f(error)
      !CALL h5fcreate_f(filename1,  H5F_ACC_TRUNC_F, file1_id, error)


      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)

      !call h5fcreate_f(filename1, H5F_ACC_TRUNC_F, file1_id, error, &
       !           access_prp = plist_id)
      CALL h5fopen_f (filename1, H5F_ACC_RDWR_F, file1_id, error,access_prp = plist_id)


       !CALL h5fopen_f (filename2, H5F_ACC_RDWR_F, file2_id, error,access_prp = plist_id)
       CALL h5pclose_f(plist_id, error)
       dimsf(1) = s
       dimsf(2) = q
       CALL h5screate_simple_f(2, dimsf, dataspace1, error)

       CALL h5dcreate_f(file1_id, dsetname1,H5T_NATIVE_DOUBLE,dataspace1,dset1_id, error)
       CALL h5sclose_f(dataspace1, error)
     !
       !CALL h5dopen_f(file2_id, dsetname2, dset2_id, error)
     !
     ! Get  dataspace identifier.
     !
      count(1) = nrows_p
       count(2) = q
       dimsm(1)=nrows_p
       dimsm(2)=q
       dimsf(1) = s
       dimsf(2) = q
       !call h5screate_simple_f(2, dimsm, memspace, error)
       !allocate(data_out(count(1),count(2)))
       !data_out = real(H12)
       offset(1) = rank*nrows_p
       offset(2) = 0
       call h5screate_simple_f(2, count, memspace, error)
       CALL h5dget_space_f(dset1_id, dataspace1, error)
       CALL h5sselect_hyperslab_f (dataspace1, H5S_SELECT_SET_F, offset,&
               count, error)


       CALL h5pcreate_f(H5P_DATASET_XFER_F, fapl_id, error)
       CALL h5pset_dxpl_mpio_f(fapl_id, H5FD_MPIO_COLLECTIVE_F, error)
       call h5screate_simple_f(2, dimsm, memspace, error)
       CALL h5dwrite_f(dset1_id, H5T_NATIVE_DOUBLE, real(H14), dimsf , error, &
                    mem_space_id=memspace, file_space_id=dataspace1)
       call h5pclose_f(fapl_id, error)

       CALL h5sclose_f(memspace, error)
       !call h5pclose_f(fapl_id, error)

       CALL h5dclose_f(dset1_id, error)
       CALL h5fclose_f(file1_id, error)

       CALL h5close_f(error)
       call cpu_time(t3)

       print*,t2-t1,t3-t2



 
 end subroutine create_H14
  subroutine create_H24(nv,nc,mf_v,mf_c,wcoul0,s,vol,q,n4,chunk)
      double precision :: fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16

      integer, intent(in) :: nv,nc,s,q,n4,chunk
      integer :: nb,i,j,alpha,beta,gama,eta,m,n,imatrix,keyword,c1,c2,c,x,A
      double precision            :: element
      double precision ::wcoul0,vol,a1,b1
      complex(kind=8) ,dimension(:,:),allocatable                          :: H24
      complex(kind=8), allocatable :: bse_mat_e(:,:,:,:)
      complex(kind=8), allocatable :: bse_mat_h(:,:,:,:)
      complex(kind=8), allocatable :: bse_mat(:,:,:,:)
      CHARACTER(LEN=11), PARAMETER :: filename1 = "singlet.h5" ! File name
      CHARACTER(LEN=5), PARAMETER :: dsetname1 = "CIS24"    ! Dataset name
      CHARACTER(LEN=11), PARAMETER :: filename2 = "singletX.h5" ! File name
      CHARACTER(LEN=5), PARAMETER :: dsetname2 = "HCISD"    ! Dataset name

      INTEGER, PARAMETER   :: DRANK = 3 ! Dataset rank 1 dimension for real or complex 2 is row-index 3 is col-index
      INTEGER(SIZE_T), PARAMETER :: NUMP = 2 ! Number of points selected
      complex(kind=8) ,DIMENSION(nv) ,intent(in) :: mf_v
      complex(kind=8),DiMENSION(nc),intent(in) :: mf_c
      INTEGER(HID_T) :: file1_id       ! File1 identifier
      INTEGER(HID_T) :: dset1_id       ! Dataset1 identifier
      INTEGER(HID_T) :: file2_id       ! File2 identifier
      INTEGER(HID_T) :: dset2_id       ! Dataset2 identifier
      INTEGER(HID_T) :: dataspace1     ! Dataspace identifier
          ! memspace identifier
      INTEGER(HID_T) :: dataspace2     ! Dataspace identifier
            ! memspace identifier
      INTEGER(HID_T) :: fapl_id

      INTEGER(HSIZE_T), DIMENSION(2) :: dimsm
                                                   ! Memory dataspace dimensions
      INTEGER(HSIZE_T), DIMENSION(2) :: dimsf
                                                   ! File dataspace dimensions
     
      INTEGER :: memrank = 1  ! Rank of the dataset in memory
      DOUBLE PRECISION,DIMENSION(:,:,:),allocatable :: bufnew
      integer(HSIZE_T) :: count(2), offset(2)


      INTEGER :: error,info  ! Error flag
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims

      integer(HID_T) :: filespace     ! Dataspace identifier in file
      integer(HID_T) :: memspace      ! Dataspace identifier in mem
      integer(HID_T) :: dset_id       ! Dataset identifier
      INTEGER(HID_T) :: plist_id      ! Property list identifier
       !!!!!!!!!!!mpivariables!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer :: nproc, chunk_size,itr_length
      integer , dimension(:,:) , allocatable :: ij_mat
      integer, dimension(:),allocatable :: row_ij
      integer :: comm, rank, numproc, ierror
      integer :: cij,cab
      !call MPI_INIT(ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)

       keyword = 9
       fac_hd=0
       fac_hx=0
       fac_ed=0
       fac_ex=0
      fac_1=0
      fac_2=0
      fac_3=0
      fac_4=0
      fac_5=0
      fac_6=0
      fac_7=0
      fac_8=0
      fac_9=0
      fac_10=0
      fac_11=0
      fac_12=0
      fac_13=0
      fac_14=0
      fac_15=0
      fac_16=0
       itr_length = nc*(nc-1)/2
      chunk_size = chunk
      print*,chunk_size,nproc

      allocate(bse_mat_h(nv,nv,nv,nv))
      allocate(bse_mat_e(nc,1,1,nc))
      allocate(bse_mat(nv,nv,2,nc))

      allocate(row_ij(chunk_size))
      A = (nv*(nv-1))/2
      allocate(ij_mat(2,chunk_size))
      c = rank+1
      nrows_p= chunk_size*A
      allocate(H24(nrows_p,q))

      H24 = cmplx(0.0,0.0)
      

      !print*,"proc_id+1",c
      call distribute_ij(rank,nc,A,chunk_size,ij_mat,row_ij)
      c1=0
      c2=0
      cij=0
      cab=0
      do imatrix = 3,4
       call fac_values(wcoul0,vol,keyword,imatrix,fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16)
         if (imatrix==3) then
            call load_Hhh(nv,bse_mat_h)
         end if
         do x = 1,chunk_size
            print*,x,ij_mat(1,x),ij_mat(2,x)

            i = ij_mat(1,x)
            j = ij_mat(2,x)
            cij = (x-1)*A
            cab =0

            if (imatrix .eq. 3) then
                call load_Hee(i,j,nc,nv,bse_mat_e)
                call load_Heh(i,j,nc,nv,imatrix,bse_mat)
            else if (imatrix .eq. 4) then
                call load_Heh(i,j,nc,nv,imatrix,bse_mat)
            end if
            do beta = 1,nv-1
               do alpha = beta+1,nv
                  cab = cab +1
                  c1 = cij+cab
                  c2 = 0
                   do m = 1,nc
                      do eta = 1,nv-1
                         do gama = eta+1,nv 
                            c2 = c2 + 1
                            if(alpha==gama .and. beta==eta) then
                               if (imatrix .ne. 4) then
                                     !print*,"2",H(c1,c2),bse_mat(n,2,1,m)
                                   H24(c1,c2)=H24(c1,c2) + (fac_ed *bse_mat_e(m,1,1,m)) &
                                       + (fac_ex * bse_mat_e(m,1,1,m))
                                    !print*,"2", i,j,alpha,beta,m,n,gama,eta
                                    !print*,c1,c2,H(c1,c2)
                               end if
                            end if 
                            if(j==m .and. beta==gama) then
                               if (imatrix .ne. 4) then
                                   H24(c1,c2)= H24(c1,c2) + fac_1 * bse_mat(alpha,eta,1,m)
                                   H24(c1,c2)=H24(c1,c2) + fac_2 * bse_mat(alpha,eta,1,m)
                                    !print*,"3", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                               end if
                               if (imatrix == 4) then 
                                   H24(c1,c2)= H24(c1,c2) + fac_1 * bse_mat(alpha,eta,1,m)
                                   H24(c1,c2)=H24(c1,c2) + fac_2 * bse_mat(alpha,eta,1,m)
                               end if
                            end if
                            if(i==m .and. beta==gama) then
                               if (imatrix .ne. 4) then
                                   H24(c1,c2)=H24(c1,c2) + fac_3 * bse_mat(alpha,eta,2,m)
                                   H24(c1,c2)=H24(c1,c2) + fac_4 * bse_mat(alpha,eta,2,m)
                                      !print*,"5", i,j,alpha,beta,m,n,gama,eta
                                      !print*, H(c1,c2),c1,c2
                               end if 
                               if (imatrix == 4) then 
                                   H24(c1,c2)=H24(c1,c2) + fac_3 * bse_mat(alpha,eta,2,m)
                                   H24(c1,c2)=H24(c1,c2) + fac_4 * bse_mat(alpha,eta,2,m)                               
                               end if
                            end if 
                            if(j==m .and. alpha==eta) then
                               if (imatrix .ne. 4) then
                                   H24(c1,c2)=H24(c1,c2) + fac_5 * bse_mat(beta,gama,1,m)
                                   H24(c1,c2)=H24(c1,c2) + fac_6 * bse_mat(beta,gama,1,m)
                                    !print*, H(c1,c2),c1,c2
                               end if
                               if (imatrix == 4) then 
                                   H24(c1,c2)=H24(c1,c2) + fac_5 * bse_mat(beta,gama,1,m)
                                   H24(c1,c2)=H24(c1,c2) + fac_6 * bse_mat(beta,gama,1,m)                                
                               end if 
                            end if 
                            if(i==m .and. alpha==eta) then
                               if (imatrix .ne. 4) then
                                   H24(c1,c2)=H24(c1,c2) + fac_7 * bse_mat(beta,gama,2,m)
                                   H24(c1,c2)=H24(c1,c2) + fac_8 * bse_mat(beta,gama,2,m)
                                    !print*,"9", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                               end if
                               if (imatrix == 4) then 
                                  H24(c1,c2)=H24(c1,c2) + fac_7 * bse_mat(beta,gama,2,m)
                                  H24(c1,c2)=H24(c1,c2) + fac_8 * bse_mat(beta,gama,2,m)                                 
                               end if  
                            end if 
                            if(j==m .and. beta==eta) then
                               if (imatrix .ne. 4) then
                                    !print*,"11",H(c1,c2),bse_mat(alpha,gama,1,m)
                                  H24(c1,c2)=H24(c1,c2) +fac_9 * bse_mat(alpha,gama,1,m)
                                  H24(c1,c2)=H24(c1,c2) + fac_10 * bse_mat(alpha,gama,1,m)
                                    !print*,"11", i,j,alpha,beta,m,n,gama,eta
                                    !print*,c1,c2,H(c1,c2)
                               end if 
                               if (imatrix == 4) then 
                                  H24(c1,c2)=H24(c1,c2) +fac_9 * bse_mat(alpha,gama,1,m)
                                  H24(c1,c2)=H24(c1,c2) + fac_10 * bse_mat(alpha,gama,1,m)                               
                               end if  
                                  !print*,c1,c2,H(1,8),
                            end if 
                            if(i==m .and. beta==eta) then
                               if (imatrix .ne. 4) then
                                  H24(c1,c2)=H24(c1,c2) + fac_11 * bse_mat(alpha,gama,2,m)
                                  H24(c1,c2)=H24(c1,c2) + fac_12 * bse_mat(alpha,gama,2,m)
                                    !print*,"13", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                               end if 
                               if (imatrix == 4) then 
                                   H24(c1,c2)=H24(c1,c2) + fac_11 * bse_mat(alpha,gama,2,m)
                                   H24(c1,c2)=H24(c1,c2) + fac_12 * bse_mat(alpha,gama,2,m)                               
                               end if 
                            end if 
                            if(j==m .and. alpha==gama) then
                               if (imatrix .ne. 4) then
                                   H24(c1,c2)=H24(c1,c2) + fac_13 * bse_mat(beta,eta,1,m)
                                   H24(c1,c2)=H24(c1,c2) + fac_14 * bse_mat(beta,eta,1,m)
                                    !print*,"15", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                               end if 
                               if (imatrix == 4) then 
                                   H24(c1,c2)=H24(c1,c2) + fac_13 * bse_mat(beta,eta,1,m)
                                   H24(c1,c2)=H24(c1,c2) + fac_14 * bse_mat(beta,eta,1,m)                                  
                               end if 
                             end if 
                            if(i==m.and. alpha==gama) then
                               if (imatrix .ne. 4) then
                                   H24(c1,c2)=H24(c1,c2) + fac_15 * bse_mat(beta,eta,2,m)
                                   H24(c1,c2)=H24(c1,c2) + fac_16 * bse_mat(beta,eta,2,m)
                                    !print*,"17", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                               end if
                               if (imatrix == 4) then 
                                   H24(c1,c2)=H24(c1,c2) + fac_15 * bse_mat(beta,eta,2,m)
                                   H24(c1,c2)=H24(c1,c2) + fac_16 * bse_mat(beta,eta,2,m)
                                                                  
                               end if  
                            end if 
                            
                         end do
                      end do
                   end do
               end do
            end do
         end do
      end do
      call cpu_time(t2)
      CALL h5open_f(error)
      !CALL h5fcreate_f(filename1,  H5F_ACC_TRUNC_F, file1_id, error)


      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)

      !call h5fcreate_f(filename1, H5F_ACC_TRUNC_F, file1_id, error, &
       !           access_prp = plist_id)
      CALL h5fopen_f (filename1, H5F_ACC_RDWR_F, file1_id, error,access_prp = plist_id)


       !CALL h5fopen_f (filename2, H5F_ACC_RDWR_F, file2_id, error,access_prp = plist_id)
       CALL h5pclose_f(plist_id, error)
       dimsf(1) = s
       dimsf(2) = q
       CALL h5screate_simple_f(2, dimsf, dataspace1, error)

       CALL h5dcreate_f(file1_id, dsetname1,H5T_NATIVE_DOUBLE,dataspace1,dset1_id, error)
       CALL h5sclose_f(dataspace1, error)
     !
       !CALL h5dopen_f(file2_id, dsetname2, dset2_id, error)
     !
     ! Get  dataspace identifier.
     !
      count(1) = nrows_p
       count(2) = q
       dimsm(1)=nrows_p
       dimsm(2)=q
       dimsf(1) = s
       dimsf(2) = q
       !call h5screate_simple_f(2, dimsm, memspace, error)
       !allocate(data_out(count(1),count(2)))
       !data_out = real(H12)
       offset(1) = rank*nrows_p
       offset(2) = 0
       call h5screate_simple_f(2, count, memspace, error)
       CALL h5dget_space_f(dset1_id, dataspace1, error)
       CALL h5sselect_hyperslab_f (dataspace1, H5S_SELECT_SET_F, offset,&
               count, error)


       CALL h5pcreate_f(H5P_DATASET_XFER_F, fapl_id, error)
       CALL h5pset_dxpl_mpio_f(fapl_id, H5FD_MPIO_COLLECTIVE_F, error)
       call h5screate_simple_f(2, dimsm, memspace, error)
       CALL h5dwrite_f(dset1_id, H5T_NATIVE_DOUBLE, real(H24), dimsf , error, &
                    mem_space_id=memspace, file_space_id=dataspace1)
       call h5pclose_f(fapl_id, error)

       CALL h5sclose_f(memspace, error)
       !call h5pclose_f(fapl_id, error)

       CALL h5dclose_f(dset1_id, error)
       CALL h5fclose_f(file1_id, error)

       CALL h5close_f(error)
       call cpu_time(t3)

       print*,t2-t1,t3-t2



 
 end subroutine create_H24

    subroutine create_H34(nv,nc,mf_v,mf_c,wcoul0,s,vol,q,n4)
      double precision :: fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16

      integer, intent(in) :: nv,nc,s,q,n4
      integer :: nb,i,j,alpha,beta,gama,eta,m,n,imatrix,keyword,c1,c2,c,x,A
      double precision            :: element
      double precision ::wcoul0,vol,a1,b1
      complex(kind=8),dimension(:,:),allocatable                          :: H34
      complex(kind=8), allocatable :: bse_mat_e(:,:,:,:)
      complex(kind=8), allocatable :: bse_mat_h(:,:,:,:)
      complex(kind=8), allocatable :: bse_mat(:,:,:,:)
      CHARACTER(LEN=11), PARAMETER :: filename1 = "singlet.h5" ! File name
      CHARACTER(LEN=5), PARAMETER :: dsetname1 = "CIS34"    ! Dataset name
      CHARACTER(LEN=11), PARAMETER :: filename2 = "singletX.h5" ! File name
      CHARACTER(LEN=5), PARAMETER :: dsetname2 = "HCISD"    ! Dataset name
      INTEGER, PARAMETER   :: DRANK = 3 ! Dataset rank 1 dimension for real or complex 2 is row-index 3 is col-index
      INTEGER(SIZE_T), PARAMETER :: NUMP = 2 ! Number of points selected
      complex(kind=8) ,DIMENSION(nv) ,intent(in) :: mf_v
      complex(kind=8),DiMENSION(nc),intent(in) :: mf_c
      INTEGER(HID_T) :: file1_id       ! File1 identifier
      INTEGER(HID_T) :: dset1_id       ! Dataset1 identifier
      INTEGER(HID_T) :: file2_id       ! File2 identifier
      INTEGER(HID_T) :: dset2_id       ! Dataset2 identifier
      INTEGER(HID_T) :: dataspace1     ! Dataspace identifier
          ! memspace identifier
      INTEGER(HID_T) :: dataspace2     ! Dataspace identifier
            ! memspace identifier
      INTEGER(HID_T) :: fapl_id

      INTEGER(HSIZE_T), DIMENSION(2) :: dimsm
                                                   ! Memory dataspace dimensions
      INTEGER(HSIZE_T), DIMENSION(2) :: dimsf
                                                   ! File dataspace dimensions
      INTEGER(HSIZE_T), DIMENSION(DRANK,NUMP) :: coord ! Elements coordinates
                                                      ! in the file
      DOUBLE PRECISION, DIMENSION(2) :: val =(/3.14,1.414/)  ! Values to write
      INTEGER :: memrank = 1  ! Rank of the dataset in memory
      DOUBLE PRECISION,DIMENSION(:,:,:),allocatable :: bufnew

      INTEGER :: error,info  ! Error flag
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims
      integer(hsize_t) :: count(2),offset(2)

      integer(HID_T) :: filespace     ! Dataspace identifier in file
      integer(HID_T) :: memspace      ! Dataspace identifier in mem
      integer(HID_T) :: dset_id       ! Dataset identifier
      INTEGER(HID_T) :: plist_id      ! Property list identifier
       !!!!!!!!!!!mpivariables!!!!!!!!!!!!!!!!!!!!!!!!!!!
       integer :: nproc, chunk_size,itr_length
      integer , dimension(:,:) , allocatable :: ij_mat
      integer, dimension(:),allocatable :: row_ij
      integer :: comm, rank, numproc, ierror
      integer :: cij,cab
      !call MPI_INIT(ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)

     
       keyword = 10
         fac_hd=0
       fac_hx=0
       fac_ed=0
       fac_ex=0
      fac_1=0
      fac_2=0
      fac_3=0
      fac_4=0
      fac_5=0
      fac_6=0
      fac_7=0
      fac_8=0
      fac_9=0
      fac_10=0
      fac_11=0
      fac_12=0
      fac_13=0
      fac_14=0
      fac_15=0
      fac_16=0

         itr_length = nc*(nc-1)/2
      if (mod(nc,2)==0) then
         nproc = nc-1
         chunk_size = nc/2
      else
         nproc = nc
         chunk_size = (nc-1)/2
      end if
      print*,chunk_size,nproc

      allocate(bse_mat_h(nv,nv,nv,nv))
      allocate(bse_mat_e(nc,1,1,nc))
      allocate(bse_mat(nv,nv,2,nc))
      !allocate(row_ij(chunk_size))
      A = nv
      !allocate(ij_mat(2,chunk_size))
      c = rank+1
      nrows_p= A
      allocate(H34(nrows_p,q))
      H34 = cmplx(0.0,0.0)
      
      c1=0
      c2=0
      cij=0
      cab=0
      do imatrix = 3,4
         call fac_values(wcoul0,vol,keyword,imatrix,fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16)

         i = c
         cij = 0
         cab = 0
         if (imatrix==3) then
            call load_Hhh(nv,bse_mat_h)
         end if

         if (imatrix .eq. 3) then
                call load_Hee(i,i,nc,nv,bse_mat_e)
                call load_Heh(i,i,nc,nv,imatrix,bse_mat)
         else if (imatrix .eq. 4) then
                call load_Heh(i,i,nc,nv,imatrix,bse_mat)
         end if

         do alpha = 1,nv
            cab = cab + 1
            c1 = cij + cab
            c2 = 0
            do m = 1,nc
               do eta = 1,nv-1
                   do gama = eta+1,nv
                      c2 = c2 + 1
                      if(i==m) then
                           if (imatrix .ne. 4) then
                                     !print*,"bb",c1,c2, H(5,11),H(71,1),imatrix
                                      !print*,"1",H(c1,c2),bse_mat(eta,beta,alpha,gama), bse_mat(gama,beta,alpha,eta)
                               H34(c1,c2) = H34(c1,c2) + (fac_hd * bse_mat_h(eta,alpha,alpha,gama)) &
                                        + (fac_hx * bse_mat_h(gama,alpha,alpha,eta))   !<
                                     !print*,"1", i,j,alpha,beta,m,n,gama,eta
                                     !print*,c1,c2,H(c1,c2) 
                           end if
                      end if
                      if(i==m .and. alpha==gama) then
                           if (imatrix .ne. 4) then
                               H34(c1,c2)= H34(c1,c2) + fac_1 * bse_mat(alpha,eta,1,m)
                               H34(c1,c2)= H34(c1,c2) + fac_2 * bse_mat(alpha,eta,1,m)
                               H34(c1,c2)= H34(c1,c2) + fac_3 * bse_mat(alpha,eta,1,m)
                               H34(c1,c2)= H34(c1,c2) + fac_4 * bse_mat(alpha,eta,1,m)
                               H34(c1,c2)=H34(c1,c2) + fac_13 * bse_mat(alpha,eta,1,m)
                               H34(c1,c2)=H34(c1,c2) + fac_14 * bse_mat(alpha,eta,1,m)
                               H34(c1,c2)=H34(c1,c2) + fac_15 * bse_mat(alpha,eta,1,m)
                               H34(c1,c2)=H34(c1,c2) + fac_16 * bse_mat(alpha,eta,1,m)
                                    !print*,"3", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                           end if
                           if (imatrix == 4) then
                              H34(c1,c2)= H34(c1,c2) + fac_1 * bse_mat(alpha,eta,1,m)
                              H34(c1,c2)= H34(c1,c2) + fac_2 * bse_mat(alpha,eta,1,m)
                              H34(c1,c2)= H34(c1,c2) + fac_3 * bse_mat(alpha,eta,1,m)
                              H34(c1,c2)= H34(c1,c2) + fac_4 * bse_mat(alpha,eta,1,m)
                              H34(c1,c2)=H34(c1,c2) + fac_13 * bse_mat(alpha,eta,1,m)
                              H34(c1,c2)=H34(c1,c2) + fac_14 * bse_mat(alpha,eta,1,m)
                              H34(c1,c2)=H34(c1,c2) + fac_15 * bse_mat(alpha,eta,1,m)
                              H34(c1,c2)=H34(c1,c2) + fac_16 * bse_mat(alpha,eta,1,m)

                           end if
                      end if
                       if(i==m .and. alpha==eta) then
                           if (imatrix .ne. 4) then
                               H34(c1,c2)=H34(c1,c2) + fac_5 * bse_mat(alpha,gama,1,m)
                               H34(c1,c2)=H34(c1,c2) + fac_6 * bse_mat(alpha,gama,1,m)
                               H34(c1,c2)=H34(c1,c2) + fac_7 * bse_mat(alpha,gama,1,m)
                               H34(c1,c2)=H34(c1,c2) + fac_8 * bse_mat(alpha,gama,1,m)
                               H34(c1,c2)=H34(c1,c2) + fac_9 * bse_mat(alpha,gama,1,m)
                               H34(c1,c2)=H34(c1,c2) + fac_10 * bse_mat(alpha,gama,1,m)
                               H34(c1,c2)=H34(c1,c2) + fac_11 * bse_mat(alpha,gama,1,m)
                               H34(c1,c2)=H34(c1,c2) + fac_12 * bse_mat(alpha,gama,1,m)
                                    !print*, H(c1,c2),c1,c2
                           end if
                           if (imatrix == 4) then
                               H34(c1,c2)=H34(c1,c2) + fac_5 * bse_mat(alpha,gama,1,m)
                               H34(c1,c2)=H34(c1,c2) + fac_6 * bse_mat(alpha,gama,1,m)
                               H34(c1,c2)=H34(c1,c2) + fac_7 * bse_mat(alpha,gama,1,m)
                               H34(c1,c2)=H34(c1,c2) + fac_8 * bse_mat(alpha,gama,1,m)
                               H34(c1,c2)=H34(c1,c2) + fac_9 * bse_mat(alpha,gama,1,m)
                               H34(c1,c2)=H34(c1,c2) + fac_10 * bse_mat(alpha,gama,1,m)
                               H34(c1,c2)=H34(c1,c2) + fac_11 * bse_mat(alpha,gama,1,m)
                               H34(c1,c2)=H34(c1,c2) + fac_12 * bse_mat(alpha,gama,1,m)
                           end if
                      end if
                     
                   end do
                end do
             end do
          end do
       end do
       call cpu_time(t2)
      CALL h5open_f(error)
      !CALL h5fcreate_f(filename1,  H5F_ACC_TRUNC_F, file1_id, error)


      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)

      !call h5fcreate_f(filename1, H5F_ACC_TRUNC_F, file1_id, error, &
       !           access_prp = plist_id)
      CALL h5fopen_f (filename1, H5F_ACC_RDWR_F, file1_id, error,access_prp = plist_id)


       !CALL h5fopen_f (filename2, H5F_ACC_RDWR_F, file2_id, error,access_prp = plist_id)
       CALL h5pclose_f(plist_id, error)
       dimsf(1) = nc*nv
       dimsf(2) = q
       CALL h5screate_simple_f(2, dimsf, dataspace1, error)

       CALL h5dcreate_f(file1_id, dsetname1,H5T_NATIVE_DOUBLE,dataspace1,dset1_id, error)
       CALL h5sclose_f(dataspace1, error)
     !
       !CALL h5dopen_f(file2_id, dsetname2, dset2_id, error)
     !
     ! Get  dataspace identifier.
     !
      count(1) = nrows_p
       count(2) = q
       dimsm(1)=nrows_p
       dimsm(2)=q
       dimsf(1) = nc*nv
       dimsf(2) = q
       !print*,"rank,count(1),count(2),offset(1),offset(2),nc*nv",rank,count(1),count(2),offset(1),offset(2),nc*nv

       !call h5screate_simple_f(2, dimsm, memspace, error)
       !allocate(data_out(count(1),count(2)))
       !data_out = real(H12)
       offset(1) = rank*nrows_p
       offset(2) = 0
       print*,"rank,count(1),count(2),offset(1),offset(2),nc*nv",rank,count(1),count(2),offset(1),offset(2),nc*nv

       call h5screate_simple_f(2, count, memspace, error)
       CALL h5dget_space_f(dset1_id, dataspace1, error)
       CALL h5sselect_hyperslab_f (dataspace1, H5S_SELECT_SET_F, offset,&
               count, error)


       CALL h5pcreate_f(H5P_DATASET_XFER_F, fapl_id, error)
       CALL h5pset_dxpl_mpio_f(fapl_id, H5FD_MPIO_COLLECTIVE_F, error)
       call h5screate_simple_f(2, dimsm, memspace, error)
       CALL h5dwrite_f(dset1_id, H5T_NATIVE_DOUBLE, real(H34), dimsf , error, &
                    mem_space_id=memspace, file_space_id=dataspace1)
       call h5pclose_f(fapl_id, error)
       CALL h5sclose_f(memspace, error)
       CALL h5dclose_f(dset1_id, error)

       CALL h5fclose_f(file1_id, error)

       CALL h5close_f(error)
       call cpu_time(t3)

       print*,t2-t1,t3-t2




  end subroutine create_H34
 subroutine create_H5(nv,nc,mf_v,mf_c,wcoul0,s,vol,p,n5,chunk)
      double precision :: fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, & 
                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16
 
      integer, intent(in) :: nv,nc,s,p,n5,chunk
      integer :: nb,i,j,alpha,beta,gama,eta,m,n,imatrix,keyword,c1,c2,c,x,A
      double precision            :: element
      double precision ::wcoul0,vol,a1,b1
      complex(kind=8),dimension(:,:),allocatable                            :: H5
      complex(kind=8), allocatable :: bse_mat_e(:,:,:,:)
      complex(kind=8), allocatable :: bse_mat_h(:,:,:,:)
      complex(kind=8), allocatable :: bse_mat(:,:,:,:)
      CHARACTER(LEN=11), PARAMETER :: filename1 = "singlet.h5" ! File name
      CHARACTER(LEN=5), PARAMETER :: dsetname1 = "CISD5"    ! Dataset name
      CHARACTER(LEN=11), PARAMETER :: filename2 = "singletX.h5" ! File name
      CHARACTER(LEN=5), PARAMETER :: dsetname2 = "HCISD"    ! Dataset name

      INTEGER, PARAMETER   :: DRANK = 3 ! Dataset rank 1 dimension for real or complex 2 is row-index 3 is col-index
      INTEGER(SIZE_T), PARAMETER :: NUMP = 2 ! Number of points selected
      complex(kind=8) ,DIMENSION(nv) ,intent(in) :: mf_v
      complex(kind=8),DiMENSION(nc),intent(in) :: mf_c
      INTEGER(HID_T) :: file1_id       ! File1 identifier
      INTEGER(HID_T) :: dset1_id       ! Dataset1 identifier
      INTEGER(HID_T) :: file2_id       ! File2 identifier
      INTEGER(HID_T) :: dset2_id       ! Dataset2 identifier
      INTEGER(HID_T) :: dataspace1     ! Dataspace identifier
          ! memspace identifier
      INTEGER(HID_T) :: dataspace2     ! Dataspace identifier
            ! memspace identifier
      INTEGER(HID_T) :: fapl_id

      INTEGER(HSIZE_T), DIMENSION(2) :: dimsm
                                                   ! Memory dataspace dimensions
      INTEGER(HSIZE_T), DIMENSION(2) :: dimsf
                                                   ! File dataspace dimensions
      INTEGER(HSIZE_T):: count(2),offset(2) ! Elements coordinates
                                                      ! in the file
      DOUBLE PRECISION, DIMENSION(2) :: val =(/3.14,1.414/)  ! Values to write
      INTEGER :: memrank = 1  ! Rank of the dataset in memory
      DOUBLE PRECISION,DIMENSION(:,:,:),allocatable :: bufnew
      

      INTEGER :: error,info  ! Error flag
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims

      integer(HID_T) :: filespace     ! Dataspace identifier in file
      integer(HID_T) :: memspace      ! Dataspace identifier in mem
      integer(HID_T) :: dset_id       ! Dataset identifier
      INTEGER(HID_T) :: plist_id      ! Property list identifier

      !!!!!!!!!!!mpivariables!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer :: nproc, chunk_size,itr_length
      integer , dimension(:,:) , allocatable :: ij_mat
      integer, dimension(:),allocatable :: row_ij
      integer :: comm, rank, numproc, ierror
      integer :: cij,cab,cijp,c1p
      !call MPI_INIT(ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)

       keyword = 11
       fac_hd=0
       fac_hx=0
       fac_ed=0
       fac_ex=0
      fac_1=0
      fac_2=0
      fac_3=0
      fac_4=0
      fac_5=0
      fac_6=0
      fac_7=0
      fac_8=0
      fac_9=0
      fac_10=0
      fac_11=0
      fac_12=0
      fac_13=0
      fac_14=0
      fac_15=0
      fac_16=0


      chunk_size = chunk
      allocate(bse_mat_h(nv,nv,nv,nv))
      allocate(bse_mat_e(nc,1,1,nc))
      allocate(bse_mat(nv,nv,2,nc))

      allocate(row_ij(chunk_size))
      A = nv
      allocate(ij_mat(2,chunk_size))
      c = rank+1
      nrows_p= chunk_size*A
      allocate(H5(nrows_p,p))
      H5 = cmplx(0.0,0.0)
      !CALL h5screate_simple_f(1, dimsm, memspace, error)

      !print*,"proc_id+1",c
      call distribute_ij(rank,nc,A,chunk_size,ij_mat,row_ij)
      c1=0
      c2=0
      cij=0
      cijp = 0
      cab=0
      c1p = 0
       do imatrix = 3,4
       call fac_values(wcoul0,vol,keyword,imatrix,fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16)
         if (imatrix==3) then
            call load_Hhh(nv,bse_mat_h)
         end if
         do x = 1,chunk_size
            print*,x,ij_mat(1,x),ij_mat(2,x)

            i = ij_mat(1,x)
            j = ij_mat(2,x)
            cijp = row_ij(x)
            cij = (x-1)*A
            cab =0

            if (imatrix .eq. 3) then
                call load_Hee(i,j,nc,nv,bse_mat_e)
                call load_Heh(i,j,nc,nv,imatrix,bse_mat)
            else if (imatrix .eq. 4) then
                call load_Heh(i,j,nc,nv,imatrix,bse_mat)
            end if
            do alpha = 1,nv
                  cab= cab+1
                  c1 = cij+cab
                  c1p = cijp + cab
                  c2 = 0
                  !print*,i,j,alpha,beta,c1,imatrix
                   do n = 1,nc-1
                      do m = n+1,nc
                          do gama = 1,nv
                             c2 = c2 + 1
                             !print*,"i,j,alpha,m,n,gama,c1,c2",i,j,alpha,m,n,gama,c1,c2
                             if (imatrix==3 .and. c1p==c2) then
                                   H5(c1,c2)=H5(c1,c2)+ (mf_c(j) +mf_c(i)  -(2*mf_v(alpha)))
                             end if

                             if(i==m .and. j==n) then
                                if (imatrix .ne. 4) then
                                     !print*,"bb",c1,c2, H(5,11),H(71,1),imatrix
                                      !print*,"1",H(c1,c2),bse_mat(eta,beta,alpha,gama), bse_mat(gama,beta,alpha,eta)
                                   H5(c1,c2) = H5(c1,c2) + (fac_hd * bse_mat_h(gama,alpha,alpha,gama)) &
                                       + (fac_hx * bse_mat_h(gama,alpha,alpha,gama))   !<
                                     !print*,"1", i,j,alpha,beta,m,n,gama,eta
                                     !print*,c1,c2,H(c1,c2)
                               end if
                            end if
                            if(alpha==gama) then
                               if (imatrix .ne. 4) then
                                     !print*,"2",H(c1,c2),bse_mat(n,2,1,m)
                                  H5(c1,c2)=H5(c1,c2) + (fac_ed *bse_mat_e(n,1,1,m)) &
                                       + (fac_ex * bse_mat_e(m,1,1,n))
                                    !print*,"2", i,j,alpha,beta,m,n,gama,eta
                                    !print*,c1,c2,H(c1,c2)
                              end if
                            end if
                            if(j==n .and. alpha==gama) then
                              if (imatrix .ne. 4) then
                                 H5(c1,c2)= H5(c1,c2) + fac_1 * bse_mat(alpha,gama,1,m)
                                 H5(c1,c2)= H5(c1,c2) + fac_5 * bse_mat(alpha,gama,1,m)
                                 H5(c1,c2)= H5(c1,c2) + fac_9 * bse_mat(alpha,gama,1,m)
                                 H5(c1,c2)= H5(c1,c2) + fac_13 * bse_mat(alpha,gama,1,m)
                                    !print*,"3", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                               end if
                               if (imatrix == 4) then
                                 H5(c1,c2)= H5(c1,c2) + fac_1 * bse_mat(alpha,gama,1,m)
                                 H5(c1,c2)= H5(c1,c2) + fac_5 * bse_mat(alpha,gama,1,m)
                                 H5(c1,c2)= H5(c1,c2) + fac_9 * bse_mat(alpha,gama,1,m)
                                 H5(c1,c2)= H5(c1,c2) + fac_13 * bse_mat(alpha,gama,1,m)

                              end if
                            end if
                            if(j==m .and. alpha==gama) then
                              if (imatrix .ne. 4) then
                                     !print*,"4",H(c1,c2),bse_mat(alpha,eta,1,n)
                                H5(c1,c2)=H5(c1,c2) + fac_2 * bse_mat(alpha,gama,1,n)
                                H5(c1,c2)=H5(c1,c2) + fac_6 * bse_mat(alpha,gama,1,n)
                                H5(c1,c2)=H5(c1,c2) + fac_10 * bse_mat(alpha,gama,1,n)
                                H5(c1,c2)=H5(c1,c2) + fac_14 * bse_mat(alpha,gama,1,n)

                                     !print*,"4", i,j,alpha,beta,m,n,gama,eta
                                     !print*, H(c1,c2),c1,c2
                                     !print*,"fac",fac
                              end if
                              if (imatrix == 4) then
                                H5(c1,c2)=H5(c1,c2) + fac_2 * bse_mat(alpha,gama,1,n)
                                H5(c1,c2)=H5(c1,c2) + fac_6 * bse_mat(alpha,gama,1,n)
                                H5(c1,c2)=H5(c1,c2) + fac_10 * bse_mat(alpha,gama,1,n)
                                H5(c1,c2)=H5(c1,c2) + fac_14 * bse_mat(alpha,gama,1,n)

                              end if
                            end if
                            if(i==n .and. alpha==gama) then
                              if (imatrix .ne. 4) then
                                H5(c1,c2)=H5(c1,c2) + fac_3 * bse_mat(alpha,gama,2,m)
                                H5(c1,c2)=H5(c1,c2) + fac_7 * bse_mat(alpha,gama,2,m)
                                H5(c1,c2)=H5(c1,c2) + fac_11 * bse_mat(alpha,gama,2,m)
                                H5(c1,c2)=H5(c1,c2) + fac_15 * bse_mat(alpha,gama,2,m)
                                      !print*,"5", i,j,alpha,beta,m,n,gama,eta
                                      !print*, H(c1,c2),c1,c2
                              end if
                              if (imatrix == 4) then
                                H5(c1,c2)=H5(c1,c2) + fac_3 * bse_mat(alpha,gama,2,m)
                                H5(c1,c2)=H5(c1,c2) + fac_7 * bse_mat(alpha,gama,2,m)
                                H5(c1,c2)=H5(c1,c2) + fac_11 * bse_mat(alpha,gama,2,m)
                                H5(c1,c2)=H5(c1,c2) + fac_15 * bse_mat(alpha,gama,2,m)

                              end if
                            end if
                            if(i==m .and. alpha==gama) then
                              if (imatrix .ne. 4) then
                                 H5(c1,c2)=H5(c1,c2) + fac_4 * bse_mat(alpha,gama,2,n)
                                 H5(c1,c2)=H5(c1,c2) + fac_8 * bse_mat(alpha,gama,2,n)
                                 H5(c1,c2)=H5(c1,c2) + fac_12 * bse_mat(alpha,gama,2,n)
                                 H5(c1,c2)=H5(c1,c2) + fac_16 * bse_mat(alpha,gama,2,n)
                                    !print*,"6", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                               end if
                               if (imatrix == 4) then
                                H5(c1,c2)=H5(c1,c2) + fac_4 * bse_mat(alpha,gama,2,n)
                                H5(c1,c2)=H5(c1,c2) + fac_8 * bse_mat(alpha,gama,2,n)
                                H5(c1,c2)=H5(c1,c2) + fac_12 * bse_mat(alpha,gama,2,n)
                                H5(c1,c2)=H5(c1,c2) + fac_16 * bse_mat(alpha,gama,2,n)

                               end if
                            end if
                           
                           end do
                        end do
                    end do
                 end do
             end do
          end do
           call cpu_time(t2)
      CALL h5open_f(error)
      !CALL h5fcreate_f(filename1,  H5F_ACC_TRUNC_F, file1_id, error)


      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)

      !call h5fcreate_f(filename1, H5F_ACC_TRUNC_F, file1_id, error, &
       !           access_prp = plist_id)
      CALL h5fopen_f (filename1, H5F_ACC_RDWR_F, file1_id, error,access_prp = plist_id)


       !CALL h5fopen_f (filename2, H5F_ACC_RDWR_F, file2_id, error,access_prp = plist_id)
       CALL h5pclose_f(plist_id, error)
       dimsf(1) = p
       dimsf(2) = p
       CALL h5screate_simple_f(2, dimsf, dataspace1, error)

       CALL h5dcreate_f(file1_id, dsetname1,H5T_NATIVE_DOUBLE,dataspace1,dset1_id, error)
       CALL h5sclose_f(dataspace1, error)
     !
       !CALL h5dopen_f(file2_id, dsetname2, dset2_id, error)
     !
     ! Get  dataspace identifier.
     !
      count(1) = nrows_p
       count(2) = p
       dimsm(1)=nrows_p
       dimsm(2)=p
       dimsf(1) = p
       dimsf(2) = p
       !call h5screate_simple_f(2, dimsm, memspace, error)
       !allocate(data_out(count(1),count(2)))
       !data_out = real(H12)
       offset(1) = rank*nrows_p
       offset(2) = 0
       call h5screate_simple_f(2, count, memspace, error)
       CALL h5dget_space_f(dset1_id, dataspace1, error)
       CALL h5sselect_hyperslab_f (dataspace1, H5S_SELECT_SET_F, offset,&
               count, error)


       CALL h5pcreate_f(H5P_DATASET_XFER_F, fapl_id, error)
       CALL h5pset_dxpl_mpio_f(fapl_id, H5FD_MPIO_COLLECTIVE_F, error)
       call h5screate_simple_f(2, dimsm, memspace, error)
       CALL h5dwrite_f(dset1_id, H5T_NATIVE_DOUBLE, real(H5), dimsf , error, &
                    mem_space_id=memspace, file_space_id=dataspace1)
       call h5pclose_f(fapl_id, error)

       CALL h5sclose_f(memspace, error)
       !call h5pclose_f(fapl_id, error)

       CALL h5dclose_f(dset1_id, error)
       CALL h5fclose_f(file1_id, error)

       CALL h5close_f(error)
       call cpu_time(t3)

       print*,t2-t1,t3-t2



  end subroutine create_H5
  subroutine create_H15(nv,nc,mf_v,mf_c,wcoul0,s,vol,p,n5,chunk)
      double precision :: fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16

      integer, intent(in) :: nv,nc,s,p,n5,chunk
      integer :: nb,i,j,alpha,beta,gama,eta,m,n,imatrix,keyword,c1,c2,c,x,A
      double precision            :: element
      double precision ::wcoul0,vol,a1,b1
      complex(kind=8),dimension(:,:),allocatable            :: H15
      complex(kind=8), allocatable :: bse_mat_e(:,:,:,:)
      complex(kind=8), allocatable :: bse_mat_h(:,:,:,:)
      complex(kind=8), allocatable :: bse_mat(:,:,:,:)
      CHARACTER(LEN=11), PARAMETER :: filename1 = "singlet.h5" ! File name
      CHARACTER(LEN=5), PARAMETER :: dsetname1 = "CIS15"    ! Dataset name
      CHARACTER(LEN=11), PARAMETER :: filename2 = "singletX.h5" ! File name
      CHARACTER(LEN=5), PARAMETER :: dsetname2 = "HCISD"    ! Dataset name

      INTEGER, PARAMETER   :: DRANK = 3 ! Dataset rank 1 dimension for real or complex 2 is row-index 3 is col-index
      INTEGER(SIZE_T), PARAMETER :: NUMP = 2 ! Number of points selected
      complex(kind=8) ,DIMENSION(nv) ,intent(in) :: mf_v
      complex(kind=8),DiMENSION(nc),intent(in) :: mf_c
      INTEGER(HID_T) :: file1_id       ! File1 identifier
      INTEGER(HID_T) :: dset1_id       ! Dataset1 identifier
      INTEGER(HID_T) :: file2_id       ! File2 identifier
      INTEGER(HID_T) :: dset2_id       ! Dataset2 identifier
      INTEGER(HID_T) :: dataspace1     ! Dataspace identifier
          ! memspace identifier
      INTEGER(HID_T) :: dataspace2     ! Dataspace identifier
            ! memspace identifier
      INTEGER(HID_T) :: fapl_id

      INTEGER(HSIZE_T), DIMENSION(2) :: dimsm
                                                   ! Memory dataspace dimensions
      INTEGER(HSIZE_T), DIMENSION(2) :: dimsf
                                                   ! File dataspace dimensions
      INTEGER(HSIZE_T) :: count(2),offset(2) ! Elements coordinates
                                                      ! in the file
      DOUBLE PRECISION, DIMENSION(2) :: val =(/3.14,1.414/)  ! Values to write
      INTEGER :: memrank = 1  ! Rank of the dataset in memory
      DOUBLE PRECISION,DIMENSION(:,:,:),allocatable :: bufnew

      INTEGER :: error,info  ! Error flag
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims

      integer(HID_T) :: filespace     ! Dataspace identifier in file
      integer(HID_T) :: memspace      ! Dataspace identifier in mem
      integer(HID_T) :: dset_id       ! Dataset identifier
      INTEGER(HID_T) :: plist_id      ! Property list identifier
       !!!!!!!!!!!mpivariables!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer :: nproc, chunk_size,itr_length
      integer , dimension(:,:) , allocatable :: ij_mat
      integer, dimension(:),allocatable :: row_ij
      integer :: comm, rank, numproc, ierror
      integer :: cij,cab
      !call MPI_INIT(ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)


       keyword = 12
       fac_hd=0
       fac_hx=0
       fac_ed=0
       fac_ex=0
      fac_1=0
      fac_2=0
      fac_3=0
      fac_4=0
      fac_5=0
      fac_6=0
      fac_7=0
      fac_8=0
      fac_9=0
      fac_10=0
      fac_11=0
      fac_12=0
      fac_13=0
      fac_14=0
      fac_15=0
      fac_16=0

       itr_length = nc*(nc-1)/2
       chunk_size = chunk
      print*,chunk_size,nproc

      allocate(bse_mat_h(nv,nv,nv,nv))
      allocate(bse_mat_e(nc,1,1,nc))
      allocate(bse_mat(nv,nv,2,nc))

      allocate(row_ij(chunk_size))
      A = ((nv)*(nv-1))/2
      allocate(ij_mat(2,chunk_size))
      c = rank+1
      nrows_p= chunk_size*A
      allocate(H15(nrows_p,p))
      H15 = cmplx(0.0,0.0)


      !print*,"proc_id+1",c
      call distribute_ij(rank,nc,A,chunk_size,ij_mat,row_ij)
      c1=0
      c2=0
      cij=0
      cab=0
      do imatrix = 3,4
       call fac_values(wcoul0,vol,keyword,imatrix,fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16)
         if (imatrix==3) then
            call load_Hhh(nv,bse_mat_h)
         end if
         do x = 1,chunk_size
            print*,x,ij_mat(1,x),ij_mat(2,x)

            i = ij_mat(1,x)
            j = ij_mat(2,x)
            cij = (x-1)*A
            cab =0

            if (imatrix .eq. 3) then
                call load_Hee(i,j,nc,nv,bse_mat_e)
                call load_Heh(i,j,nc,nv,imatrix,bse_mat)
            else if (imatrix .eq. 4) then
                call load_Heh(i,j,nc,nv,imatrix,bse_mat)
            end if
            do beta = 1,nv-1
               do alpha = beta+1,nv
                  cab = cab +1
                  c1 = cij+cab
                  c2 = 0
                   do n = 1,nc-1
                      do m = n+1,nc
                         do gama = 1,nv
                             c2 = c2 +1
                              !print*,"i,j,alpha,beta,m,n,gama,c1,c2",i,j,alpha,beta,m,n,gama,c1,c2

                             !if(i==m .and. j==n) then
                              !  if (imatrix .ne. 4) then
                                     !print*,"bb",c1,c2, H(5,11),H(71,1),imatrix
                                      !print*,"1",H(c1,c2),bse_mat(eta,beta,alpha,gama), bse_mat(gama,beta,alpha,eta)
                               !    H15(c1,c2) = H15(c1,c2) + (fac_hd * bse_mat_h(gama,beta,alpha,gama)) &
                                !       + (fac_hx * bse_mat_h(gama,beta,alpha,gama))   !<
                                     !print*,"1", i,j,alpha,beta,m,n,gama,eta
                                     !print*,c1,c2,H15(c1,c2),"12"
                                !end if
                             !end if
                             if(j==n .and. beta==gama) then
                               if (imatrix .ne. 4) then
                                   H15(c1,c2)= H15(c1,c2) + fac_1 * bse_mat(alpha,gama,1,m)
                                   H15(c1,c2)=H15(c1,c2) +fac_9 * bse_mat(alpha,gama,1,m)
                                    !print*,"3", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                               end if
                               if (imatrix == 4) then
                                   H15(c1,c2)= H15(c1,c2) + fac_1 * bse_mat(alpha,gama,1,m)
                                   H15(c1,c2)=H15(c1,c2) +fac_9 * bse_mat(alpha,gama,1,m)
                                   !print*,c1,c2,H15(c1,c2), bse_mat(alpha,gama,1,m),fac_9
                               end if
                            end if
                            if(j==m .and. beta==gama) then
                               if (imatrix .ne. 4) then
                                     !print*,"4",H(c1,c2),bse_mat(alpha,eta,1,n)
                                   H15(c1,c2)=H15(c1,c2) + fac_2 * bse_mat(alpha,gama,1,n)
                                   H15(c1,c2)=H15(c1,c2) + fac_10 * bse_mat(alpha,gama,1,n)
                                     !print*,"4", i,j,alpha,beta,m,n,gama,eta
                                     !print*, H(c1,c2),c1,c2
                                     !print*,"fac",fac
                               end if
                               if (imatrix == 4) then
                                   H15(c1,c2)=H15(c1,c2) + fac_2 * bse_mat(alpha,gama,1,n)
                                   H15(c1,c2)=H15(c1,c2) + fac_10 * bse_mat(alpha,gama,1,n)
                               end if
                            end if
                            if(i==n .and. beta==gama) then
                               if (imatrix .ne. 4) then
                                   H15(c1,c2)=H15(c1,c2) + fac_3 * bse_mat(alpha,gama,2,m)
                                   H15(c1,c2)=H15(c1,c2) + fac_11 * bse_mat(alpha,gama,2,m)
                                      !print*,"5", i,j,alpha,beta,m,n,gama,eta
                                      !print*, H(c1,c2),c1,c2
                               end if
                               if (imatrix == 4) then
                                   H15(c1,c2)=H15(c1,c2) + fac_3 * bse_mat(alpha,gama,2,m)
                                   H15(c1,c2)=H15(c1,c2) + fac_11 * bse_mat(alpha,gama,2,m)
                               end if
                            end if
                             if(i==m .and. beta==gama) then
                               if (imatrix .ne. 4) then
                                  H15(c1,c2)=H15(c1,c2) + fac_4 * bse_mat(alpha,gama,2,n)
                                  H15(c1,c2)=H15(c1,c2) + fac_12 * bse_mat(alpha,gama,2,n)
                                    !print*,"6", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                               end if
                               if (imatrix == 4) then
                                  H15(c1,c2)=H15(c1,c2) + fac_4 * bse_mat(alpha,gama,2,n)
                                  H15(c1,c2)=H15(c1,c2) + fac_12 * bse_mat(alpha,gama,2,n)
                               end if
                            end if
                            if(j==n .and. alpha==gama) then
                               if (imatrix .ne. 4) then
                                  H15(c1,c2)=H15(c1,c2) + fac_5 * bse_mat(beta,gama,1,m)
                                  H15(c1,c2)=H15(c1,c2) + fac_13 * bse_mat(beta,gama,1,m)
                                    !print*, H(c1,c2),c1,c2
                               end if
                               if (imatrix == 4) then
                                  H15(c1,c2)=H15(c1,c2) + fac_5 * bse_mat(beta,gama,1,m)
                                  H15(c1,c2)=H15(c1,c2) + fac_13 * bse_mat(beta,gama,1,m)
                               end if
                            end if
                            if(j==m .and. alpha==gama) then
                               if (imatrix .ne. 4) then
                                  H15(c1,c2)=H15(c1,c2) + fac_6 * bse_mat(beta,gama,1,n)
                                  H15(c1,c2)=H15(c1,c2) + fac_14 * bse_mat(beta,gama,1,n)
                                    !print*,"8", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                               end if
                               if (imatrix == 4) then
                                  H15(c1,c2)=H15(c1,c2) + fac_6 * bse_mat(beta,gama,1,n)
                                  H15(c1,c2)=H15(c1,c2) + fac_14 * bse_mat(beta,gama,1,n)
                               end if
                            end if
                            if(i==n .and. alpha==gama) then
                               if (imatrix .ne. 4) then
                                  H15(c1,c2)=H15(c1,c2) + fac_7 * bse_mat(beta,gama,2,m)
                                  H15(c1,c2)=H15(c1,c2) + fac_15 * bse_mat(beta,gama,2,m)
                                    !print*,"9", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                               end if
                               if (imatrix == 4) then
                                  H15(c1,c2)=H15(c1,c2) + fac_7 * bse_mat(beta,gama,2,m)
                                  H15(c1,c2)=H15(c1,c2) + fac_15 * bse_mat(beta,gama,2,m)
                               end if
                            end if
                            if(i==m .and. alpha==gama) then
                               if (imatrix .ne. 4) then
                                  H15(c1,c2)=H15(c1,c2) + fac_8 * bse_mat(beta,gama,2,n)
                                  H15(c1,c2)=H15(c1,c2) + fac_16 * bse_mat(beta,gama,2,n)
                                    !print*,"10", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                               end if
                               if (imatrix == 4) then
                                  H15(c1,c2)=H15(c1,c2) + fac_8 * bse_mat(beta,gama,2,n)
                                  H15(c1,c2)=H15(c1,c2) + fac_16 * bse_mat(beta,gama,2,n)
                               end if
                            end if

                         end do
                      end do
                   end do
               end do
            end do
         end do
      end do
       call cpu_time(t2)
      CALL h5open_f(error)
      !CALL h5fcreate_f(filename1,  H5F_ACC_TRUNC_F, file1_id, error)


      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)

      !call h5fcreate_f(filename1, H5F_ACC_TRUNC_F, file1_id, error, &
       !           access_prp = plist_id)
      CALL h5fopen_f (filename1, H5F_ACC_RDWR_F, file1_id, error,access_prp = plist_id)


       !CALL h5fopen_f (filename2, H5F_ACC_RDWR_F, file2_id, error,access_prp = plist_id)
       CALL h5pclose_f(plist_id, error)
       dimsf(1) = s
       dimsf(2) = p
       CALL h5screate_simple_f(2, dimsf, dataspace1, error)

       CALL h5dcreate_f(file1_id, dsetname1,H5T_NATIVE_DOUBLE,dataspace1,dset1_id, error)
       CALL h5sclose_f(dataspace1, error)
     !
       !CALL h5dopen_f(file2_id, dsetname2, dset2_id, error)
     !
     ! Get  dataspace identifier.
     !
      count(1) = nrows_p
       count(2) = p
       dimsm(1)=nrows_p
       dimsm(2)=p
       dimsf(1) = s
       dimsf(2) = p
       !call h5screate_simple_f(2, dimsm, memspace, error)
       !allocate(data_out(count(1),count(2)))
       !data_out = real(H12)
       offset(1) = rank*nrows_p
       offset(2) = 0
       call h5screate_simple_f(2, count, memspace, error)
       CALL h5dget_space_f(dset1_id, dataspace1, error)
       CALL h5sselect_hyperslab_f (dataspace1, H5S_SELECT_SET_F, offset,&
               count, error)


       CALL h5pcreate_f(H5P_DATASET_XFER_F, fapl_id, error)
       CALL h5pset_dxpl_mpio_f(fapl_id, H5FD_MPIO_COLLECTIVE_F, error)
       call h5screate_simple_f(2, dimsm, memspace, error)
       CALL h5dwrite_f(dset1_id, H5T_NATIVE_DOUBLE, real(H15), dimsf , error, &
                    mem_space_id=memspace, file_space_id=dataspace1)
       call h5pclose_f(fapl_id, error)

       CALL h5sclose_f(memspace, error)
        CALL h5dclose_f(dset1_id, error)

       CALL h5fclose_f(file1_id, error)

       CALL h5close_f(error)
       call cpu_time(t3)

       print*,t2-t1,t3-t2


  end subroutine create_H15
  subroutine create_H25(nv,nc,mf_v,mf_c,wcoul0,s,vol,p,n5,chunk)
      double precision :: fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16

      integer, intent(in) :: nv,nc,s,p,n5,chunk
      integer :: nb,i,j,alpha,beta,gama,eta,m,n,imatrix,keyword,c1,c2,c,x,A
      double precision            :: element
      double precision ::wcoul0,vol,a1,b1
      complex(kind=8),dimension(:,:),allocatable            :: H25
      complex(kind=8), allocatable :: bse_mat_e(:,:,:,:)
      complex(kind=8), allocatable :: bse_mat_h(:,:,:,:)
      complex(kind=8), allocatable :: bse_mat(:,:,:,:)
      CHARACTER(LEN=11), PARAMETER :: filename1 = "singlet.h5" ! File name
      CHARACTER(LEN=5), PARAMETER :: dsetname1 = "CIS25"    ! Dataset name
      CHARACTER(LEN=11), PARAMETER :: filename2 = "singletX.h5" ! File name
      CHARACTER(LEN=5), PARAMETER :: dsetname2 = "HCISD"    ! Dataset name

      INTEGER, PARAMETER   :: DRANK = 3 ! Dataset rank 1 dimension for real or complex 2 is row-index 3 is col-index
      INTEGER(SIZE_T), PARAMETER :: NUMP = 2 ! Number of points selected
      complex(kind=8) ,DIMENSION(nv) ,intent(in) :: mf_v
      complex(kind=8),DiMENSION(nc),intent(in) :: mf_c
      INTEGER(HID_T) :: file1_id       ! File1 identifier
      INTEGER(HID_T) :: dset1_id       ! Dataset1 identifier
      INTEGER(HID_T) :: file2_id       ! File2 identifier
      INTEGER(HID_T) :: dset2_id       ! Dataset2 identifier
      INTEGER(HID_T) :: dataspace1     ! Dataspace identifier
          ! memspace identifier
      INTEGER(HID_T) :: dataspace2     ! Dataspace identifier
            ! memspace identifier
      INTEGER(HID_T) :: fapl_id

      INTEGER(HSIZE_T), DIMENSION(2) :: dimsm
                                                   ! Memory dataspace dimensions
      INTEGER(HSIZE_T), DIMENSION(2) :: dimsf
                                                   ! File dataspace dimensions
      INTEGER(HSIZE_T) :: count(2),offset(2) ! Elements coordinates
                                                      ! in the file
      DOUBLE PRECISION, DIMENSION(2) :: val =(/3.14,1.414/)  ! Values to write
      INTEGER :: memrank = 1  ! Rank of the dataset in memory
      DOUBLE PRECISION,DIMENSION(:,:,:),allocatable :: bufnew

      INTEGER :: error,info  ! Error flag
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims

      integer(HID_T) :: filespace     ! Dataspace identifier in file
      integer(HID_T) :: memspace      ! Dataspace identifier in mem
      integer(HID_T) :: dset_id       ! Dataset identifier
      INTEGER(HID_T) :: plist_id      ! Property list identifier
       !!!!!!!!!!!mpivariables!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer :: nproc, chunk_size,itr_length
      integer , dimension(:,:) , allocatable :: ij_mat
      integer, dimension(:),allocatable :: row_ij
      integer :: comm, rank, numproc, ierror
      integer :: cij,cab
      !call MPI_INIT(ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)


       keyword = 13
       fac_hd=0
       fac_hx=0
       fac_ed=0
       fac_ex=0
      fac_1=0
      fac_2=0
      fac_3=0
      fac_4=0
      fac_5=0
      fac_6=0
      fac_7=0
      fac_8=0
      fac_9=0
      fac_10=0
      fac_11=0
      fac_12=0
      fac_13=0
      fac_14=0
      fac_15=0
      fac_16=0

       itr_length = nc*(nc-1)/2
       chunk_size = chunk
      print*,chunk_size,nproc

      allocate(bse_mat_h(nv,nv,nv,nv))
      allocate(bse_mat_e(nc,1,1,nc))
      allocate(bse_mat(nv,nv,2,nc))

      allocate(row_ij(chunk_size))
      A = ((nv)*(nv-1))/2
      allocate(ij_mat(2,chunk_size))
      c = rank+1
      nrows_p= chunk_size*A
      allocate(H25(nrows_p,p))
      H25 = cmplx(0.0,0.0)


      !print*,"proc_id+1",c
      call distribute_ij(rank,nc,A,chunk_size,ij_mat,row_ij)
      c1=0
      c2=0
      cij=0
      cab=0
      do imatrix = 3,4
       call fac_values(wcoul0,vol,keyword,imatrix,fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16)
         if (imatrix==3) then
            call load_Hhh(nv,bse_mat_h)
         end if
         do x = 1,chunk_size
            print*,x,ij_mat(1,x),ij_mat(2,x)

            i = ij_mat(1,x)
            j = ij_mat(2,x)
            cij = (x-1)*A
            cab =0

            if (imatrix .eq. 3) then
                call load_Hee(i,j,nc,nv,bse_mat_e)
                call load_Heh(i,j,nc,nv,imatrix,bse_mat)
            else if (imatrix .eq. 4) then
                call load_Heh(i,j,nc,nv,imatrix,bse_mat)
            end if
            do beta = 1,nv-1
               do alpha = beta+1,nv
                  cab = cab +1
                  c1 = cij+cab
                  c2 = 0
                   do n = 1,nc-1
                      do m = n+1,nc
                         do gama = 1,nv
                             c2 = c2 +1
                              !print*,"i,j,alpha,beta,m,n,gama,c1,c2",i,j,alpha,beta,m,n,gama,c1,c2

                             !if(i==m .and. j==n) then
                              !  if (imatrix .ne. 4) then
                                     !print*,"bb",c1,c2, H(5,11),H(71,1),imatrix
                                      !print*,"1",H(c1,c2),bse_mat(eta,beta,alpha,gama), bse_mat(gama,beta,alpha,eta)
                               !    H25(c1,c2) = H25(c1,c2) + (fac_hd * bse_mat_h(gama,beta,alpha,gama)) &
                                !       + (fac_hx * bse_mat_h(gama,beta,alpha,gama))   !<
                                     !print*,"1", i,j,alpha,beta,m,n,gama,eta
                                     !print*,c1,c2,H25(c1,c2),"12"
                                !end if
                             !end if
                             if(j==n .and. beta==gama) then
                               if (imatrix .ne. 4) then
                                   H25(c1,c2)= H25(c1,c2) + fac_1 * bse_mat(alpha,gama,1,m)
                                   H25(c1,c2)=H25(c1,c2) +fac_9 * bse_mat(alpha,gama,1,m)
                                    !print*,"3", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                               end if
                               if (imatrix == 4) then
                                   H25(c1,c2)= H25(c1,c2) + fac_1 * bse_mat(alpha,gama,1,m)
                                   H25(c1,c2)=H25(c1,c2) +fac_9 * bse_mat(alpha,gama,1,m)
                                   !print*,c1,c2,H25(c1,c2), bse_mat(alpha,gama,1,m),fac_9
                               end if
                            end if
                            if(j==m .and. beta==gama) then
                               if (imatrix .ne. 4) then
                                     !print*,"4",H(c1,c2),bse_mat(alpha,eta,1,n)
                                   H25(c1,c2)=H25(c1,c2) + fac_2 * bse_mat(alpha,gama,1,n)
                                   H25(c1,c2)=H25(c1,c2) + fac_10 * bse_mat(alpha,gama,1,n)
                                     !print*,"4", i,j,alpha,beta,m,n,gama,eta
                                     !print*, H(c1,c2),c1,c2
                                     !print*,"fac",fac
                               end if
                               if (imatrix == 4) then
                                   H25(c1,c2)=H25(c1,c2) + fac_2 * bse_mat(alpha,gama,1,n)
                                   H25(c1,c2)=H25(c1,c2) + fac_10 * bse_mat(alpha,gama,1,n)
                               end if
                            end if
                            if(i==n .and. beta==gama) then
                               if (imatrix .ne. 4) then
                                   H25(c1,c2)=H25(c1,c2) + fac_3 * bse_mat(alpha,gama,2,m)
                                   H25(c1,c2)=H25(c1,c2) + fac_11 * bse_mat(alpha,gama,2,m)
                                      !print*,"5", i,j,alpha,beta,m,n,gama,eta
                                      !print*, H(c1,c2),c1,c2
                               end if
                               if (imatrix == 4) then
                                   H25(c1,c2)=H25(c1,c2) + fac_3 * bse_mat(alpha,gama,2,m)
                                   H25(c1,c2)=H25(c1,c2) + fac_11 * bse_mat(alpha,gama,2,m)
                               end if
                            end if
                             if(i==m .and. beta==gama) then
                               if (imatrix .ne. 4) then
                                  H25(c1,c2)=H25(c1,c2) + fac_4 * bse_mat(alpha,gama,2,n)
                                  H25(c1,c2)=H25(c1,c2) + fac_12 * bse_mat(alpha,gama,2,n)
                                    !print*,"6", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                               end if
                               if (imatrix == 4) then
                                  H25(c1,c2)=H25(c1,c2) + fac_4 * bse_mat(alpha,gama,2,n)
                                  H25(c1,c2)=H25(c1,c2) + fac_12 * bse_mat(alpha,gama,2,n)
                               end if
                            end if
                            if(j==n .and. alpha==gama) then
                               if (imatrix .ne. 4) then
                                  H25(c1,c2)=H25(c1,c2) + fac_5 * bse_mat(beta,gama,1,m)
                                  H25(c1,c2)=H25(c1,c2) + fac_13 * bse_mat(beta,gama,1,m)
                                    !print*, H(c1,c2),c1,c2
                               end if
                               if (imatrix == 4) then
                                  H25(c1,c2)=H25(c1,c2) + fac_5 * bse_mat(beta,gama,1,m)
                                  H25(c1,c2)=H25(c1,c2) + fac_13 * bse_mat(beta,gama,1,m)
                               end if
                            end if
                            if(j==m .and. alpha==gama) then
                               if (imatrix .ne. 4) then
                                  H25(c1,c2)=H25(c1,c2) + fac_6 * bse_mat(beta,gama,1,n)
                                  H25(c1,c2)=H25(c1,c2) + fac_14 * bse_mat(beta,gama,1,n)
                                    !print*,"8", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                               end if
                               if (imatrix == 4) then
                                  H25(c1,c2)=H25(c1,c2) + fac_6 * bse_mat(beta,gama,1,n)
                                  H25(c1,c2)=H25(c1,c2) + fac_14 * bse_mat(beta,gama,1,n)
                               end if
                            end if
                            if(i==n .and. alpha==gama) then
                               if (imatrix .ne. 4) then
                                  H25(c1,c2)=H25(c1,c2) + fac_7 * bse_mat(beta,gama,2,m)
                                  H25(c1,c2)=H25(c1,c2) + fac_15 * bse_mat(beta,gama,2,m)
                                    !print*,"9", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                               end if
                               if (imatrix == 4) then
                                  H25(c1,c2)=H25(c1,c2) + fac_7 * bse_mat(beta,gama,2,m)
                                  H25(c1,c2)=H25(c1,c2) + fac_15 * bse_mat(beta,gama,2,m)
                               end if
                            end if
                            if(i==m .and. alpha==gama) then
                               if (imatrix .ne. 4) then
                                  H25(c1,c2)=H25(c1,c2) + fac_8 * bse_mat(beta,gama,2,n)
                                  H25(c1,c2)=H25(c1,c2) + fac_16 * bse_mat(beta,gama,2,n)
                                    !print*,"10", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                               end if
                               if (imatrix == 4) then
                                  H25(c1,c2)=H25(c1,c2) + fac_8 * bse_mat(beta,gama,2,n)
                                  H25(c1,c2)=H25(c1,c2) + fac_16 * bse_mat(beta,gama,2,n)
                               end if
                            end if

                         end do
                      end do
                   end do
               end do
            end do
         end do
      end do
       call cpu_time(t2)
      CALL h5open_f(error)
      !CALL h5fcreate_f(filename1,  H5F_ACC_TRUNC_F, file1_id, error)


      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)

      !call h5fcreate_f(filename1, H5F_ACC_TRUNC_F, file1_id, error, &
       !           access_prp = plist_id)
      CALL h5fopen_f (filename1, H5F_ACC_RDWR_F, file1_id, error,access_prp = plist_id)


       !CALL h5fopen_f (filename2, H5F_ACC_RDWR_F, file2_id, error,access_prp = plist_id)
       CALL h5pclose_f(plist_id, error)
       dimsf(1) = s
       dimsf(2) = p
       CALL h5screate_simple_f(2, dimsf, dataspace1, error)

       CALL h5dcreate_f(file1_id, dsetname1,H5T_NATIVE_DOUBLE,dataspace1,dset1_id, error)
       CALL h5sclose_f(dataspace1, error)
     !
       !CALL h5dopen_f(file2_id, dsetname2, dset2_id, error)
     !
     ! Get  dataspace identifier.
     !
      count(1) = nrows_p
       count(2) = p
       dimsm(1)=nrows_p
       dimsm(2)=p
       dimsf(1) = s
       dimsf(2) = p
       !call h5screate_simple_f(2, dimsm, memspace, error)
       !allocate(data_out(count(1),count(2)))
       !data_out = real(H12)
       offset(1) = rank*nrows_p
       offset(2) = 0
       call h5screate_simple_f(2, count, memspace, error)
       CALL h5dget_space_f(dset1_id, dataspace1, error)
       CALL h5sselect_hyperslab_f (dataspace1, H5S_SELECT_SET_F, offset,&
               count, error)


       CALL h5pcreate_f(H5P_DATASET_XFER_F, fapl_id, error)
       CALL h5pset_dxpl_mpio_f(fapl_id, H5FD_MPIO_COLLECTIVE_F, error)
       call h5screate_simple_f(2, dimsm, memspace, error)
       CALL h5dwrite_f(dset1_id, H5T_NATIVE_DOUBLE, real(H25), dimsf , error, &
                    mem_space_id=memspace, file_space_id=dataspace1)
       call h5pclose_f(fapl_id, error)

       CALL h5sclose_f(memspace, error)
        CALL h5dclose_f(dset1_id, error)

       CALL h5fclose_f(file1_id, error)

       CALL h5close_f(error)
       call cpu_time(t3)

       print*,t2-t1,t3-t2


  end subroutine create_H25
   subroutine create_H35(nv,nc,mf_v,mf_c,wcoul0,s,vol,p,n5)
      double precision :: fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16

      integer, intent(in) :: nv,nc,s,p,n5
      integer :: nb,i,j,alpha,beta,gama,eta,m,n,imatrix,keyword,c1,c2,c,x,A
      double precision            :: element
      double precision ::wcoul0,vol,a1,b1
      complex(kind=8),dimension(:,:),allocatable                            :: H35
      complex(kind=8), allocatable :: bse_mat_e(:,:,:,:)
      complex(kind=8), allocatable :: bse_mat_h(:,:,:,:)
      complex(kind=8), allocatable :: bse_mat(:,:,:,:)
      CHARACTER(LEN=11), PARAMETER :: filename1 = "singlet.h5" ! File name
      CHARACTER(LEN=5), PARAMETER :: dsetname1 = "CIS35"    ! Dataset name
      CHARACTER(LEN=11), PARAMETER :: filename2 = "singletX.h5" ! File name
      CHARACTER(LEN=5), PARAMETER :: dsetname2 = "HCISD"    ! Dataset name
      INTEGER, PARAMETER   :: DRANK = 3 ! Dataset rank 1 dimension for real or complex 2 is row-index 3 is col-index
      INTEGER(SIZE_T), PARAMETER :: NUMP = 2 ! Number of points selected
      complex(kind=8) ,DIMENSION(nv) ,intent(in) :: mf_v
      complex(kind=8),DiMENSION(nc),intent(in) :: mf_c
      INTEGER(HID_T) :: file1_id       ! File1 identifier
      INTEGER(HID_T) :: dset1_id       ! Dataset1 identifier
      INTEGER(HID_T) :: file2_id       ! File2 identifier
      INTEGER(HID_T) :: dset2_id       ! Dataset2 identifier
      INTEGER(HID_T) :: dataspace1     ! Dataspace identifier
          ! memspace identifier
      INTEGER(HID_T) :: dataspace2     ! Dataspace identifier
            ! memspace identifier
      INTEGER(HID_T) :: fapl_id

      INTEGER(HSIZE_T), DIMENSION(2) :: dimsm
                                                   ! Memory dataspace dimensions
      INTEGER(HSIZE_T), DIMENSION(2) :: dimsf
                                                   ! File dataspace dimensions
      INTEGER(HSIZE_T) :: offset(2),count(2) ! Elements coordinates
                                                      ! in the file
      DOUBLE PRECISION, DIMENSION(2) :: val =(/3.14,1.414/)  ! Values to write
      INTEGER :: memrank = 1  ! Rank of the dataset in memory
      DOUBLE PRECISION,DIMENSION(:,:,:),allocatable :: bufnew

      INTEGER :: error,info  ! Error flag
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims

      integer(HID_T) :: filespace     ! Dataspace identifier in file
      integer(HID_T) :: memspace      ! Dataspace identifier in mem
      integer(HID_T) :: dset_id       ! Dataset identifier
      INTEGER(HID_T) :: plist_id      ! Property list identifier
       !!!!!!!!!!!mpivariables!!!!!!!!!!!!!!!!!!!!!!!!!!!
       integer :: nproc, chunk_size,itr_length
      integer , dimension(:,:) , allocatable :: ij_mat
      integer, dimension(:),allocatable :: row_ij
      integer :: comm, rank, numproc, ierror
      integer :: cij,cab
      !call MPI_INIT(ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)

     

       keyword = 14
         fac_hd=0
       fac_hx=0
       fac_ed=0
       fac_ex=0
      fac_1=0
      fac_2=0
      fac_3=0
      fac_4=0
      fac_5=0
      fac_6=0
      fac_7=0
      fac_8=0
      fac_9=0
      fac_10=0
      fac_11=0
      fac_12=0
      fac_13=0
      fac_14=0
      fac_15=0
      fac_16=0
         itr_length = nc*(nc-1)/2
      if (mod(nc,2)==0) then
         nproc = nc-1
         chunk_size = nc/2
      else
         nproc = nc
         chunk_size = (nc-1)/2
      end if
      print*,chunk_size,nproc

      allocate(bse_mat_h(nv,nv,nv,nv))
      allocate(bse_mat_e(nc,1,1,nc))
      allocate(bse_mat(nv,nv,2,nc))
      !allocate(row_ij(chunk_size))
      A = nv
      !allocate(ij_mat(2,chunk_size))
      c = rank+1
      nrows_p= A
      allocate(H35(nrows_p,p))
      H35 = cmplx(0.0,0.0)
   

      c1=0
      c2=0
      cij=0
      cab=0
      do imatrix = 3,4
         call fac_values(wcoul0,vol,keyword,imatrix,fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16)

         i = c
         cij = 0
         cab = 0
         if (imatrix==3) then
            call load_Hhh(nv,bse_mat_h)
         end if

         if (imatrix .eq. 3) then
                call load_Hee(i,i,nc,nv,bse_mat_e)
                call load_Heh(i,i,nc,nv,imatrix,bse_mat)
         else if (imatrix .eq. 4) then
                call load_Heh(i,i,nc,nv,imatrix,bse_mat)
         end if

         do alpha = 1,nv
            cab = cab + 1
            c1 = cij + cab
            c2 = 0
            do n = 1,nc-1
               do m = n+1,nc
                   do gama = 1,nv
                      c2 = c2 + 1
                       if(alpha==gama) then
                        if (imatrix .ne. 4) then
                                     !print*,"2",H(c1,c2),bse_mat(n,j+nv,1,m)
                           H35(c1,c2)=H35(c1,c2) + (fac_ed *bse_mat_e(n,1,1,m)) &
                                       + (fac_ex * bse_mat_e(m,1,1,n))
                                    !print*,"2", i,j,alpha,beta,m,n,gama,eta
                                    !print*,c1,c2,H(c1,c2)
                        end if
                     end if
                     if(i==n .and. alpha==gama) then
                        if (imatrix .ne. 4) then
                            H35(c1,c2)= H35(c1,c2) + fac_1 * bse_mat(alpha,gama,1,m)
                            H35(c1,c2)=H35(c1,c2) + fac_3 * bse_mat(alpha,gama,1,m)
                            H35(c1,c2)=H35(c1,c2) + fac_5 * bse_mat(alpha,gama,1,m)
                            H35(c1,c2)=H35(c1,c2) + fac_7 * bse_mat(alpha,gama,1,m)
                            H35(c1,c2)=H35(c1,c2) + fac_9 * bse_mat(alpha,gama,1,m)
                            H35(c1,c2)=H35(c1,c2) + fac_11 * bse_mat(alpha,gama,1,m)
                            H35(c1,c2)=H35(c1,c2) + fac_13 * bse_mat(alpha,gama,1,m)
                            H35(c1,c2)=H35(c1,c2) + fac_15 * bse_mat(alpha,gama,1,m)
                             !       print*,"imatrix, i,alpha,m,n,gama,bse_mat(alpha,gama,1,m),c1,c2",&
                              !      imatrix, i,alpha,m,n,gama,bse_mat(alpha,gama,1,m),c1,c2
                                    !print*, H(c1,c2),c1,c2
                        end if
                        if (imatrix == 4) then
                            H35(c1,c2)= H35(c1,c2) + fac_1 * bse_mat(alpha,gama,1,m)
                            H35(c1,c2)=H35(c1,c2) + fac_3 * bse_mat(alpha,gama,1,m)
                            H35(c1,c2)=H35(c1,c2) + fac_5 * bse_mat(alpha,gama,1,m)
                            H35(c1,c2)=H35(c1,c2) + fac_7 * bse_mat(alpha,gama,1,m)
                            H35(c1,c2)=H35(c1,c2) + fac_9 * bse_mat(alpha,gama,1,m)
                            H35(c1,c2)=H35(c1,c2) + fac_11 * bse_mat(alpha,gama,1,m)
                            H35(c1,c2)=H35(c1,c2) + fac_13 * bse_mat(alpha,gama,1,m)
                            H35(c1,c2)=H35(c1,c2) + fac_15 * bse_mat(alpha,gama,1,m)
                            !print*,"imatrix, i,alpha,m,n,gama,bse_mat(alpha,gama,1,m),c1,c2",&
                            !        imatrix, i,alpha,m,n,gama,bse_mat(alpha,gama,1,m),c1,c2
                        end if
                     end if
                     if(i==m .and. alpha==gama) then
                        if (imatrix .ne. 4) then
                                     !print*,"4",H(c1,c2),bse_mat(alpha,eta,1,n)
                            H35(c1,c2)=H35(c1,c2) + fac_2 * bse_mat(alpha,gama,1,n)
                            H35(c1,c2)=H35(c1,c2) + fac_4 * bse_mat(alpha,gama,1,n)
                            H35(c1,c2)=H35(c1,c2) + fac_6 * bse_mat(alpha,gama,1,n)
                            H35(c1,c2)=H35(c1,c2) + fac_8 * bse_mat(alpha,gama,1,n)
                            H35(c1,c2)=H35(c1,c2) + fac_10 * bse_mat(alpha,gama,1,n)
                            H35(c1,c2)=H35(c1,c2) + fac_12 * bse_mat(alpha,gama,1,n)
                            H35(c1,c2)=H35(c1,c2) + fac_14 * bse_mat(alpha,gama,1,n)
                            H35(c1,c2)=H35(c1,c2) + fac_16 * bse_mat(alpha,gama,1,n)
                                     !print*,"4", i,j,alpha,beta,m,n,gama,eta
                                     !print*, H(c1,c2),c1,c2
                                     !print*,"fac",fac
                           !print*,"imatrix, i,alpha,m,n,gama,bse_mat(alpha,gama,1,m),c1,c2",&
                           !         imatrix, i,alpha,m,n,gama,bse_mat(alpha,gama,1,n),c1,c2
                        end if
                        if (imatrix == 4) then
                            H35(c1,c2)=H35(c1,c2) + fac_2 * bse_mat(alpha,gama,1,n)
                            H35(c1,c2)=H35(c1,c2) + fac_4 * bse_mat(alpha,gama,1,n)
                            H35(c1,c2)=H35(c1,c2) + fac_6 * bse_mat(alpha,gama,1,n)
                            H35(c1,c2)=H35(c1,c2) + fac_8 * bse_mat(alpha,gama,1,n)
                            H35(c1,c2)=H35(c1,c2) + fac_10 * bse_mat(alpha,gama,1,n)
                            H35(c1,c2)=H35(c1,c2) + fac_12 * bse_mat(alpha,gama,1,n)
                            H35(c1,c2)=H35(c1,c2) + fac_14 * bse_mat(alpha,gama,1,n)
                            H35(c1,c2)=H35(c1,c2) + fac_16 * bse_mat(alpha,gama,1,n)
                            !print*,"imatrix, i,alpha,m,n,gama,bse_mat(alpha,gama,1,m),c1,c2",&
                             !       imatrix, i,alpha,m,n,gama,bse_mat(alpha,gama,1,n),c1,c2
                        end if
                     end if
                     
                   end do
                end do
             end do
          end do
       end do
        call cpu_time(t2)
      CALL h5open_f(error)
      !CALL h5fcreate_f(filename1,  H5F_ACC_TRUNC_F, file1_id, error)


      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)

      !call h5fcreate_f(filename1, H5F_ACC_TRUNC_F, file1_id, error, &
       !           access_prp = plist_id)
      CALL h5fopen_f (filename1, H5F_ACC_RDWR_F, file1_id, error,access_prp = plist_id)


       !CALL h5fopen_f (filename2, H5F_ACC_RDWR_F, file2_id, error,access_prp = plist_id)
       CALL h5pclose_f(plist_id, error)
       dimsf(1) = nc*nv
       dimsf(2) = p
       CALL h5screate_simple_f(2, dimsf, dataspace1, error)

       CALL h5dcreate_f(file1_id, dsetname1,H5T_NATIVE_DOUBLE,dataspace1,dset1_id, error)
       CALL h5sclose_f(dataspace1, error)
     !
       !CALL h5dopen_f(file2_id, dsetname2, dset2_id, error)
     !
     ! Get  dataspace identifier.
     !
      count(1) = nrows_p
       count(2) = p
       dimsm(1)=nrows_p
       dimsm(2)=p
       dimsf(1) = nc*nv
       dimsf(2) = p
       !print*,"rank,count(1),count(2),offset(1),offset(2),nc*nv",rank,count(1),count(2),offset(1),offset(2),nc*nv

       !call h5screate_simple_f(2, dimsm, memspace, error)
       !allocate(data_out(count(1),count(2)))
       !data_out = real(H12)
       offset(1) = rank*nrows_p
       offset(2) = 0
       print*,"rank,count(1),count(2),offset(1),offset(2),nc*nv",rank,count(1),count(2),offset(1),offset(2),nc*nv

       call h5screate_simple_f(2, count, memspace, error)
       CALL h5dget_space_f(dset1_id, dataspace1, error)
       CALL h5sselect_hyperslab_f (dataspace1, H5S_SELECT_SET_F, offset,&
               count, error)


       CALL h5pcreate_f(H5P_DATASET_XFER_F, fapl_id, error)
       CALL h5pset_dxpl_mpio_f(fapl_id, H5FD_MPIO_COLLECTIVE_F, error)
       call h5screate_simple_f(2, dimsm, memspace, error)
       CALL h5dwrite_f(dset1_id, H5T_NATIVE_DOUBLE, real(H35), dimsf , error, &
                    mem_space_id=memspace, file_space_id=dataspace1)
       call h5pclose_f(fapl_id, error)
        CALL h5sclose_f(memspace, error)
       CALL h5dclose_f(dset1_id, error)

       CALL h5fclose_f(file1_id, error)

       CALL h5close_f(error)
       call cpu_time(t3)

       print*,t2-t1,t3-t2



 end subroutine create_H35
  subroutine create_H45(nv,nc,mf_v,mf_c,wcoul0,s,vol,q,p,n4,n5)
      double precision :: fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16

      integer, intent(in) :: nv,nc,s,q,n4,p,n5
      integer :: nb,i,j,alpha,beta,gama,eta,m,n,imatrix,keyword,c1,c2,c,x,A
      double precision            :: element
      double precision ::wcoul0,vol,a1,b1
      complex(kind=8) , dimension(:,:),allocatable                           :: H45
      complex(kind=8), allocatable :: bse_mat_e(:,:,:,:)
      complex(kind=8), allocatable :: bse_mat_h(:,:,:,:)
      complex(kind=8), allocatable :: bse_mat(:,:,:,:)
      CHARACTER(LEN=11), PARAMETER :: filename1 = "singlet.h5" ! File name
      CHARACTER(LEN=5), PARAMETER :: dsetname1 = "CIS45"    ! Dataset name
      CHARACTER(LEN=11), PARAMETER :: filename2 = "singletX.h5" ! File name
      CHARACTER(LEN=5), PARAMETER :: dsetname2 = "HCISD"    ! Dataset name
      INTEGER, PARAMETER   :: DRANK = 3 ! Dataset rank 1 dimension for real or complex 2 is row-index 3 is col-index
      INTEGER(SIZE_T), PARAMETER :: NUMP = 2 ! Number of points selected
      complex(kind=8) ,DIMENSION(nv) ,intent(in) :: mf_v
      complex(kind=8),DiMENSION(nc),intent(in) :: mf_c
      INTEGER(HID_T) :: file1_id       ! File1 identifier
      INTEGER(HID_T) :: dset1_id       ! Dataset1 identifier
      INTEGER(HID_T) :: file2_id       ! File2 identifier
      INTEGER(HID_T) :: dset2_id       ! Dataset2 identifier
      INTEGER(HID_T) :: dataspace1     ! Dataspace identifier
          ! memspace identifier
      INTEGER(HID_T) :: dataspace2     ! Dataspace identifier
            ! memspace identifier
      INTEGER(HID_T) :: fapl_id

      INTEGER(HSIZE_T), DIMENSION(2) :: dimsm
                                                   ! Memory dataspace dimensions
      INTEGER(HSIZE_T), DIMENSION(2) :: dimsf
                                                   ! File dataspace dimensions
      INTEGER(HSIZE_T) :: count(2),offset(2) ! Elements coordinates
                                                      ! in the file
      DOUBLE PRECISION, DIMENSION(2) :: val =(/3.14,1.414/)  ! Values to write
      INTEGER :: memrank = 1  ! Rank of the dataset in memory
      DOUBLE PRECISION,DIMENSION(:,:,:),allocatable :: bufnew

      INTEGER :: error,info  ! Error flag
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims

      integer(HID_T) :: filespace     ! Dataspace identifier in file
      integer(HID_T) :: memspace      ! Dataspace identifier in mem
      integer(HID_T) :: dset_id       ! Dataset identifier
      INTEGER(HID_T) :: plist_id      ! Property list identifier
      !!!!!!!!!!!mpivariables!!!!!!!!!!!!!!!!!!!!!!!!!!!
       integer :: nproc, chunk_size,itr_length
      integer , dimension(:,:) , allocatable :: ij_mat
      integer, dimension(:),allocatable :: row_ij
      integer :: comm, rank, numproc, ierror
      integer :: cij,cab
      !call MPI_INIT(ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)


       keyword = 15
       fac_hd=0
       fac_hx=0
       fac_ed=0
       fac_ex=0
      fac_1=0
      fac_2=0
      fac_3=0
      fac_4=0
      fac_5=0
      fac_6=0
      fac_7=0
      fac_8=0
      fac_9=0
      fac_10=0
      fac_11=0
      fac_12=0
      fac_13=0
      fac_14=0
      fac_15=0
      fac_16=0
       itr_length = nc*(nc-1)/2
      if (mod(nc,2)==0) then
         nproc = nc-1
         chunk_size = nc/2
      else
         nproc = nc
         chunk_size = (nc-1)/2
      end if
      print*,chunk_size,nproc

      allocate(bse_mat_h(nv,nv,nv,nv))
      allocate(bse_mat_e(nc,1,1,nc))
      allocate(bse_mat(nv,nv,2,nc))
      !allocate(row_ij(chunk_size))
      A = (nv*(nv-1))/2
      !allocate(ij_mat(2,chunk_size))
      c = rank+1
      nrows_p= A
      allocate(H45(nrows_p,p))

      H45 = cmplx(0.0,0.0)


      c1=0
      c2=0
      cij=0
      cab=0
      do imatrix = 3,4
         call fac_values(wcoul0,vol,keyword,imatrix,fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16)

         i = c
         cij = 0
         cab = 0
         if (imatrix==3) then
            call load_Hhh(nv,bse_mat_h)
         end if

         if (imatrix .eq. 3) then
                call load_Hee(i,i,nc,nv,bse_mat_e)
                call load_Heh(i,i,nc,nv,imatrix,bse_mat)
         else if (imatrix .eq. 4) then
                call load_Heh(i,i,nc,nv,imatrix,bse_mat)
         end if
         do beta = 1,nv-1
            do alpha = beta+1,nv
               cab = cab + 1
               c1 = cij + cab
               c2=0
               do n = 1,nc-1
                  do m = n+1,nc
                     do gama = 1,nv
                         c2 = c2 + 1
                         if(i==n .and. beta==gama) then
                            if (imatrix .ne. 4) then
                                H45(c1,c2)= H45(c1,c2) + fac_1 * bse_mat(alpha,gama,1,m)
                                H45(c1,c2)=H45(c1,c2) + fac_3 * bse_mat(alpha,gama,1,m)
                                H45(c1,c2)=H45(c1,c2) +fac_9 * bse_mat(alpha,gama,1,m)
                                H45(c1,c2)=H45(c1,c2) + fac_11 * bse_mat(alpha,gama,1,m)
                                    !print*,"3", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                            end if
                            if (imatrix == 4) then
                                H45(c1,c2)= H45(c1,c2) + fac_1 * bse_mat(alpha,gama,1,m)
                                H45(c1,c2)=H45(c1,c2) + fac_3 * bse_mat(alpha,gama,1,m)
                                H45(c1,c2)=H45(c1,c2) +fac_9 * bse_mat(alpha,gama,1,m)
                                H45(c1,c2)=H45(c1,c2) + fac_11 * bse_mat(alpha,gama,1,m)

                            end if
                         end if
                         if(i==m .and. beta==gama) then
                            if (imatrix .ne. 4) then
                                     !print*,"4",H(c1,c2),bse_mat(alpha,eta,1,n)
                               H45(c1,c2)=H45(c1,c2) + fac_2 * bse_mat(alpha,gama,1,n)
                               H45(c1,c2)=H45(c1,c2) + fac_4 * bse_mat(alpha,gama,1,n)
                               H45(c1,c2)=H45(c1,c2) + fac_10 * bse_mat(alpha,gama,1,n)
                               H45(c1,c2)=H45(c1,c2) + fac_12 * bse_mat(alpha,gama,1,n)

                                     !print*,"4", i,j,alpha,beta,m,n,gama,eta
                                     !print*, H(c1,c2),c1,c2
                                     !print*,"fac",fac
                            end if
                            if (imatrix == 4) then
                               H45(c1,c2)=H45(c1,c2) + fac_2 * bse_mat(alpha,gama,1,n)
                               H45(c1,c2)=H45(c1,c2) + fac_4 * bse_mat(alpha,gama,1,n)
                               H45(c1,c2)=H45(c1,c2) + fac_10 * bse_mat(alpha,gama,1,n)
                               H45(c1,c2)=H45(c1,c2) + fac_12 * bse_mat(alpha,gama,1,n)
                            end if
                         end if
                          if(i==n .and. alpha==gama) then
                            if (imatrix .ne. 4) then
                               H45(c1,c2)=H45(c1,c2) + fac_5 * bse_mat(beta,gama,1,m)
                               H45(c1,c2)=H45(c1,c2) + fac_7 * bse_mat(beta,gama,1,m)
                               H45(c1,c2)=H45(c1,c2) + fac_13 * bse_mat(beta,gama,1,m)
                               H45(c1,c2)=H45(c1,c2) + fac_15 * bse_mat(beta,gama,1,m)
                                    !print*, H(c1,c2),c1,c2
                            end if
                            if (imatrix == 4) then
                               H45(c1,c2)=H45(c1,c2) + fac_5 * bse_mat(beta,gama,1,m)
                               H45(c1,c2)=H45(c1,c2) + fac_7 * bse_mat(beta,gama,1,m)
                               H45(c1,c2)=H45(c1,c2) + fac_13 * bse_mat(beta,gama,1,m)
                               H45(c1,c2)=H45(c1,c2) + fac_15 * bse_mat(beta,gama,1,m)
                            end if
                         end if
                         if(i==m .and. alpha==gama) then
                            if (imatrix .ne. 4) then
                                H45(c1,c2)=H45(c1,c2) + fac_6 * bse_mat(beta,gama,1,n)
                                H45(c1,c2)=H45(c1,c2) + fac_8 * bse_mat(beta,gama,1,n)
                                H45(c1,c2)=H45(c1,c2) + fac_14 * bse_mat(beta,gama,1,n)
                                H45(c1,c2)=H45(c1,c2) + fac_16 * bse_mat(beta,gama,1,n)
                                    !print*,"8", i,j,alpha,beta,m,n,gama,eta
                                    !print*, H(c1,c2),c1,c2
                             end if
                             if (imatrix == 4) then
                                H45(c1,c2)=H45(c1,c2) + fac_6 * bse_mat(beta,gama,1,n)
                                H45(c1,c2)=H45(c1,c2) + fac_8 * bse_mat(beta,gama,1,n)
                                H45(c1,c2)=H45(c1,c2) + fac_14 * bse_mat(beta,gama,1,n)
                                H45(c1,c2)=H45(c1,c2) + fac_16 * bse_mat(beta,gama,1,n)
                             end if
                         end if


                      end do
                   end do
                end do
             end do
          end do
       end do
       call cpu_time(t2)
      CALL h5open_f(error)
      !CALL h5fcreate_f(filename1,  H5F_ACC_TRUNC_F, file1_id, error)


      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)

      !call h5fcreate_f(filename1, H5F_ACC_TRUNC_F, file1_id, error, &
       !           access_prp = plist_id)
      CALL h5fopen_f (filename1, H5F_ACC_RDWR_F, file1_id, error,access_prp = plist_id)


       !CALL h5fopen_f (filename2, H5F_ACC_RDWR_F, file2_id, error,access_prp = plist_id)
       CALL h5pclose_f(plist_id, error)
       dimsf(1) = q
       dimsf(2) = p
       CALL h5screate_simple_f(2, dimsf, dataspace1, error)

       CALL h5dcreate_f(file1_id, dsetname1,H5T_NATIVE_DOUBLE,dataspace1,dset1_id, error)
       CALL h5sclose_f(dataspace1, error)
     !
       !CALL h5dopen_f(file2_id, dsetname2, dset2_id, error)
     !
     ! Get  dataspace identifier.
     !
      count(1) = nrows_p
       count(2) = p
       dimsm(1)=nrows_p
       dimsm(2)=p
       dimsf(1) = q
       dimsf(2) = p
       !print*,"rank,count(1),count(2),offset(1),offset(2),nc*nv",rank,count(1),count(2),offset(1),offset(2),nc*nv

       !call h5screate_simple_f(2, dimsm, memspace, error)
       !allocate(data_out(count(1),count(2)))
       !data_out = real(H12)
       offset(1) = rank*nrows_p
       offset(2) = 0
       print*,"rank,count(1),count(2),offset(1),offset(2),nc*nv",rank,count(1),count(2),offset(1),offset(2),nc*nv

       call h5screate_simple_f(2, count, memspace, error)
       CALL h5dget_space_f(dset1_id, dataspace1, error)
       CALL h5sselect_hyperslab_f (dataspace1, H5S_SELECT_SET_F, offset,&
               count, error)


       CALL h5pcreate_f(H5P_DATASET_XFER_F, fapl_id, error)
       CALL h5pset_dxpl_mpio_f(fapl_id, H5FD_MPIO_COLLECTIVE_F, error)
       call h5screate_simple_f(2, dimsm, memspace, error)
       CALL h5dwrite_f(dset1_id, H5T_NATIVE_DOUBLE, real(H45), dimsf , error, &
                    mem_space_id=memspace, file_space_id=dataspace1)
       call h5pclose_f(fapl_id, error)
       CALL h5sclose_f(memspace, error)
       CALL h5dclose_f(dset1_id, error)

       CALL h5fclose_f(file1_id, error)

       CALL h5close_f(error)
       call cpu_time(t3)

       print*,t2-t1,t3-t2


 end subroutine create_H45
  subroutine create_H6(nv,nc,mf_v,mf_c,wcoul0,s,vol,n6)
      double precision :: fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16

      integer, intent(in) :: nv,nc,s
      integer :: nb,i,j,alpha,beta,gama,eta,m,n,imatrix,keyword,c1,c2,c,x,A
      double precision            :: element
      double precision ::wcoul0,vol,a1,b1
      complex(kind=8),dimension(:,:),allocatable                            :: H6
      !complex(kind=8), allocatable :: bse_mat_e(:,:,:,:)
      !complex(kind=8), allocatable :: bse_mat_h(:,:,:,:)
      complex(kind=8), allocatable :: bse_mat(:,:,:,:)
      CHARACTER(LEN=11), PARAMETER :: filename1 = "singlet.h5" ! File name
      CHARACTER(LEN=5), PARAMETER :: dsetname1 = "CISD6"    ! Dataset name
      CHARACTER(LEN=11), PARAMETER :: filename2 = "singletX.h5" ! File name
      CHARACTER(LEN=5), PARAMETER :: dsetname2 = "HCISD"    ! Dataset name
      INTEGER, PARAMETER   :: DRANK = 3 ! Dataset rank 1 dimension for real or complex 2 is row-index 3 is col-index
      INTEGER(SIZE_T), PARAMETER :: NUMP = 2 ! Number of points selected
      complex(kind=8) ,DIMENSION(nv) ,intent(in) :: mf_v
      complex(kind=8),DiMENSION(nc),intent(in) :: mf_c
      INTEGER(HID_T) :: file1_id       ! File1 identifier
      INTEGER(HID_T) :: dset1_id       ! Dataset1 identifier
      INTEGER(HID_T) :: file2_id       ! File2 identifier
      INTEGER(HID_T) :: dset2_id       ! Dataset2 identifier
      INTEGER(HID_T) :: dataspace1     ! Dataspace identifier
          ! memspace identifier
      INTEGER(HID_T) :: dataspace2     ! Dataspace identifier
            ! memspace identifier
      INTEGER(HID_T) :: fapl_id

      INTEGER(HSIZE_T), DIMENSION(2) :: dimsm
                                                   ! Memory dataspace dimensions
      INTEGER(HSIZE_T), DIMENSION(2) :: dimsf
                                                   ! File dataspace dimensions
      INTEGER(HSIZE_T) :: offset(2),count(2) ! Elements coordinates
                                                      ! in the file
      DOUBLE PRECISION, DIMENSION(2) :: val =(/3.14,1.414/)  ! Values to write
      INTEGER :: memrank = 1  ! Rank of the dataset in memory
      DOUBLE PRECISION,DIMENSION(:,:,:),allocatable :: bufnew

      INTEGER :: error,info  ! Error flag
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims

      integer(HID_T) :: filespace     ! Dataspace identifier in file
      integer(HID_T) :: memspace      ! Dataspace identifier in mem
      integer(HID_T) :: dset_id       ! Dataset identifier
      INTEGER(HID_T) :: plist_id      ! Property list identifier
        !!!!!!!!!!!mpivariables!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer :: nproc, chunk_size,itr_length
      integer , dimension(:,:) , allocatable :: ij_mat
      integer, dimension(:),allocatable :: row_ij
      integer :: comm, rank, numproc, ierror
      integer :: cij,cab,c1p,cijp
      !call MPI_INIT(ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)

     
       keyword = 16

       fac_hd=0
       fac_hx=0
       fac_ed=0
       fac_ex=0
      fac_1=0
      fac_2=0
      fac_3=0
      fac_4=0
      fac_5=0
      fac_6=0
      fac_7=0
      fac_8=0
      fac_9=0
      fac_10=0
      fac_11=0
      fac_12=0
      fac_13=0
      fac_14=0
      fac_15=0
      fac_16=0
       itr_length = nc*(nc-1)/2
      if (mod(nc,2)==0) then
         nproc = nc-1
         chunk_size = nc/2
      else
         nproc = nc
         chunk_size = (nc-1)/2
      end if
      print*,chunk_size,nproc

      !allocate(bse_mat_h(nv,nv,nv,nv))
      !allocate(bse_mat_e(nc,1,1,nc))
      allocate(bse_mat(nv,nv,2,nc))
      !allocate(row_ij(chunk_size))
      A = nv
      !allocate(ij_mat(2,chunk_size))
      
      c = rank+1
      nrows_p= A
      allocate(H6(nrows_p,nc*nv))
      H6 = cmplx(0.0,0.0)
      CALL h5screate_simple_f(1, dimsm, memspace, error)

      c1=0
      c2=0
      cij=0
      cab=0
      c1p = 0
      cijp =0
      do imatrix = 3,4
         call fac_values(wcoul0,vol,keyword,imatrix,fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16)

         i = c
         cij = 0
         cijp = (c-1)*A
         cab = 0
         if (imatrix==3) then
            !call load_Hhh(nv,bse_mat_h)
         end if

         if (imatrix .eq. 3) then
                !call load_Hee(i,i,nc,nv,bse_mat_e)
                call load_Heh(i,i,nc,nv,imatrix,bse_mat)
         else if (imatrix .eq. 4) then
                call load_Heh(i,i,nc,nv,imatrix,bse_mat)
         end if

         do alpha = 1,nv
            cab = cab + 1
            c1 = cij + cab
            c1p = cijp + cab
            c2 = 0
            do j = 1,nc
               do beta = 1,nv
                   c2 = c2+1
                   if (imatrix==3 .and. c1p==c2) then
                      H6(c1,c2)=H6(c1,c2)+ (mf_c(i) - mf_v(alpha))
                  end if
                  if (imatrix .ne. 4) then
                      H6(c1,c2)= H6(c1,c2) + fac_1 * bse_mat(alpha,beta,1,j)  ! direct term
                  end if
                  if (imatrix == 4) then
                      H6(c1,c2)= H6(c1,c2) + fac_1 * bse_mat(alpha,beta,1,j)  !exchange term
                  end if
                  !print*,H6(c1,c2),c1,c2,"start of bse"
                 
               end do
            end do
         end do
      end do
      call cpu_time(t2)
      CALL h5open_f(error)
      !CALL h5fcreate_f(filename1,  H5F_ACC_TRUNC_F, file1_id, error)


      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)

      !call h5fcreate_f(filename1, H5F_ACC_TRUNC_F, file1_id, error, &
       !           access_prp = plist_id)
      CALL h5fopen_f (filename1, H5F_ACC_RDWR_F, file1_id, error,access_prp = plist_id)


       !CALL h5fopen_f (filename2, H5F_ACC_RDWR_F, file2_id, error,access_prp = plist_id)
       CALL h5pclose_f(plist_id, error)
       dimsf(1) = nc*nv
       dimsf(2) = nc*nv
       CALL h5screate_simple_f(2, dimsf, dataspace1, error)

       CALL h5dcreate_f(file1_id, dsetname1,H5T_NATIVE_DOUBLE,dataspace1,dset1_id, error)
       CALL h5sclose_f(dataspace1, error)
     !
       !CALL h5dopen_f(file2_id, dsetname2, dset2_id, error)
     !
     ! Get  dataspace identifier.
     !
      count(1) = nrows_p
       count(2) = nc*nv
       dimsm(1)=nrows_p
       dimsm(2)=nc*nv
       dimsf(1) = nc*nv
       dimsf(2) = nc*nv
       !print*,"rank,count(1),count(2),offset(1),offset(2),nc*nv",rank,count(1),count(2),offset(1),offset(2),nc*nv

       !call h5screate_simple_f(2, dimsm, memspace, error)
       !allocate(data_out(count(1),count(2)))
       !data_out = real(H12)
       offset(1) = rank*nrows_p
       offset(2) = 0
       print*,"rank,count(1),count(2),offset(1),offset(2),nc*nv",rank,count(1),count(2),offset(1),offset(2),nc*nv

       call h5screate_simple_f(2, count, memspace, error)
       CALL h5dget_space_f(dset1_id, dataspace1, error)
       CALL h5sselect_hyperslab_f (dataspace1, H5S_SELECT_SET_F, offset,&
               count, error)
        CALL h5pcreate_f(H5P_DATASET_XFER_F, fapl_id, error)
       CALL h5pset_dxpl_mpio_f(fapl_id, H5FD_MPIO_COLLECTIVE_F, error)
       call h5screate_simple_f(2, dimsm, memspace, error)
       CALL h5dwrite_f(dset1_id, H5T_NATIVE_DOUBLE, real(H6), dimsf , error, &
                    mem_space_id=memspace, file_space_id=dataspace1)
       call h5pclose_f(fapl_id, error)

       CALL h5sclose_f(memspace, error)
       CALL h5dclose_f(dset1_id, error)

       CALL h5fclose_f(file1_id, error)

       CALL h5close_f(error)
       call cpu_time(t3)

       print*,t2-t1,t3-t2


     
 end subroutine create_H6
  subroutine create_H16(nv,nc,mf_v,mf_c,wcoul0,s,vol,n6,chunk)
      double precision :: fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16

      integer, intent(in) :: nv,nc,s,n6,chunk
      integer :: nb,i,j,alpha,beta,gama,eta,m,n,imatrix,keyword,c1,c2,c,x,A
      double precision            :: element
      double precision ::wcoul0,vol,a1,b1
      complex(kind=8),dimension(:,:),allocatable                            :: H16
      complex(kind=8), allocatable :: bse_mat_e(:,:,:,:)
      complex(kind=8), allocatable :: bse_mat_h(:,:,:,:)
      !complex(kind=8), allocatable :: bse_mat(:,:,:,:)
      CHARACTER(LEN=11), PARAMETER :: filename1 = "singlet.h5" ! File name
      CHARACTER(LEN=5), PARAMETER :: dsetname1 = "CIS16"    ! Dataset name
      CHARACTER(LEN=11), PARAMETER :: filename2 = "singletX.h5" ! File name
      CHARACTER(LEN=5), PARAMETER :: dsetname2 = "HCISD"    ! Dataset name

      INTEGER, PARAMETER   :: DRANK = 3 ! Dataset rank 1 dimension for real or complex 2 is row-index 3 is col-index
      INTEGER(SIZE_T), PARAMETER :: NUMP = 2 ! Number of points selected
      complex(kind=8) ,DIMENSION(nv) ,intent(in) :: mf_v
      complex(kind=8),DiMENSION(nc),intent(in) :: mf_c
      INTEGER(HID_T) :: file1_id       ! File1 identifier
      INTEGER(HID_T) :: dset1_id       ! Dataset1 identifier
      INTEGER(HID_T) :: file2_id       ! File2 identifier
      INTEGER(HID_T) :: dset2_id       ! Dataset2 identifier
      INTEGER(HID_T) :: dataspace1     ! Dataspace identifier
          ! memspace identifier
      INTEGER(HID_T) :: dataspace2     ! Dataspace identifier
            ! memspace identifier
      INTEGER(HID_T) :: fapl_id

      INTEGER(HSIZE_T), DIMENSION(2) :: dimsm
                                                   ! Memory dataspace dimensions
      INTEGER(HSIZE_T), DIMENSION(2) :: dimsf
                                                   ! File dataspace dimensions
      INTEGER(HSIZE_T) :: count(2),offset(2) ! Elements coordinates
                                                      ! in the file
      DOUBLE PRECISION, DIMENSION(2) :: val =(/3.14,1.414/)  ! Values to write
      INTEGER :: memrank = 1  ! Rank of the dataset in memory
      DOUBLE PRECISION,DIMENSION(:,:,:),allocatable :: bufnew

      INTEGER :: error,info  ! Error flag
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims

      integer(HID_T) :: filespace     ! Dataspace identifier in file
      integer(HID_T) :: memspace      ! Dataspace identifier in mem
      integer(HID_T) :: dset_id       ! Dataset identifier
      INTEGER(HID_T) :: plist_id      ! Property list identifier
       !!!!!!!!!!!mpivariables!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer :: nproc, chunk_size,itr_length
      integer , dimension(:,:) , allocatable :: ij_mat
      integer, dimension(:),allocatable :: row_ij
      integer :: comm, rank, numproc, ierror
      integer :: cij,cab
      !call MPI_INIT(ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)

       keyword = 17

       fac_hd=0
       fac_hx=0
       fac_ed=0
       fac_ex=0
      fac_1=0
      fac_2=0
      fac_3=0
      fac_4=0
      fac_5=0
      fac_6=0
      fac_7=0
      fac_8=0
      fac_9=0
      fac_10=0
      fac_11=0
      fac_12=0
      fac_13=0
      fac_14=0
      fac_15=0
      fac_16=0
       itr_length = nc*(nc-1)/2
       chunk_size = chunk
      print*,chunk_size,nproc
      nb = nv + nc

      allocate(bse_mat_h(nv,nv,2,nv))
      allocate(bse_mat_e(nb,1,1,nb))
      !allocate(bse_mat(nv,nv,2,nc))

      allocate(row_ij(chunk_size))
      A = (nv*(nv-1))/2
      allocate(ij_mat(2,chunk_size))
      c = rank+1
      nrows_p= chunk_size*A
      allocate(H16(nrows_p,nc*nv))
      H16 = cmplx(0.0,0.0)
      


      !print*,"proc_id+1",c
      call distribute_ij(rank,nc,A,chunk_size,ij_mat,row_ij)
      c1=0
      c2=0
      cij=0
      cab=0
      do imatrix = 3,4
       call fac_values(wcoul0,vol,keyword,imatrix,fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16)
         !if (imatrix==3) then
            !call load_s_Hhh(i,j,nv,bse_mat_h)
         !end if
         do x = 1,chunk_size
            print*,x,ij_mat(1,x),ij_mat(2,x)

            i = ij_mat(1,x)
            j = ij_mat(2,x)
            cij = (x-1)*A
            cab =0
             if (imatrix==3) then
                call load_s_Hhh(i,j,nv,bse_mat_h)
             end if


            if (imatrix .eq. 3) then
                call load_s_Hee(i,j,nb,nv,bse_mat_e)
            end if     
            do beta = 1,nv-1
               do alpha = beta+1,nv  
                  cab = cab +1
                  c1 = cij + cab
                  c2 = 0
                  do m = 1,nc
                     do gama = 1,nv
                        c2 = c2 + 1
                        
                        if (j==m) then
                            if (imatrix .ne. 4) then
                                     !print*,"bb",c1,c2, H(5,11),H(71,1),imatrix
                                      !print*,"1",H(c1,c2),bse_mat(eta,beta,alpha,gama), bse_mat(gama,beta,alpha,eta)
                               H16(c1,c2) =H16(c1,c2) + (fac_1 * bse_mat_h(alpha,gama,1,beta)) &          !Fac 1234 stands for direct terms fac 5678 stands for exchange terms             
                                       + (fac_5 * bse_mat_h(beta,gama,1,alpha))
                            end if
                         end if
                         if (i==m) then
                            if (imatrix .ne. 4) then
                                     !print*,"bb",c1,c2, H(5,11),H(71,1),imatrix
                                      !print*,"1",H(c1,c2),bse_mat(eta,beta,alpha,gama), bse_mat(gama,beta,alpha,eta)
                               H16(c1,c2) =H16(c1,c2) + (fac_2 * bse_mat_h(alpha,gama,2,beta)) &          !Fac 1234 stands for direct terms fac 5678 stands for exchange terms             
                                       + (fac_6 * bse_mat_h(beta,gama,2,alpha))
                            end if
                         end if
                         if (beta == gama ) then
                            if (imatrix .ne. 4) then
                                     !print*,"bb",c1,c2, H(5,11),H(71,1),imatrix
                                      !print*,"1",H(c1,c2),bse_mat(eta,beta,alpha,gama), bse_mat(gama,beta,alpha,eta)
                               H16(c1,c2) =H16(c1,c2) + (fac_3 * bse_mat_e(alpha,1,1,m+nv)) &          !Fac 1234 stands for direct terms fac 5678 stands for exchange terms             
                                       + (fac_7 * bse_mat_e(m+nv,1,1,alpha))
                            end if
                         end if
                         if (alpha == gama ) then
                            if (imatrix .ne. 4) then
                                     !print*,"bb",c1,c2, H(5,11),H(71,1),imatrix
                                      !print*,"1",H(c1,c2),bse_mat(eta,beta,alpha,gama), bse_mat(gama,beta,alpha,eta)
                               H16(c1,c2) =H16(c1,c2) + (fac_4 * bse_mat_e(beta,1,1,m+nv)) &          !Fac 1234 stands for direct terms fac 5678 stands for exchange terms             
                                       + (fac_8 * bse_mat_e(m+nv,1,1,beta))
                            end if
                         end if                        
                      end do
                   end do
                end do
             end do
         end do
      end do
      CALL h5open_f(error)
      !CALL h5fopen_f(filename1,  H5F_ACC_TRUNC_F, file1_id, error)


      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)

      !call h5fcreate_f(filename1, H5F_ACC_TRUNC_F, file1_id, error, &
       !           access_prp = plist_id)
      CALL h5fopen_f (filename1, H5F_ACC_RDWR_F, file1_id, error,access_prp = plist_id)


       !CALL h5fopen_f (filename2, H5F_ACC_RDWR_F, file2_id, error,access_prp = plist_id)
       CALL h5pclose_f(plist_id, error)
       dimsf(1) = s
       dimsf(2) = nc*nv
       CALL h5screate_simple_f(2, dimsf, dataspace1, error)

       CALL h5dcreate_f(file1_id, dsetname1,H5T_NATIVE_DOUBLE,dataspace1,dset1_id, error)
       CALL h5sclose_f(dataspace1, error)
     !
       !CALL h5dopen_f(file2_id, dsetname2, dset2_id, error)
     !
     ! Get  dataspace identifier.
     !
       !CALL h5dget_space_f(dset1_id, dataspace1, error)
       !CALL h5dget_space_f(dset2_id, dataspace2, error)
      !call h5screate_simple_f(2, count, memspace, error)
       count(1) = nrows_p
       count(2) = nc*nv
       dimsm(1)=nrows_p
       dimsm(2)=nc*nv
       dimsf(1) = s
       dimsf(2) = nc*nv
       !call h5screate_simple_f(2, dimsm, memspace, error)
       !allocate(data_out(count(1),count(2)))
       !data_out = real(H1)
       offset(1) = rank*nrows_p
       offset(2) = 0
       print*,"rank,count(1),count(2),offset(1),offset(2),s",rank,count(1),count(2),offset(1),offset(2),s
       !CALL h5dget_space_f(dset1_id, dataspace1, error)

       !CALL h5dget_space_f(dset1_id, dataspace1, error)
       call h5screate_simple_f(2, count, memspace, error)
      ! CALL h5sselect_hyperslab_f (memspace, H5S_SELECT_SET_F, offset,&
       !        count, error)
       CALL h5dget_space_f(dset1_id, dataspace1, error)
        CALL h5sselect_hyperslab_f (dataspace1, H5S_SELECT_SET_F, offset,&
               count, error)


      ! CALL h5pcreate_f(H5P_DATASET_XFER_F, fapl_id, error)
       !CALL h5pset_dxpl_mpio_f(fapl_id, H5FD_MPIO_COLLECTIVE_F, error)
       !call h5screate_simple_f(2, dimsm, memspace, error)
       CALL h5dwrite_f(dset1_id, H5T_NATIVE_DOUBLE, real(H16), dimsf , error, &
                    mem_space_id=memspace, file_space_id=dataspace1)
       !call h5pclose_f(fapl_id, error)

       CALL h5sclose_f(memspace, error)
       !call h5pclose_f(fapl_id, error)

       CALL h5dclose_f(dset1_id, error)

       CALL h5fclose_f(file1_id, error)

       CALL h5close_f(error)
                                                         
 end subroutine create_H16
  subroutine create_H26(nv,nc,mf_v,mf_c,wcoul0,s,vol,n6,chunk)
      double precision :: fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16

      integer, intent(in) :: nv,nc,s,n6,chunk
      integer :: nb,i,j,alpha,beta,gama,eta,m,n,imatrix,keyword,c1,c2,c,x,A
      double precision            :: element
      double precision ::wcoul0,vol,a1,b1
      complex(kind=8),dimension(:,:),allocatable                            :: H26
      complex(kind=8), allocatable :: bse_mat_e(:,:,:,:)
      complex(kind=8), allocatable :: bse_mat_h(:,:,:,:)
      !complex(kind=8), allocatable :: bse_mat(:,:,:,:)
      CHARACTER(LEN=11), PARAMETER :: filename1 = "singlet.h5" ! File name
      CHARACTER(LEN=5), PARAMETER :: dsetname1 = "CIS26"    ! Dataset name
      CHARACTER(LEN=11), PARAMETER :: filename2 = "singletX.h5" ! File name
      CHARACTER(LEN=5), PARAMETER :: dsetname2 = "HCISD"    ! Dataset name

      INTEGER, PARAMETER   :: DRANK = 3 ! Dataset rank 1 dimension for real or complex 2 is row-index 3 is col-index
      INTEGER(SIZE_T), PARAMETER :: NUMP = 2 ! Number of points selected
      complex(kind=8) ,DIMENSION(nv) ,intent(in) :: mf_v
      complex(kind=8),DiMENSION(nc),intent(in) :: mf_c
      INTEGER(HID_T) :: file1_id       ! File1 identifier
      INTEGER(HID_T) :: dset1_id       ! Dataset1 identifier
      INTEGER(HID_T) :: file2_id       ! File2 identifier
      INTEGER(HID_T) :: dset2_id       ! Dataset2 identifier
      INTEGER(HID_T) :: dataspace1     ! Dataspace identifier
          ! memspace identifier
      INTEGER(HID_T) :: dataspace2     ! Dataspace identifier
            ! memspace identifier
      INTEGER(HID_T) :: fapl_id

      INTEGER(HSIZE_T), DIMENSION(2) :: dimsm
                                                   ! Memory dataspace dimensions
      INTEGER(HSIZE_T), DIMENSION(2) :: dimsf
                                                   ! File dataspace dimensions
      INTEGER(HSIZE_T) :: count(2),offset(2) ! Elements coordinates
                                                      ! in the file
      DOUBLE PRECISION, DIMENSION(2) :: val =(/3.14,1.414/)  ! Values to write
      INTEGER :: memrank = 1  ! Rank of the dataset in memory
      DOUBLE PRECISION,DIMENSION(:,:,:),allocatable :: bufnew

      INTEGER :: error,info  ! Error flag
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims

      integer(HID_T) :: filespace     ! Dataspace identifier in file
      integer(HID_T) :: memspace      ! Dataspace identifier in mem
      integer(HID_T) :: dset_id       ! Dataset identifier
      INTEGER(HID_T) :: plist_id      ! Property list identifier
       !!!!!!!!!!!mpivariables!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer :: nproc, chunk_size,itr_length
      integer , dimension(:,:) , allocatable :: ij_mat
      integer, dimension(:),allocatable :: row_ij
      integer :: comm, rank, numproc, ierror
      integer :: cij,cab
      !call MPI_INIT(ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)

       keyword = 18

       fac_hd=0
       fac_hx=0
       fac_ed=0
       fac_ex=0
      fac_1=0
      fac_2=0
      fac_3=0
      fac_4=0
      fac_5=0
      fac_6=0
      fac_7=0
      fac_8=0
      fac_9=0
      fac_10=0
      fac_11=0
      fac_12=0
      fac_13=0
      fac_14=0
      fac_15=0
      fac_16=0
       itr_length = nc*(nc-1)/2
       chunk_size = chunk 
      print*,chunk_size,nproc
      nb = nv + nc

      allocate(bse_mat_h(nv,nv,2,nv))
      allocate(bse_mat_e(nb,1,1,nb))
      !allocate(bse_mat(nv,nv,2,nc))

      allocate(row_ij(chunk_size))
      A = (nv*(nv-1))/2
      allocate(ij_mat(2,chunk_size))
      c = rank+1
      nrows_p= chunk_size*A
      allocate(H26(nrows_p,nc*nv))
      H26 = cmplx(0.0,0.0)
      


      !print*,"proc_id+1",c
      call distribute_ij(rank,nc,A,chunk_size,ij_mat,row_ij)
      c1=0
      c2=0
      cij=0
      cab=0
      do imatrix = 3,4
       call fac_values(wcoul0,vol,keyword,imatrix,fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16)
         !if (imatrix==3) then
            !call load_s_Hhh(i,j,nv,bse_mat_h)
         !end if
         do x = 1,chunk_size
            print*,x,ij_mat(1,x),ij_mat(2,x)

            i = ij_mat(1,x)
            j = ij_mat(2,x)
            cij = (x-1)*A
            cab =0
             if (imatrix==3) then
                call load_s_Hhh(i,j,nv,bse_mat_h)
             end if


            if (imatrix .eq. 3) then
                call load_s_Hee(i,j,nb,nv,bse_mat_e)
            end if     
            do beta = 1,nv-1
               do alpha = beta+1,nv  
                  cab = cab +1
                  c1 = cij + cab
                  c2 = 0
                  do m = 1,nc
                     do gama = 1,nv
                        c2 = c2 + 1
                        
                        if (j==m) then
                            if (imatrix .ne. 4) then
                                     !print*,"bb",c1,c2, H(5,11),H(71,1),imatrix
                                      !print*,"1",H(c1,c2),bse_mat(eta,beta,alpha,gama), bse_mat(gama,beta,alpha,eta)
                               H26(c1,c2) =H26(c1,c2) + (fac_1 * bse_mat_h(alpha,gama,1,beta)) &          !Fac 1234 stands for direct terms fac 5678 stands for exchange terms             
                                       + (fac_5 * bse_mat_h(beta,gama,1,alpha))
                            end if
                         end if
                         if (i==m) then
                            if (imatrix .ne. 4) then
                                     !print*,"bb",c1,c2, H(5,11),H(71,1),imatrix
                                      !print*,"1",H(c1,c2),bse_mat(eta,beta,alpha,gama), bse_mat(gama,beta,alpha,eta)
                               H26(c1,c2) =H26(c1,c2) + (fac_2 * bse_mat_h(alpha,gama,2,beta)) &          !Fac 1234 stands for direct terms fac 5678 stands for exchange terms             
                                       + (fac_6 * bse_mat_h(beta,gama,2,alpha))
                            end if
                         end if
                         if (beta == gama ) then
                            if (imatrix .ne. 4) then
                                     !print*,"bb",c1,c2, H(5,11),H(71,1),imatrix
                                      !print*,"1",H(c1,c2),bse_mat(eta,beta,alpha,gama), bse_mat(gama,beta,alpha,eta)
                               H26(c1,c2) =H26(c1,c2) + (fac_3 * bse_mat_e(alpha,1,1,m+nv)) &          !Fac 1234 stands for direct terms fac 5678 stands for exchange terms             
                                       + (fac_7 * bse_mat_e(m+nv,1,1,alpha))
                            end if
                         end if
                         if (alpha == gama ) then
                            if (imatrix .ne. 4) then
                                     !print*,"bb",c1,c2, H(5,11),H(71,1),imatrix
                                      !print*,"1",H(c1,c2),bse_mat(eta,beta,alpha,gama), bse_mat(gama,beta,alpha,eta)
                               H26(c1,c2) =H26(c1,c2) + (fac_4 * bse_mat_e(beta,1,1,m+nv)) &          !Fac 1234 stands for direct terms fac 5678 stands for exchange terms             
                                       + (fac_8 * bse_mat_e(m+nv,1,1,beta))
                            end if
                         end if                        
                      end do
                   end do
                end do
             end do
         end do
      end do
      CALL h5open_f(error)
      !CALL h5fcreate_f(filename1,  H5F_ACC_TRUNC_F, file1_id, error)


      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)

      !call h5fcreate_f(filename1, H5F_ACC_TRUNC_F, file1_id, error, &
      !            access_prp = plist_id)
      CALL h5fopen_f (filename1, H5F_ACC_RDWR_F, file1_id, error,access_prp = plist_id)


       !CALL h5fopen_f (filename2, H5F_ACC_RDWR_F, file2_id, error,access_prp = plist_id)
       CALL h5pclose_f(plist_id, error)
       dimsf(1) = s
       dimsf(2) = nc*nv
       CALL h5screate_simple_f(2, dimsf, dataspace1, error)

       CALL h5dcreate_f(file1_id, dsetname1,H5T_NATIVE_DOUBLE,dataspace1,dset1_id, error)
       CALL h5sclose_f(dataspace1, error)
     !
       !CALL h5dopen_f(file2_id, dsetname2, dset2_id, error)
     !
     ! Get  dataspace identifier.
     !
       !CALL h5dget_space_f(dset1_id, dataspace1, error)
       !CALL h5dget_space_f(dset2_id, dataspace2, error)
      !call h5screate_simple_f(2, count, memspace, error)
       count(1) = nrows_p
       count(2) = nc*nv
       dimsm(1)=nrows_p
       dimsm(2)=nc*nv
       dimsf(1) = s
       dimsf(2) = nc*nv
       !call h5screate_simple_f(2, dimsm, memspace, error)
       !allocate(data_out(count(1),count(2)))
       !data_out = real(H1)
       offset(1) = rank*nrows_p
       offset(2) = 0
       print*,"rank,count(1),count(2),offset(1),offset(2),s",rank,count(1),count(2),offset(1),offset(2),s
       !CALL h5dget_space_f(dset1_id, dataspace1, error)

       !CALL h5dget_space_f(dset1_id, dataspace1, error)
       call h5screate_simple_f(2, count, memspace, error)
      ! CALL h5sselect_hyperslab_f (memspace, H5S_SELECT_SET_F, offset,&
       !        count, error)
       CALL h5dget_space_f(dset1_id, dataspace1, error)
        CALL h5sselect_hyperslab_f (dataspace1, H5S_SELECT_SET_F, offset,&
               count, error)


      ! CALL h5pcreate_f(H5P_DATASET_XFER_F, fapl_id, error)
       !CALL h5pset_dxpl_mpio_f(fapl_id, H5FD_MPIO_COLLECTIVE_F, error)
       !call h5screate_simple_f(2, dimsm, memspace, error)
       CALL h5dwrite_f(dset1_id, H5T_NATIVE_DOUBLE, real(H26), dimsf , error, &
                    mem_space_id=memspace, file_space_id=dataspace1)
       !call h5pclose_f(fapl_id, error)

       CALL h5sclose_f(memspace, error)
       !call h5pclose_f(fapl_id, error)

       CALL h5dclose_f(dset1_id, error)

       CALL h5fclose_f(file1_id, error)

       CALL h5close_f(error)
                                                         
 end subroutine create_H26
  subroutine create_H36(nv,nc,mf_v,mf_c,wcoul0,s,vol,n6)
      double precision :: fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16

      integer, intent(in) :: nv,nc,s,n6
      integer :: nb,i,j,alpha,beta,gama,eta,m,n,imatrix,keyword,c1,c2,c,x,A
      double precision            :: element
      double precision ::wcoul0,vol,a1,b1
      complex(kind=8),dimension(:,:),allocatable   :: H36
      complex(kind=8), allocatable :: bse_mat_e(:,:,:,:)
      complex(kind=8), allocatable :: bse_mat_h(:,:,:,:)
      !complex(kind=8), allocatable :: bse_mat(:,:,:,:)
      CHARACTER(LEN=11), PARAMETER :: filename1 = "singlet.h5" ! File name
      CHARACTER(LEN=5), PARAMETER :: dsetname1 = "CIS36"    ! Dataset name
      CHARACTER(LEN=11), PARAMETER :: filename2 = "singletX.h5" ! File name
      CHARACTER(LEN=5), PARAMETER :: dsetname2 = "HCISD"    ! Dataset name
      INTEGER, PARAMETER   :: DRANK = 3 ! Dataset rank 1 dimension for real or complex 2 is row-index 3 is col-index
      INTEGER(SIZE_T), PARAMETER :: NUMP = 2 ! Number of points selected
      complex(kind=8) ,DIMENSION(nv) ,intent(in) :: mf_v
      complex(kind=8),DiMENSION(nc),intent(in) :: mf_c
      INTEGER(HID_T) :: file1_id       ! File1 identifier
      INTEGER(HID_T) :: dset1_id       ! Dataset1 identifier
      INTEGER(HID_T) :: file2_id       ! File2 identifier
      INTEGER(HID_T) :: dset2_id       ! Dataset2 identifier
      INTEGER(HID_T) :: dataspace1     ! Dataspace identifier
          ! memspace identifier
      INTEGER(HID_T) :: dataspace2     ! Dataspace identifier
            ! memspace identifier
      INTEGER(HID_T) :: fapl_id

      INTEGER(HSIZE_T), DIMENSION(2) :: dimsm
                                                   ! Memory dataspace dimensions
      INTEGER(HSIZE_T), DIMENSION(2) :: dimsf
                                                   ! File dataspace dimensions
      INTEGER(HSIZE_T) :: offset(2),count(2) ! Elements coordinates
                                                      ! in the file
      DOUBLE PRECISION, DIMENSION(2) :: val =(/3.14,1.414/)  ! Values to write
      INTEGER :: memrank = 1  ! Rank of the dataset in memory
      DOUBLE PRECISION,DIMENSION(:,:,:),allocatable :: bufnew

      INTEGER :: error,info  ! Error flag
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims

      integer(HID_T) :: filespace     ! Dataspace identifier in file
      integer(HID_T) :: memspace      ! Dataspace identifier in mem
      integer(HID_T) :: dset_id       ! Dataset identifier
      INTEGER(HID_T) :: plist_id      ! Property list identifier
       !!!!!!!!!!!mpivariables!!!!!!!!!!!!!!!!!!!!!!!!!!!
       integer :: nproc, chunk_size,itr_length
      integer , dimension(:,:) , allocatable :: ij_mat
      integer, dimension(:),allocatable :: row_ij
      integer :: comm, rank, numproc, ierror
      integer :: cij,cab
      !call MPI_INIT(ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)     

       keyword = 19
       fac_hd=0
       fac_hx=0
       fac_ed=0
       fac_ex=0
      fac_1=0
      fac_2=0
      fac_3=0
      fac_4=0
      fac_5=0
      fac_6=0
      fac_7=0
      fac_8=0
      fac_9=0
      fac_10=0
      fac_11=0
      fac_12=0
      fac_13=0
      fac_14=0
      fac_15=0
      fac_16=0
      nb = nc+nv
         itr_length = nc*(nc-1)/2
      if (mod(nc,2)==0) then
         nproc = nc-1
         chunk_size = nc/2
      else
         nproc = nc
         chunk_size = (nc-1)/2
      end if
      print*,chunk_size,nproc
      allocate(bse_mat_h(nv,nv,2,nv))
      allocate(bse_mat_e(nb,1,1,nb))
      A = nv
      !allocate(ij_mat(2,chunk_size))
      c = rank+1
      nrows_p= A
      allocate(H36(nrows_p,nc*nv))

      H36 = cmplx(0.0,0.0)
  

      c1=0
      c2=0
      cij=0
      cab=0
      do imatrix = 3,4
         call fac_values(wcoul0,vol,keyword,imatrix,fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16)

         i = c
         cij = 0
         cab = 0
         if (imatrix==3) then
                call load_s_Hhh(i,i,nv,bse_mat_h)
          end if


         if (imatrix .eq. 3) then
               call load_s_Hee(i,i,nb,nv,bse_mat_e)
          end if


          do alpha = 1,nv
            cab = cab + 1
            c1 = cij + cab
            c2 = 0
            do m = 1,nc
               do gama = 1,nv
                  c2 = c2 + 1
                  if (i==m) then
                       if (imatrix .ne. 4) then
                                     !print*,"bb",c1,c2, H(5,11),H(71,1),imatrix
                                      !print*,"1",H(c1,c2),bse_mat(eta,beta,alpha,gama), bse_mat(gama,beta,alpha,eta)
                           H36(c1,c2) =H36(c1,c2) + (fac_1 * bse_mat_h(alpha,gama,1,alpha)) &          !Fac 1234 stands for direct terms fac 5678 stands for exchange terms
                                       + (fac_5 * bse_mat_h(alpha,gama,1,alpha))
                           H36(c1,c2) =H36(c1,c2) + (fac_2 * bse_mat_h(alpha,gama,1,alpha)) &          !Fac 1234 stands for direct terms fac 5678 stands for exchange terms
                                       + (fac_6 * bse_mat_h(alpha,gama,1,alpha))
                       end if
                   end if
                   if (alpha == gama ) then
                      if (imatrix .ne. 4) then
                                     !print*,"bb",c1,c2, H(5,11),H(71,1),imatrix
                                      !print*,"1",H(c1,c2),bse_mat(eta,beta,alpha,gama), bse_mat(gama,beta,alpha,eta)
                           H36(c1,c2) =H36(c1,c2) + (fac_3 * bse_mat_e(alpha,1,1,m+nv)) &          !Fac 1234 stands for direct terms fac 5678 stands for exchange terms
                                       + (fac_7 * bse_mat_e(m+nv,1,1,alpha))
                           H36(c1,c2) =H36(c1,c2) + (fac_4 * bse_mat_e(alpha,1,1,m+nv)) &          !Fac 1234 stands for direct terms fac 5678 stands for exchange terms
                                       + (fac_8 * bse_mat_e(m+nv,1,1,alpha))
                      end if
                   end if                  
                 end do
              end do
           end do
       end do
        call cpu_time(t2)
      CALL h5open_f(error)
      !CALL h5fcreate_f(filename1,  H5F_ACC_TRUNC_F, file1_id, error)


      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)

      !call h5fcreate_f(filename1, H5F_ACC_TRUNC_F, file1_id, error, &
       !           access_prp = plist_id)
      CALL h5fopen_f (filename1, H5F_ACC_RDWR_F, file1_id, error,access_prp = plist_id)


       !CALL h5fopen_f (filename2, H5F_ACC_RDWR_F, file2_id, error,access_prp = plist_id)
       CALL h5pclose_f(plist_id, error)
       dimsf(1) = nc*nv
       dimsf(2) = nc*nv
       CALL h5screate_simple_f(2, dimsf, dataspace1, error)

       CALL h5dcreate_f(file1_id, dsetname1,H5T_NATIVE_DOUBLE,dataspace1,dset1_id, error)
       CALL h5sclose_f(dataspace1, error)
     !
       !CALL h5dopen_f(file2_id, dsetname2, dset2_id, error)
     !
     ! Get  dataspace identifier.
     !
      count(1) = nrows_p
       count(2) = nc*nv
       dimsm(1)=nrows_p
       dimsm(2)=nc*nv
       dimsf(1) = nc*nv
       dimsf(2) = nc*nv
       !print*,"rank,count(1),count(2),offset(1),offset(2),nc*nv",rank,count(1),count(2),offset(1),offset(2),nc*nv

       !call h5screate_simple_f(2, dimsm, memspace, error)
       !allocate(data_out(count(1),count(2)))
       !data_out = real(H12)
       offset(1) = rank*nrows_p
       offset(2) = 0
       print*,"rank,count(1),count(2),offset(1),offset(2),nc*nv",rank,count(1),count(2),offset(1),offset(2),nc*nv

       call h5screate_simple_f(2, count, memspace, error)
       CALL h5dget_space_f(dset1_id, dataspace1, error)
       CALL h5sselect_hyperslab_f (dataspace1, H5S_SELECT_SET_F, offset,&
               count, error)
       CALL h5pcreate_f(H5P_DATASET_XFER_F, fapl_id, error)
       CALL h5pset_dxpl_mpio_f(fapl_id, H5FD_MPIO_COLLECTIVE_F, error)
       call h5screate_simple_f(2, dimsm, memspace, error)
       CALL h5dwrite_f(dset1_id, H5T_NATIVE_DOUBLE, real(H36), dimsf , error, &
                    mem_space_id=memspace, file_space_id=dataspace1)
       call h5pclose_f(fapl_id, error)

       CALL h5sclose_f(memspace, error)
       CALL h5dclose_f(dset1_id, error)

       CALL h5fclose_f(file1_id, error)

       CALL h5close_f(error)
       call cpu_time(t3)

       print*,t2-t1,t3-t2

 end subroutine create_H36
 subroutine create_H46(nv,nc,mf_v,mf_c,wcoul0,s,vol,q,n4,n6)
      double precision :: fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16

      integer, intent(in) :: nv,nc,s,q,n4,n6
      integer :: nb,i,j,alpha,beta,gama,eta,m,n,imatrix,keyword,c1,c2,c,x,A
      double precision            :: element
      double precision ::wcoul0,vol,a1,b1
      complex(kind=8),dimension(:,:),allocatable                            :: H46
      complex(kind=8), allocatable :: bse_mat_e(:,:,:,:)
      complex(kind=8), allocatable :: bse_mat_h(:,:,:,:)
      complex(kind=8), allocatable :: bse_mat(:,:,:,:)
      CHARACTER(LEN=11), PARAMETER :: filename1 = "singlet.h5" ! File name
      CHARACTER(LEN=5), PARAMETER :: dsetname1 = "CIS46"    ! Dataset name
      CHARACTER(LEN=11), PARAMETER :: filename2 = "singletX.h5" ! File name
      CHARACTER(LEN=5), PARAMETER :: dsetname2 = "HCISD"    ! Dataset name
      INTEGER, PARAMETER   :: DRANK = 3 ! Dataset rank 1 dimension for real or complex 2 is row-index 3 is col-index
      INTEGER(SIZE_T), PARAMETER :: NUMP = 2 ! Number of points selected
      complex(kind=8) ,DIMENSION(nv) ,intent(in) :: mf_v
      complex(kind=8),DiMENSION(nc),intent(in) :: mf_c
      INTEGER(HID_T) :: file1_id       ! File1 identifier
      INTEGER(HID_T) :: dset1_id       ! Dataset1 identifier
      INTEGER(HID_T) :: file2_id       ! File2 identifier
      INTEGER(HID_T) :: dset2_id       ! Dataset2 identifier
      INTEGER(HID_T) :: dataspace1     ! Dataspace identifier
          ! memspace identifier
      INTEGER(HID_T) :: dataspace2     ! Dataspace identifier
            ! memspace identifier
      INTEGER(HID_T) :: fapl_id

      INTEGER(HSIZE_T), DIMENSION(2) :: dimsm
                                                   ! Memory dataspace dimensions
      INTEGER(HSIZE_T), DIMENSION(2) :: dimsf
                                                   ! File dataspace dimensions
      INTEGER(HSIZE_T) :: count(2),offset(2) ! Elements coordinates
                                                      ! in the file
      DOUBLE PRECISION, DIMENSION(2) :: val =(/3.14,1.414/)  ! Values to write
      INTEGER :: memrank = 1  ! Rank of the dataset in memory
      DOUBLE PRECISION,DIMENSION(:,:,:),allocatable :: bufnew

      INTEGER :: error,info  ! Error flag
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims

      integer(HID_T) :: filespace     ! Dataspace identifier in file
      integer(HID_T) :: memspace      ! Dataspace identifier in mem
      integer(HID_T) :: dset_id       ! Dataset identifier
      INTEGER(HID_T) :: plist_id      ! Property list identifier
      !!!!!!!!!!!mpivariables!!!!!!!!!!!!!!!!!!!!!!!!!!!
       integer :: nproc, chunk_size,itr_length
      integer , dimension(:,:) , allocatable :: ij_mat
      integer, dimension(:),allocatable :: row_ij
      integer :: comm, rank, numproc, ierror
      integer :: cij,cab
      !call MPI_INIT(ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)

     

       keyword = 20
       fac_hd=0
       fac_hx=0
       fac_ed=0
       fac_ex=0
      fac_1=0
      fac_2=0
      fac_3=0
      fac_4=0
      fac_5=0
      fac_6=0
      fac_7=0
      fac_8=0
      fac_9=0
      fac_10=0
      fac_11=0
      fac_12=0
      fac_13=0
      fac_14=0
      fac_15=0
      fac_16=0
      nb = nc+nv
       itr_length = nc*(nc-1)/2
      if (mod(nc,2)==0) then
         nproc = nc-1
         chunk_size = nc/2
      else
         nproc = nc
         chunk_size = (nc-1)/2
      end if
      print*,chunk_size,nproc
      allocate(bse_mat_h(nv,nv,2,nv))
      allocate(bse_mat_e(nb,1,1,nb))
      A = (nv*(nv-1))/2
      !allocate(ij_mat(2,chunk_size))
      c = rank+1
      nrows_p= A
      allocate(H46(nrows_p,nc*nv))
      H46 = cmplx(0.0,0.0)
   

      c1=0
      c2=0
      cij=0
      cab=0
      do imatrix = 3,4
         call fac_values(wcoul0,vol,keyword,imatrix,fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16)

         i = c
         cij = 0
         cab = 0
         if (imatrix==3) then
                call load_s_Hhh(i,i,nv,bse_mat_h)
          end if


         if (imatrix .eq. 3) then
               call load_s_Hee(i,i,nb,nv,bse_mat_e)
          end if
         do beta = 1,nv-1
            do alpha = beta+1,nv
               cab = cab + 1
               c1 = cij + cab
               c2=0
               do m = 1,nc
                  do gama = 1,nv
                     c2 = c2 + 1
                     if (i==m) then
                         if (imatrix .ne. 4) then
                                     !print*,"bb",c1,c2, H(5,11),H(71,1),imatrix
                                      !print*,"1",H(c1,c2),bse_mat(eta,beta,alpha,gama), bse_mat(gama,beta,alpha,eta)
                               H46(c1,c2) =H46(c1,c2) + (fac_1 * bse_mat_h(alpha,gama,1,beta)) &          !Fac 1234 stands for direct terms fac 5678 stands for exchange terms
                                       + (fac_5 * bse_mat_h(beta,gama,1,alpha))
                               H46(c1,c2) =H46(c1,c2) + (fac_2 * bse_mat_h(alpha,gama,1,beta)) &          !Fac 1234 stands for direct terms fac 5678 stands for exchange terms
                                       + (fac_6 * bse_mat_h(beta,gama,1,alpha))
                         end if
                     end if
                     if (beta == gama) then
                         if (imatrix .ne. 4) then
                                     !print*,"bb",c1,c2, H(5,11),H(71,1),imatrix
                                      !print*,"1",H(c1,c2),bse_mat(eta,beta,alpha,gama), bse_mat(gama,beta,alpha,eta)
                               H46(c1,c2) =H46(c1,c2) + (fac_3 * bse_mat_e(alpha,1,1,m+nv)) &          !Fac 1234 stands for direct terms fac 5678 stands for exchange terms
                                       + (fac_7 * bse_mat_e(m+nv,1,1,alpha))
                         end if
                     end if
                     if (alpha == gama ) then
                         if (imatrix .ne. 4) then
                                     !print*,"bb",c1,c2, H(5,11),H(71,1),imatrix
                                      !print*,"1",H(c1,c2),bse_mat(eta,beta,alpha,gama), bse_mat(gama,beta,alpha,eta)
                               H46(c1,c2) =H46(c1,c2) + (fac_4 * bse_mat_e(beta,1,1,m+nv)) &          !Fac 1234 stands for direct terms fac 5678 stands for exchange terms
                                       + (fac_8 * bse_mat_e(m+nv,1,1,beta))
                         end if
                     end if
                     
                 end do
               end do
            end do
         end do
      end do
      call cpu_time(t2)
      CALL h5open_f(error)
      !CALL h5fcreate_f(filename1,  H5F_ACC_TRUNC_F, file1_id, error)


      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)

      !call h5fcreate_f(filename1, H5F_ACC_TRUNC_F, file1_id, error, &
       !           access_prp = plist_id)
      CALL h5fopen_f (filename1, H5F_ACC_RDWR_F, file1_id, error,access_prp = plist_id)


       !CALL h5fopen_f (filename2, H5F_ACC_RDWR_F, file2_id, error,access_prp = plist_id)
       CALL h5pclose_f(plist_id, error)
       dimsf(1) = q
       dimsf(2) = nc*nv
       CALL h5screate_simple_f(2, dimsf, dataspace1, error)

       CALL h5dcreate_f(file1_id, dsetname1,H5T_NATIVE_DOUBLE,dataspace1,dset1_id, error)
       CALL h5sclose_f(dataspace1, error)
     !
       !CALL h5dopen_f(file2_id, dsetname2, dset2_id, error)
     !
     ! Get  dataspace identifier.
     !
      count(1) = nrows_p
       count(2) = nc*nv
       dimsm(1)=nrows_p
       dimsm(2)=nc*nv
       dimsf(1) = q
       dimsf(2) = nc*nv
       !print*,"rank,count(1),count(2),offset(1),offset(2),nc*nv",rank,count(1),count(2),offset(1),offset(2),nc*nv

       !call h5screate_simple_f(2, dimsm, memspace, error)
       !allocate(data_out(count(1),count(2)))
       !data_out = real(H12)
       offset(1) = rank*nrows_p
       offset(2) = 0
       print*,"rank,count(1),count(2),offset(1),offset(2),nc*nv",rank,count(1),count(2),offset(1),offset(2),nc*nv

       call h5screate_simple_f(2, count, memspace, error)
       CALL h5dget_space_f(dset1_id, dataspace1, error)
       CALL h5sselect_hyperslab_f (dataspace1, H5S_SELECT_SET_F, offset,&
               count, error)
       CALL h5pcreate_f(H5P_DATASET_XFER_F, fapl_id, error)
       CALL h5pset_dxpl_mpio_f(fapl_id, H5FD_MPIO_COLLECTIVE_F, error)
       call h5screate_simple_f(2, dimsm, memspace, error)
       CALL h5dwrite_f(dset1_id, H5T_NATIVE_DOUBLE, real(H46), dimsf , error, &
                    mem_space_id=memspace, file_space_id=dataspace1)
       call h5pclose_f(fapl_id, error)
       CALL h5sclose_f(memspace, error)
       CALL h5dclose_f(dset1_id, error)

       CALL h5fclose_f(file1_id, error)

       CALL h5close_f(error)
       call cpu_time(t3)

       print*,t2-t1,t3-t2

end subroutine create_H46
subroutine create_H56(nv,nc,mf_v,mf_c,wcoul0,s,vol,p,n5,n6,chunk)
      double precision :: fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16

      integer, intent(in) :: nv,nc,s,p,n5,n6,chunk
      integer :: nb,i,j,alpha,beta,gama,eta,m,n,imatrix,keyword,c1,c2,c,x,A
      double precision            :: element
      double precision ::wcoul0,vol,a1,b1
      complex(kind=8),dimension(:,:),allocatable                            :: H56
      complex(kind=8), allocatable :: bse_mat_e(:,:,:,:)
      complex(kind=8), allocatable :: bse_mat_h(:,:,:,:)
      complex(kind=8), allocatable :: bse_mat(:,:,:,:)
      CHARACTER(LEN=11), PARAMETER :: filename1 = "singlet.h5" ! File name
      CHARACTER(LEN=5), PARAMETER :: dsetname1 = "CIS56"    ! Dataset name
      CHARACTER(LEN=11), PARAMETER :: filename2 = "singletX.h5" ! File name
      CHARACTER(LEN=5), PARAMETER :: dsetname2 = "HCISD"    ! Dataset name

      INTEGER, PARAMETER   :: DRANK = 3 ! Dataset rank 1 dimension for real or complex 2 is row-index 3 is col-index
      INTEGER(SIZE_T), PARAMETER :: NUMP = 2 ! Number of points selected
      complex(kind=8) ,DIMENSION(nv) ,intent(in) :: mf_v
      complex(kind=8),DiMENSION(nc),intent(in) :: mf_c
      INTEGER(HID_T) :: file1_id       ! File1 identifier
      INTEGER(HID_T) :: dset1_id       ! Dataset1 identifier
      INTEGER(HID_T) :: file2_id       ! File2 identifier
      INTEGER(HID_T) :: dset2_id       ! Dataset2 identifier
      INTEGER(HID_T) :: dataspace1     ! Dataspace identifier
          ! memspace identifier
      INTEGER(HID_T) :: dataspace2     ! Dataspace identifier
            ! memspace identifier
      INTEGER(HID_T) :: fapl_id

      INTEGER(HSIZE_T), DIMENSION(2) :: dimsm
                                                   ! Memory dataspace dimensions
      INTEGER(HSIZE_T), DIMENSION(2) :: dimsf
                                                   ! File dataspace dimensions
      INTEGER(HSIZE_T) :: offset(2),count(2) ! Elements coordinates
                                                      ! in the file
      DOUBLE PRECISION, DIMENSION(2) :: val =(/3.14,1.414/)  ! Values to write
      INTEGER :: memrank = 1  ! Rank of the dataset in memory
      DOUBLE PRECISION,DIMENSION(:,:,:),allocatable :: bufnew

      INTEGER :: error,info  ! Error flag
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims

      integer(HID_T) :: filespace     ! Dataspace identifier in file
      integer(HID_T) :: memspace      ! Dataspace identifier in mem
      integer(HID_T) :: dset_id       ! Dataset identifier
      INTEGER(HID_T) :: plist_id      ! Property list identifier

      !!!!!!!!!!!mpivariables!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer :: nproc, chunk_size,itr_length
      integer , dimension(:,:) , allocatable :: ij_mat
      integer, dimension(:),allocatable :: row_ij
      integer :: comm, rank, numproc, ierror
      integer :: cij,cab
      !call MPI_INIT(ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)


       keyword = 21
       fac_hd=0
       fac_hx=0
       fac_ed=0
       fac_ex=0
      fac_1=0
      fac_2=0
      fac_3=0
      fac_4=0
      fac_5=0
      fac_6=0
      fac_7=0
      fac_8=0
      fac_9=0
      fac_10=0
      fac_11=0
      fac_12=0
      fac_13=0
      fac_14=0
      fac_15=0
      fac_16=0

       itr_length = nc*(nc-1)/2
      chunk_size = chunk
      print*,chunk_size,nproc
      nb = nv + nc

      allocate(bse_mat_h(nv,nv,2,nv))
      allocate(bse_mat_e(nb,1,1,nb))
      allocate(row_ij(chunk_size))
      A = nv
      allocate(ij_mat(2,chunk_size))
      c = rank+1
      nrows_p= chunk_size*A
      allocate(H56(nrows_p,nc*nv))
      H56= cmplx(0.0,0.0)
      CALL h5screate_simple_f(1, dimsm, memspace, error)

      !print*,"proc_id+1",c
      call distribute_ij(rank,nc,A,chunk_size,ij_mat,row_ij)
      c1=0
      c2=0
      cij=0
      cab=0
       do imatrix = 3,4
          call fac_values(wcoul0,vol,keyword,imatrix,fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                                       fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16)
         do x = 1,chunk_size
            print*,x,ij_mat(1,x),ij_mat(2,x)

            i = ij_mat(1,x)
            j = ij_mat(2,x)
            cij = (x-1)*A
            cab =0
            if (imatrix==3) then
                call load_s_Hhh(i,j,nv,bse_mat_h)
             end if


            if (imatrix .eq. 3) then
                call load_s_Hee(i,j,nb,nv,bse_mat_e)
            end if

            do alpha = 1,nv
                  cab= cab+1
                  c1 = cij+cab
                  c2 = 0
                  do m = 1,nc
                     do  gama = 1,nv
                         c2 = c2 + 1
                         if (j==m) then
                            if (imatrix .ne. 4) then
                                     !print*,"bb",c1,c2, H(5,11),H(71,1),imatrix
                                      !print*,"1",H(c1,c2),bse_mat(eta,beta,alpha,gama), bse_mat(gama,beta,alpha,eta)
                               H56(c1,c2) =H56(c1,c2) + (fac_1 * bse_mat_h(alpha,gama,1,alpha)) &          !Fac 1234 stands for direct terms fac 5678 stands for exchange terms
                                       + (fac_5 * bse_mat_h(alpha,gama,1,alpha))
                            end if
                         end if
                         if (i==m) then
                           if (imatrix .ne. 4) then
                                     !print*,"bb",c1,c2, H(5,11),H(71,1),imatrix
                                      !print*,"1",H(c1,c2),bse_mat(eta,beta,alpha,gama), bse_mat(gama,beta,alpha,eta)
                              H56(c1,c2) =H56(c1,c2) + (fac_2 * bse_mat_h(alpha,gama,2,alpha)) &          !Fac 1234 stands for direct terms fac 5678 stands for exchange terms
                                       + (fac_6 * bse_mat_h(alpha,gama,2,alpha))
                            end if
                         end if
                         if (alpha == gama ) then
                           if (imatrix .ne. 4) then
                                     !print*,"bb",c1,c2, H(5,11),H(71,1),imatrix
                                      !print*,"1",H(c1,c2),bse_mat(eta,beta,alpha,gama), bse_mat(gama,beta,alpha,eta)
                             H56(c1,c2) =H56(c1,c2) + (fac_3 * bse_mat_e(alpha,1,1,m+nv)) &          !Fac 1234 stands for direct terms fac 5678 stands for exchange terms
                                       + (fac_7 * bse_mat_e(m+nv,1,1,alpha))
                             H56(c1,c2) =H56(c1,c2) + (fac_4 * bse_mat_e(alpha,1,1,m+nv)) &          !Fac 1234 stands for direct terms fac 5678 stands for exchange terms
                                       + (fac_8 * bse_mat_e(m+nv,1,1,alpha))
                           end if
                         end if

                      end do
                   end do
               end do
            end do
        end do
           call cpu_time(t2)
      CALL h5open_f(error)
      !CALL h5fcreate_f(filename1,  H5F_ACC_TRUNC_F, file1_id, error)


      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)

      !call h5fcreate_f(filename1, H5F_ACC_TRUNC_F, file1_id, error, &
       !           access_prp = plist_id)
      CALL h5fopen_f (filename1, H5F_ACC_RDWR_F, file1_id, error,access_prp = plist_id)


       !CALL h5fopen_f (filename2, H5F_ACC_RDWR_F, file2_id, error,access_prp = plist_id)
       CALL h5pclose_f(plist_id, error)
       dimsf(1) = p
       dimsf(2) = nc*nv
       CALL h5screate_simple_f(2, dimsf, dataspace1, error)

       CALL h5dcreate_f(file1_id, dsetname1,H5T_NATIVE_DOUBLE,dataspace1,dset1_id, error)
       CALL h5sclose_f(dataspace1, error)
     !
       !CALL h5dopen_f(file2_id, dsetname2, dset2_id, error)
     !
     ! Get  dataspace identifier.
     !
      count(1) = nrows_p
       count(2) = nc*nv
       dimsm(1)=nrows_p
       dimsm(2)=nc*nv
       dimsf(1) = p
       dimsf(2) = nc*nv
       !call h5screate_simple_f(2, dimsm, memspace, error)
       !allocate(data_out(count(1),count(2)))
       !data_out = real(H12)
       offset(1) = rank*nrows_p
       offset(2) = 0
       call h5screate_simple_f(2, count, memspace, error)
       CALL h5dget_space_f(dset1_id, dataspace1, error)
       CALL h5sselect_hyperslab_f (dataspace1, H5S_SELECT_SET_F, offset,&
               count, error)


       CALL h5pcreate_f(H5P_DATASET_XFER_F, fapl_id, error)
       CALL h5pset_dxpl_mpio_f(fapl_id, H5FD_MPIO_COLLECTIVE_F, error)
       call h5screate_simple_f(2, dimsm, memspace, error)
       CALL h5dwrite_f(dset1_id, H5T_NATIVE_DOUBLE, real(H56), dimsf , error, &
                    mem_space_id=memspace, file_space_id=dataspace1)
       call h5pclose_f(fapl_id, error)
       CALL h5sclose_f(memspace, error)
       !call h5pclose_f(fapl_id, error)

       CALL h5dclose_f(dset1_id, error)
       CALL h5fclose_f(file1_id, error)

       CALL h5close_f(error)
       call cpu_time(t3)

       print*,t2-t1,t3-t2


 end subroutine create_H56

end module create_H
