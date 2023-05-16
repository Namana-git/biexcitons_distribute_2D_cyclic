module hdf5_read_m

  
  use hdf5
  use global_variables
  include 'mpif.h'

  

  private

  public :: read_inp,load_mf,load_weff,bse_matrix,fac_values,load_Hee,load_Hhh,&
             distribute_ij,load_Heh,load_s_Hee,load_s_Hhh,calculate_chunk_size,&
             distrubute_rows_2D_cyclo_1,calculate_c2start_icol,load_larray, &
             distrubute_rows_2D_cyclo


            

contains


subroutine read_inp(nv,nc,vol)

  integer, intent(out) :: nc,nv
  double precision, intent(out) :: vol
  character(len=9), parameter :: filename = "bsemat.h5"  
  character(len=21), parameter :: dsetname_1 = "/bse_header/bands/nvb"
  character(len=21), parameter :: dsetname_2 = "/bse_header/bands/ncb"
  character(len=25), parameter :: dsetname_3 = "/mf_header/crystal/celvol"

  integer(hid_t) :: file_id  
  integer(hid_t) :: dset1_id,dset2_id,dset3_id,dset4_id
  integer     ::   error
  integer(hsize_t), DIMENSION(0) :: data1_dims,data2_dims,data3_dims,data4_dims
  
  
  
  call h5open_f(error)
  call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)
  
  call h5dopen_f(file_id, dsetname_1, dset1_id, error)
  call h5dread_f(dset1_id, H5T_NATIVE_INTEGER, nv, data1_dims, error)
  call h5dclose_f(dset1_id, error)
  
  call h5dopen_f(file_id, dsetname_2, dset2_id, error)
  call h5dread_f(dset2_id, H5T_NATIVE_INTEGER, nc, data2_dims, error)
  call h5dclose_f(dset2_id, error)  
  
  call h5dopen_f(file_id, dsetname_3, dset3_id, error)
  call h5dread_f(dset3_id, H5T_NATIVE_DOUBLE, vol, data3_dims, error)
  call h5dclose_f(dset3_id, error)
  
  call h5fclose_f(file_id, error)
  call h5close_f(error)
 
  !print*,"vol",vol
end subroutine read_inp

subroutine load_mf(nv,nc,mf_v,mf_c)

    integer, intent(in) :: nc,nv
    complex(kind=8) ,DIMENSION(nv) ,intent(out) :: mf_v
    complex(kind=8),DiMENSION(nc),intent(out) :: mf_c
    double precision :: zero,m,n
    character(len=9), parameter :: filename = "bsemat.h5"  
    character(len=25), parameter :: dsetname_1 = "/mf_header/kpoints/mnband"
    character(len=21), parameter :: dsetname_2 = "/mf_header/kpoints/el"
    integer(hid_t) :: file_id  
    integer(hid_t) :: dset1_id,dset2_id
    integer     ::   error,i,j,v,c,d
    integer :: mnband
    integer(hsize_t), DIMENSION(0) :: data1_dims
    double precision, DIMENSION(:,:,:), allocatable :: mf_b
    integer(hsize_t), DIMENSION(3) :: data2_dims
    
    !call h5open_f(error)
    !call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)
    
    !call h5dopen_f(file_id, dsetname_1, dset1_id, error)
    !call h5dread_f(dset1_id, H5T_NATIVE_INTEGER, mnband, data1_dims, error)
    !call h5dclose_f(dset1_id, error)
    
    !call h5dopen_f(file_id, dsetname_2, dset2_id, error)
    i = 1
    j = 1
    !allocate(mf_b(mnband,i,j))
    !call h5dread_f(dset2_id, H5T_NATIVE_DOUBLE, mf_b, data2_dims, error)
    !call h5dclose_f(dset2_id, error)
    !call h5fclose_f(file_id, error)
    !call h5close_f(error)
    zero = 0.0000000
    open(1, file = 'valence states', status = 'old') 
    do v = 1,nv
       read(1,*) m
       mf_v(v) = cmplx(m, zero)/13.6056980659 
       !print*,"valence",v,mf_v(v)
    end do 
    close(1)
    open(2, file = 'conduction states', status = 'old')
    do c = 1,nc
       read(2,*) n
       mf_c(c) = cmplx(n, zero)/13.6056980659 
       !print*,"conduction",c,mf_c(c)
    end do 
    close(2)
   ! deallocate(mf_b)
end subroutine load_mf

subroutine load_weff(wcoul0)
  double precision :: a,b,c,weff
  integer :: d,e,f
  double precision, intent(out) :: wcoul0
  open (2, file = 'vcoul', status = 'old')
  read(2,*) a,b,c,d,e,f,weff
  close(2)
  wcoul0 = weff
  


end subroutine load_weff
!This subroutine broadcasts H1 to all the processors
!Then the processors that contains blocks of elements of grow, copies them at the corresponding local_indices 
subroutine distrubute_rows_2D_cyclo(H1,grow_size,c1,c2_start,x,chunk_size)
       integer , intent(in) :: c1,grow_size
       integer , intent(in):: c2_start,x,chunk_size

       complex(kind=8) , dimension(grow_size), intent(in) :: H1
       complex(kind=8) , dimension(grow_size)             :: grow
       integer                                            ::  ierr,i,j,k,grow_row_index,c,l
       integer                                            :: icol_start,icol_contrib
       integer                                            :: ic2_start,ic2_start_block_size
       integer                                            :: jastart,jaend
       integer                                            :: lrindx,lcindx   
       integer                                            :: ipos
       integer                                            :: copy_me,root  
       grow = cmplx(0.0,0.0)
       grow_row_index = 1
       !ipos = 2  !debug statement
       !now each col proc has its own starting global col_index
       !some col procs may not have any elements prsent within H1
       !such procs has icol_contrib set to zero 
       
       
       do i = 0,(mpi%size_-1) !loop over all the processors
          copy_me = 0
        
          
         ! initialy grow is a zero array
         ! we need to load it with H1 array of ith proc
         !But if x .ge. chunk_size the grow should remain zero
         !in that case all the processors should have read flag set to zero

           call MPI_Barrier(MPI_COMM_WORLD, error)
         if(mpi%rank == i)then 
           if(x < chunk_size+1)then
             grow = H1
             grow_row_index = c1 
             copy_me = 1
             !do l = 1, grow_size
               ! print*,"i,grow_row_index,l ,grow(l)****",i,grow_row_index,l,grow(l)
             !end do

           end if
         end if
           call MPI_Barrier(MPI_COMM_WORLD, error)

         !Now grow of the ith processor contains H1 that has to be broad casted
        ! also the row index is broad casted
        ! Also if the grow is a valid H1 or zero grow in case where x .ge. chunk_size 
        
          call MPI_BCAST(grow_row_index, 1, MPI_INTEGER, i, MPI_COMM_WORLD, ierr)
          call MPI_BCAST(copy_me, 1, MPI_INTEGER, i, MPI_COMM_WORLD, ierr)
          call MPI_Barrier(MPI_COMM_WORLD, error)
           CALL MPI_BCAST(grow,grow_size,MPI_DOUBLE_COMPLEX, i, MPI_COMM_WORLD, ierr)

           call MPI_Barrier(MPI_COMM_WORLD, error)

         
        
           do l = 1, grow_size
              CALL MPI_BCAST(grow(l), 1, MPI_COMPLEX, i, MPI_COMM_WORLD, ierr)
           end do
          
          
           !if(mpi%rank == 3)then
            
           !      do l = 1, grow_size
            !         print*,"grow_row_index,l ,grow(l)",grow_row_index,l,grow(l)
            !     end do
             
          ! end if

          call MPI_Barrier(MPI_COMM_WORLD, error)

          !only if valid H1 is copied   
 
         if (copy_me == 1) then
             
                 
             ! the myprow that contains the row broadcasted
             irow = indxg2p(grow_row_index,pzheevx_vars%mb,0,0,grid%nprow)

             ! the starting mypcol that corresponds to c2_start
             icol_start = indxg2p( c2_start,pzheevx_vars%nb,0,0,grid%npcol)
             ! calculate starting global indices and starting block sizes for each of the mypcol
             call calculate_c2start_icol(grid%mypcol,icol_start,c2_start,grow_size,ic2_start,&
                icol_contrib,ic2_start_block_size)
             


          if (icol_contrib == 1) then

            if (grid%myprow == irow) then 

               !jastart and jaend are the global indices of elements belonging to that processor 
               jastart = ic2_start                                                          
               
              
               jaend = ic2_start+ic2_start_block_size-1 
               lrindx = INDXG2L(grow_row_index,pzheevx_vars%nb, 0, 0,grid%npcol)  
               
            
               do k = ic2_start,c2_start+grow_size-1,(grid%nprow)*pzheevx_vars%nb 
                  do j = jastart,jaend !copy elements within the row block 
               
                     lcindx = INDXG2L(j,pzheevx_vars%mb, 0, 0,grid%nprow) 
                     ipos = lrindx + (lcindx-1)*hamiltonian%lld
                     hamiltonian%mat(ipos) = real(grow(j-c2_start+1))
                     if(grid%mypcol == 1 .and. grid%myprow ==0)then
                       ! print*," hamiltonian%mat(ipos),ipos", hamiltonian%mat(ipos),ipos
                     end if
                    
                     
                 end do
                 jastart = jastart + ((grid%nprow)*pzheevx_vars%nb) 
                 jaend= min(jastart+pzheevx_vars%mb-1,c2_start+grow_size-1)
              end do

             


          
           end if 
         end if
         !call MPI_Barrier(MPI_COMM_WORLD, error)
        end if
     
        call MPI_Barrier(MPI_COMM_WORLD, error)

       end do
end subroutine distrubute_rows_2D_cyclo



subroutine distrubute_rows_2D_cyclo_1(H1,grow_size,c1,c2_start)

    
    integer , intent(in) :: c1,grow_size
    integer , intent(in):: c2_start
    complex(kind=8) , dimension(grow_size), intent(in) :: H1 
    
    complex(kind=8), dimension(:), allocatable :: larray
    integer              :: i,iblock,r,jastart,jaend
    integer              :: nblocks
    integer              :: irow,icol,icol_start
    integer              :: ia,ja_start
    integer              :: ic2_start
    integer              :: icol_contrib,ic2_start_block_size
    integer              :: larray_size
    nblocks = sys_var%n1/pzheevx_vars%mb  

    do i = 0,(mpi%size_-1) !row index is fixed with i
       
       
      
       do icol = 0, (grid%npcol -1)

          

           if(mpi%rank == i)then
             irow = indxg2p( c1,pzheevx_vars%mb,0,0,grid%nprow)  
             ia =  INDXG2L( c1,pzheevx_vars%mb,0,0,grid%nprow)
             icol_start = indxg2p( c2_start,pzheevx_vars%nb,0,0,grid%npcol)
             call calculate_c2start_icol(icol,icol_start,c2_start,grow_size,ic2_start,&
             icol_contrib,ic2_start_block_size)
             
             if (icol_contrib == 1) then
                allocate(larray(larray_size))
                call load_larray(H1,larray,larray_size,c2_start,grow_size,ic2_start,&
                ic2_start_block_size,icol)
                deallocate(larray)
                

            
             !mpi send  


             end if
          endif
          if (grid%mypcol == icol .and. grid%myprow == irow)then

              if (icol_contrib == 1) then
              ! recieve
              !copy 
              end if
          end if
       end do 
              
    end do
       
 
    
    



end subroutine distrubute_rows_2D_cyclo_1


subroutine calculate_c2start_icol(icol,icol_start,c2_start,grow_size,ic2_start,icol_contrib,ic2_start_block_size)
          integer, intent(in) :: icol,icol_start,c2_start,grow_size
          integer, intent(inout) :: ic2_start,icol_contrib,ic2_start_block_size
          integer                :: iblock
          integer                :: i_size
          integer                :: r,rp
          integer                :: start,stop_,step
          

          r = mod(c2_start,pzheevx_vars%nb)

          !icol_start is the col proc that that contains c2_start
          !icol_start block are the block of consecutive elements that belong to proc icol_start
          i_size = 0  !i_size calculates the size of the local array

          if (r .eq. 0) then
              icol_start_block_size = 1
          else
              icol_start_block_size = pzheevx_vars%nb - r +1
          end if
          if (icol==icol_start) then
              ic2_start = c2_start
              ic2_start_block_size =  icol_start_block_size
              i_size =  icol_start_block_size
          end if
          
          if (icol .ne. icol_start) then
               if ( icol > icol_start) then
                   ic2_start = (icol - icol_start -1) * pzheevx_vars%nb + (icol_start_block_size + c2_start)
                   rp = mod(ic2_start,pzheevx_vars%nb)
                   ic2_start_block_size = pzheevx_vars%nb - rp + 1


               end if
               if ( icol < icol_start) then
                   ic2_start = ((grid%npcol-1)-icol_start + icol ) * pzheevx_vars%nb + (icol_start_block_size + c2_start  )
                   rp = mod(ic2_start,pzheevx_vars%nb)
                   ic2_start_block_size = pzheevx_vars%nb - rp + 1

               end if
 
              
          end if
          i_size =  ic2_start_block_size 
          !print*,"icol_start",icol_start
          !print*,"icol,ic2_start",icol,ic2_start

          if (ic2_start > (c2_start + grow_size-1)) then
             icol_contrib = 0
             larray_size = 0
          else
             icol_contrib = 1

          end if

         !print*,"icol,icol_contrib",icol,icol_contrib

          ! calculate the starting block size

          ! print*,"icol,ic2_start_block_size",icol,ic2_start_block_size



          ! calculate size of the local_array if icol_contrib = 1
         

                  



          
         !

end subroutine calculate_c2start_icol 

subroutine load_larray(H1,larray,larray_size,c2_start,grow_size,ic2_start,ic2_start_block_size,icol)
    integer , intent(in) ::larray_size,grow_size,ic2_start_block_size  
    integer , intent(in) :: c2_start,ic2_start,icol
    complex(kind=8) , dimension(grow_size), intent(in) :: H1
    complex(kind=8) , dimension(larray_size), intent(inout) :: larray
    integer              :: i,ip,l_i,g_i
    integer              :: start,stop_,step
    integer              :: g_start,g_end,l_start,l_end
    integer              :: iblock
    !load the first block which may not be a complete block
    
    do i = 1,ic2_start_block_size
       l_i = i
       g_i = ic2_start - c2_start+i
       larray(l_i) = H1(g_i)
       !print*,"icol,l_i,g_i",icol,l_i,g_i

    end do
    !print*,"ic2_start,icol",ic2_start,icol
   

    !load the remaining blocs of the local array
     start =(ic2_start-(pzheevx_vars%nb-ic2_start_block_size) +((grid%nprow)*pzheevx_vars%nb))   !second block start
     stop_ = (c2_start+grow_size-1)   !last index of the row 
     step = ((grid%nprow)*pzheevx_vars%nb)   !
     !print*,"icol, start,stop_,step",icol, start,stop_,step
     ip = ic2_start_block_size

    do iblock = start,stop_,step
      
                !print*,"icol,iblock",icol,iblock
                !print*,"hey",icol,iblock+pzheevx_vars%nb-1,(c2_start+grow_size -1)
                !print*,"2",icol, start,stop_,step

       ! if block resides inside the grow ( not the last block, 
       ! the no of elements may be less the the no of elements in a block 
       if((iblock+pzheevx_vars%nb-1) .le. (c2_start+grow_size-1))then
                        g_start = iblock - (c2_start -1)
                        g_end = g_start+pzheevx_vars%nb -1
                        l_start = ip+1
                        l_end = l_start+pzheevx_vars%nb-1
                        larray(l_start:l_end) = H1(g_start:g_end)
                        ip = l_end
                       ! print*,"g_start,g_end,l_start,l_end",g_start,g_end,l_start,l_end
       end if
       !if its the last block
       if ((iblock+pzheevx_vars%nb-1) > (c2_start+grow_size-1))then

                    g_start = iblock - (c2_start -1)
                    g_end = grow_size
                    l_start = ip+1
                    l_end = larray_size
                    larray(l_start:l_end) = H1(g_start:g_end)
                    
                    !print*,"g_start,g_end,l_start,l_end",g_start,g_end,l_start,l_end


       end if
       !print*,"icol,g_start,g_end,l_start,l_end",icol,g_start,g_end,l_start,l_end 


   end do
   
    print*,"larray,@,icol",larray,"@",icol



end subroutine load_larray

subroutine gather_rows_test()
      integer         :: i,ierror
      complex(kind=8) :: temp_row

      temp_row = cmplx(0.0,0.0)    

      !do i = 2,grid%nprow
          if (mpi%rank == 1)then
              temp_row = cmplx(10,0.2)
              call MPI_SEND(temp_row,1, MPI_COMPLEX, 0, 1, MPI_COMM_WORLD, ierror)
              print*,"send",ierror

          end if
          if (mpi%rank == 0)then
             call MPI_RECV(temp_row,1, MPI_COMPLEX, 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE,ierror)
             print*,"recv",ierror
            ! H_b(i,:) = temp_row(:)
          end if
          
          !print*,mpi%rank,temp_row
  
      !end do




end subroutine gather_rows_test

subroutine calculate_chunk_size(chunk_size)

               integer , intent(inout)  :: chunk_size
               integer         :: tot_length,c,r


               tot_length =sys_var%nc*(sys_var%nc-1)/2

               c = tot_length/mpi%size_
               r =  mod(tot_length,mpi%size_)
               !print*,"c,r,grid%nprow,tot_length",c,r,grid%nprow,tot_length
               if (mpi%rank < r) then
                
                  chunk_size = c + 1
                 ! print*,1,grid%myprow,chunk_size
               else
                  chunk_size = c
                !  print*,1,grid%myprow,chunk_size


               endif

end subroutine calculate_chunk_size



subroutine distribute_ij(rank,nc,A,chunk_size,ij_mat,row_ij)
   integer, intent(in) :: rank,nc,A,chunk_size
   integer            :: c1,c2,c3,i,j,r,l
   integer, dimension(2,chunk_size) , intent(inout) :: ij_mat
   integer, dimension(chunk_size) , intent(inout) :: row_ij
   c1 = 0
   c2 = 0
   c3 = 0
   l = nc*(nc-1)/2 
   r = mod(l,mpi%size_) 
 
    do j = 1, nc-1
       do i = j+1, nc
           c1=c1+1
           !print*,"ij",i,j
           if ( rank < r)then
              if ((c1 .ge. ((rank*chunk_size)+1)) .and. (c1 .le.((rank+1)*chunk_size))) then
                c2 = c2+1
                ij_mat(1,c2) = i
                ij_mat(2,c2) = j
                row_ij(c2)=c3*A
                print*,"before",ij_mat(1,c2),ij_mat(2,c2),c2
              end if

           end if
           if ( rank .ge.r) then 
               if ((c1 .ge. ((rank*chunk_size)+1+r)) .and. (c1 .le.(((rank+1)*chunk_size)+r))) then  
                  c2 = c2+1
                  ij_mat(1,c2) = i
                  ij_mat(2,c2) = j
                  row_ij(c2)=c3*A
               end if

           end if

           c3 = c3+1
       end do
    end do



end subroutine distribute_ij

subroutine load_Hee(i,j,nc,nv,bse_mat_e)
   integer, intent(in) :: i,j,nv,nc
   complex(kind=8), allocatable, intent(out) :: bse_mat_e(:,:,:,:)
   integer(HID_T) :: dataspace     ! Dataspace identifier in file
   integer(HID_T) :: memspace      ! Dataspace identifier in mem
   integer(HID_T) :: dset_id
   integer(HID_T) :: file_id
   integer(HSIZE_T) :: count(7), offset(7)
   integer :: rank, error
   double precision, allocatable :: data_out(:,:,:,:,:,:,:)
   character(len=9), parameter :: filename = "bsemat.h5"
   integer :: m,n
   allocate(bse_mat_e(nc,1,1,nc))
   CALL h5open_f(error)
   call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)
   call h5dopen_f(file_id, 'mats/body', dset_id, error)
   print*,nv,nc
   
   rank = 7
   count(1) = 2
   count(2) =nc
   count(3) = 1
   count(4) = 1
   count(5) = nc
   count(6) = 1
   count(7) = 1

   CALL h5dget_space_f(dset_id, dataspace, error)
   call h5screate_simple_f(rank, count, memspace, error)

   allocate(data_out(count(1),count(2),count(3),count(4),count(5),count(6),count(7)))


   offset(1)=0
   offset(2)=nv
   offset(3)=j+nv-1
   offset(4)=i+nv-1
   offset(5)=nv
   offset(6)=0
   offset(7)=0

   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                offset, count, error)
   call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data_out, count, error,&
         memspace, dataspace)
   bse_mat_e(:,1,1,:) = cmplx(data_out(1,:,1,1,:,1,1),data_out(2,:,1,1,:,1,1))
   !print*,"hey",i,j,(offset(2)+count(2)),(offset(3)+count(3)),(offset(4)+count(4)),(offset(5)+count(5)),bse_mat_e(:,1,1,:)
   
   deallocate(data_out) 
   


end subroutine load_Hee
subroutine load_s_Hee(i,j,nb,nv,bse_mat_e)
   integer, intent(in) :: i,j,nb,nv
   complex(kind=8), allocatable, intent(out) :: bse_mat_e(:,:,:,:)
   integer(HID_T) :: dataspace     ! Dataspace identifier in file
   integer(HID_T) :: memspace      ! Dataspace identifier in mem
   integer(HID_T) :: dset_id
   integer(HID_T) :: file_id
   integer(HSIZE_T) :: count(7), offset(7)
   integer :: rank, error
   double precision, allocatable :: data_out(:,:,:,:,:,:,:)
   character(len=9), parameter :: filename = "bsemat.h5"
   integer :: m,n
   allocate(bse_mat_e(nb,1,1,nb))
   CALL h5open_f(error)
   call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)
   call h5dopen_f(file_id, 'mats/body', dset_id, error)
  

   rank = 7
   count(1) = 2
   count(2) =nb
   count(3) = 1
   count(4) = 1
   count(5) = nb
   count(6) = 1
   count(7) = 1

   CALL h5dget_space_f(dset_id, dataspace, error)
   call h5screate_simple_f(rank, count, memspace, error)

   allocate(data_out(count(1),count(2),count(3),count(4),count(5),count(6),count(7)))


   offset(1)=0
   offset(2)=0
   offset(3)=j+nv-1
   offset(4)=i+nv-1
   offset(5)=0
   offset(6)=0
   offset(7)=0

   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                offset, count, error)
   call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data_out, count, error,&
         memspace, dataspace)
   bse_mat_e(:,1,1,:) = cmplx(data_out(1,:,1,1,:,1,1),data_out(2,:,1,1,:,1,1))
   deallocate(data_out)



end subroutine load_s_Hee



subroutine load_Hhh(nv,bse_mat_h)
   integer, intent(in) ::nv
   complex(kind=8), allocatable, intent(out) :: bse_mat_h(:,:,:,:)
   integer(HID_T) :: dataspace     ! Dataspace identifier in file
   integer(HID_T) :: memspace      ! Dataspace identifier in mem
   integer(HID_T) :: dset_id
   integer(HID_T) :: file_id
   integer(HSIZE_T) :: count(7), offset(7)
   integer :: rank, error
   double precision, allocatable :: data_out(:,:,:,:,:,:,:)
   character(len=9), parameter :: filename = "bsemat.h5"
   
   allocate(bse_mat_h(nv,nv,nv,nv))
   CALL h5open_f(error)
   call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)
   call h5dopen_f(file_id, 'mats/body', dset_id, error)
   print*,nv

   rank = 7
   count(1) = 2
   count(2) =nv
   count(3) = nv
   count(4) = nv
   count(5) = nv
   count(6) = 1
   count(7) = 1

   CALL h5dget_space_f(dset_id, dataspace, error)
   call h5screate_simple_f(rank, count, memspace, error)

   allocate(data_out(count(1),count(2),count(3),count(4),count(5),count(6),count(7)))


   offset(1)=0
   offset(2)=0
   offset(3)=0
   offset(4)=0
   offset(5)=0
   offset(6)=0
   offset(7)=0

   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                offset, count, error)
   call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data_out, count, error,&
         memspace, dataspace)
   bse_mat_h(:,:,:,:) = cmplx(data_out(1,:,:,:,:,1,1),data_out(2,:,:,:,:,1,1))

   deallocate(data_out)

  ! print *,bse_mat_h(:,:,:,:)




end subroutine load_Hhh
subroutine load_s_Hhh(i,j,nv,bse_mat_h)
   integer, intent(in) ::nv,i,j
   complex(kind=8), allocatable, intent(out) :: bse_mat_h(:,:,:,:)
   integer(HID_T) :: dataspace     ! Dataspace identifier in file
   integer(HID_T) :: memspace      ! Dataspace identifier in mem
   integer(HID_T) :: dset_id
   integer(HID_T) :: file_id
   integer(HSIZE_T) :: count(7), offset(7)
   integer :: rank, error
   double precision, allocatable :: data_out(:,:,:,:,:,:,:)
   character(len=9), parameter :: filename = "bsemat.h5"

   allocate(bse_mat_h(nv,nv,2,nv))
   CALL h5open_f(error)
   call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)
   call h5dopen_f(file_id, 'mats/body', dset_id, error)
   print*,nv

   rank = 7
   count(1) = 2
   count(2) =nv
   count(3) = nv
   count(4) = 1
   count(5) = nv
   count(6) = 1
   count(7) = 1

   CALL h5dget_space_f(dset_id, dataspace, error)
   call h5screate_simple_f(rank, count, memspace, error)

   allocate(data_out(count(1),count(2),count(3),count(4),count(5),count(6),count(7)))


   offset(1)=0
   offset(2)=0
   offset(3)=0
   offset(4)=i+nv-1
   offset(5)=0
   offset(6)=0
   offset(7)=0

   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                offset, count, error)
   call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data_out, count, error,&
         memspace, dataspace)
  
   bse_mat_h(:,:,1,:) = cmplx(data_out(1,:,:,1,:,1,1),data_out(2,:,:,1,:,1,1))

   deallocate(data_out)
   rank = 7
   count(1) = 2
   count(2) =nv
   count(3) = nv
   count(4) = 1
   count(5) = nv
   count(6) = 1
   count(7) = 1

   CALL h5dget_space_f(dset_id, dataspace, error)
   call h5screate_simple_f(rank, count, memspace, error)

   allocate(data_out(count(1),count(2),count(3),count(4),count(5),count(6),count(7)))


   offset(1)=0
   offset(2)=0
   offset(3)=0
   offset(4)=j+nv-1
   offset(5)=0
   offset(6)=0
   offset(7)=0

   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                offset, count, error)
   call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data_out, count, error,&
         memspace, dataspace)
   bse_mat_h(:,:,2,:) = cmplx(data_out(1,:,:,1,:,1,1),data_out(2,:,:,1,:,1,1))

   deallocate(data_out)

end subroutine load_s_Hhh


subroutine load_Heh(i,j,nc,nv,imatrix,bse_mat)
   integer, intent(in) :: i,j,nc,nv,imatrix
   complex(kind=8), allocatable, intent(out) :: bse_mat(:,:,:,:)
   integer(HID_T) :: dataspace     ! Dataspace identifier in file
   integer(HID_T) :: memspace      ! Dataspace identifier in mem
   integer(HID_T) :: dset_id
   integer(HID_T) :: file_id
   integer(HSIZE_T) :: count(7), offset(7)
   integer :: rank, error
   double precision, allocatable :: data_out(:,:,:,:,:,:,:)
   character(len=9), parameter :: filename = "bsemat.h5"

   allocate(bse_mat(nv,nv,2,nc))
   CALL h5open_f(error)
   call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)
   if (imatrix .eq. 3) then
        call h5dopen_f(file_id, 'mats/body', dset_id, error)
   else if (imatrix .eq. 4) then
        call h5dopen_f(file_id, 'mats/exchange', dset_id, error)
   end if

   print*,"imatrix",imatrix

   rank = 7
   count(1) = 2
   count(2) =nv
   count(3) = nv
   count(4) = 1
   count(5) = nc
   count(6) = 1
   count(7) = 1

   CALL h5dget_space_f(dset_id, dataspace, error)
   call h5screate_simple_f(rank, count, memspace, error)

   allocate(data_out(count(1),count(2),count(3),count(4),count(5),count(6),count(7)))


   offset(1)=0
   offset(2)=0
   offset(3)=0
   offset(4)=i+nv-1
   offset(5)=nv
   offset(6)=0
   offset(7)=0

   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                offset, count, error)
   call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data_out, count, error,&
         memspace, dataspace)
   bse_mat(:,:,1,:) = cmplx(data_out(1,:,:,1,:,1,1),data_out(2,:,:,1,:,1,1))
   offset(1)=0
   offset(2)=0
   offset(3)=0
   offset(4)=j+nv-1
   offset(5)=nv
   offset(6)=0
   offset(7)=0

   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                offset, count, error)
   call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data_out, count, error,&
         memspace, dataspace)
   bse_mat(:,:,2,:) = cmplx(data_out(1,:,:,1,:,1,1),data_out(2,:,:,1,:,1,1))

   deallocate(data_out)



end subroutine load_Heh 
subroutine bse_matrix(imatrix,bse_mat,nb)
  integer, intent(in) :: imatrix,nb
  complex(kind=8),dimension(nb,nb,nb,nb),intent(inout) :: bse_mat
  character(len=9), parameter :: filename = "bsemat.h5"  
  character(len=10), parameter :: name_1 = "/mats/head"
  character(len=10), parameter :: name_2 = "/mats/wing"
  character(len=10), parameter :: name_3 = "/mats/body"
  character(len=14), parameter :: name_4 = "/mats/exchange"
  double precision, DIMENSION(:,:,:,:,:,:,:), allocatable :: data1_out

  integer(hid_t) :: file_id  
  integer(hid_t) :: dset1_id,dset2_id
  integer     ::   error,iv,ivp,ic,icp
  integer(hsize_t) :: data1_dims(7)  
  
  
  !call h5open_f(error)
  call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)
  if (imatrix== 1) then
     call h5dopen_f(file_id, name_1, dset1_id, error)
  elseif (imatrix== 2) then
     call h5dopen_f(file_id, name_2, dset1_id, error)
  elseif (imatrix== 3) then
     call h5dopen_f(file_id, name_3, dset1_id, error)
  elseif (imatrix== 4) then
     call h5dopen_f(file_id, name_4, dset1_id, error)
  end if
  allocate(data1_out(2,nb,nb,nb,nb,1,1))
  call h5dread_f(dset1_id, H5T_NATIVE_DOUBLE, data1_out, data1_dims, error)
  call h5dclose_f(dset1_id, error)
  do icp=1,nb
     do ic=1,nb
        do ivp=1,nb
            do iv=1,nb
                bse_mat(iv,ivp,ic,icp)=cmplx(data1_out(1,iv,ivp,ic,icp,1,1),data1_out(2,iv,ivp,ic,icp,1,1))
                !print*,bse_mat(iv,ivp,ic,icp),iv,ivp,ic,icp
            enddo
        enddo
     enddo
  enddo
  !print*,bse_mat(1,1,2,1),imatrix
  call h5fclose_f(file_id, error)
  !call h5close_f(error)
  deallocate(data1_out)
end subroutine bse_matrix

subroutine fac_values(wcoul0,vol,keyword,imatrix,fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                         fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16)     
  double precision, intent(in):: wcoul0,vol
  integer, intent(in) :: keyword,imatrix
  double precision, intent(inout) :: fac_ed,fac_ex,fac_hd,fac_hx,fac_1,fac_2,fac_3, &
                                     fac_4,fac_5,fac_6,fac_7,fac_8,fac_9,fac_10,fac_11,fac_12,fac_13,fac_14,fac_15,fac_16
  double precision :: w,f,r3t,r2t
  w = 0
  f = 2 
  !print*,"fac",f,w
  r3t = 1.732050808
  r2t = 1.414213562
  if (keyword==1) then !singlet
     if (imatrix==1) then
        fac_ed = w
        fac_ex = -w
        fac_hd = w
        fac_hx = -w
        fac_1 = w
        fac_2 = -w 
        fac_3 = -w
        fac_4 = w
        fac_5 = w
        fac_6 = -w
        fac_7= -w
        fac_8 = w
        fac_9 = -w
        fac_10 = w
        fac_11 = w
        fac_12 = -w
        fac_13 = -w
        fac_14 = w
        fac_15 = w
        fac_16 = -w
     elseif (imatrix == 4) then 
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = -1.5*f
        fac_2 = 1.5*f
        fac_3 = 1.5*f
        fac_4 = -1.5*f
        fac_5 = -1.5*f
        fac_6 = 1.5*f
        fac_7= 1.5*f
        fac_8 = -1.5*f
        fac_9 = 1.5*f
        fac_10 = -1.5*f
        fac_11 = -1.5*f
        fac_12 = 1.5*f
        fac_13 = 1.5*f
        fac_14 = -1.5*f
        fac_15 = -1.5*f
        fac_16 = 1.5*f
     else
        fac_ed = f
        fac_ex = -f
        fac_hd = f
        fac_hx = -f
        fac_1 = f
        fac_2 = -f 
        fac_3 = -f
        fac_4 = f
        fac_5 = f
        fac_6 = -f
        fac_7= -f
        fac_8 = f
        fac_9 = -f
        fac_10 = f
        fac_11 = f
        fac_12 = -f
        fac_13 = -f
        fac_14 = f
        fac_15 = f
        fac_16 = -f
     end if
  endif
  if (keyword==2) then !singlet2
     if (imatrix==1) then
        fac_ed = w
        fac_ex = w
        fac_hd = w
        fac_hx = w
        fac_1 = -w
        fac_2 = -w 
        fac_3 = -w
        fac_4 = -w
        fac_5 = -w
        fac_6 = -w
        fac_7= -w
        fac_8 = -w
        fac_9 = -w
        fac_10 = -w
        fac_11 = -w
        fac_12 = -w
        fac_13 = -w
        fac_14 = -w
        fac_15 = -w
        fac_16 = -w
     elseif (imatrix == 4) then 
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = 0.5*f
        fac_2 = 0.5*f
        fac_3 = 0.5*f
        fac_4 = 0.5*f
        fac_5 = 0.5*f
        fac_6 =0.5*f
        fac_7= 0.5*f
        fac_8 = 0.5*f
        fac_9 = 0.5*f
        fac_10 = 0.5*f
        fac_11 = 0.5*f
        fac_12 = 0.5*f
        fac_13 = 0.5*f
        fac_14 = 0.5*f
        fac_15 = 0.5*f
        fac_16 = 0.5*f
     else
        fac_ed = f
        fac_ex = f
        fac_hd = f
        fac_hx = f
        fac_1 = -f
        fac_2 = -f 
        fac_3 = -f
        fac_4 = -f
        fac_5 = -f
        fac_6 = -f
        fac_7= -f
        fac_8 = -f
        fac_9 = -f
        fac_10 = -f
        fac_11 = -f
        fac_12 = -f
        fac_13 = -f
        fac_14 = -f
        fac_15 = -f
        fac_16 = -f
    end if
  endif
  if (keyword==3) then !singlet
     if (imatrix==1) then
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = 0
        fac_2 = 0 
        fac_3 = 0
        fac_4 = 0
        fac_5 = 0
        fac_6 = 0
        fac_7= 0
        fac_8 = 0
        fac_9 = 0
        fac_10 = 0
        fac_11 = 0
        fac_12 = 0
        fac_13 = 0
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0
     elseif (imatrix == 4) then 
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = -r3t*f/2
        fac_2 = -r3t*f/2
        fac_3 = r3t*f/2
        fac_4 = r3t*f/2
        fac_5 = r3t*f/2
        fac_6 = r3t*f/2
        fac_7= -r3t*f/2
        fac_8 = -r3t*f/2
        fac_9 = -r3t*f/2
        fac_10 = -r3t*f/2
        fac_11 = r3t*f/2
        fac_12 = r3t*f/2
        fac_13 = r3t*f/2
        fac_14 = r3t*f/2
        fac_15 = -r3t*f/2
        fac_16 = -r3t*f/2
     else
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = 0
        fac_2 = 0 
        fac_3 = 0
        fac_4 = 0
        fac_5 = 0
        fac_6 = 0
        fac_7= 0
        fac_8 = 0
        fac_9 = 0
        fac_10 = 0
        fac_11 = 0
        fac_12 = 0
        fac_13 = 0
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0

    end if
  endif
  if (keyword==4) then !singlet_double
    if (imatrix==1) then
        fac_ed = w
        fac_ex = 0
        fac_hd = w
        fac_hx = 0
        fac_1 = 0
        fac_2 = 0 
        fac_3 = 0
        fac_4 = 0
        fac_5 = 0
        fac_6 = 0
        fac_7= 0
        fac_8 = 0
        fac_9 = -w
        fac_10 = 0
        fac_11 = 0
        fac_12 = -w
        fac_13 = -w
        fac_14 = 0
        fac_15 = 0
        fac_16 = -w
    elseif (imatrix == 4) then 
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = 0
        fac_2 = 0
        fac_3 = 0
        fac_4 = 0
        fac_5 = 0
        fac_6 = 0
        fac_7=  0
        fac_8 = 0
        fac_9 = 0 
        fac_10 = 0
        fac_11 = 0
        fac_12 = f
        fac_13 = f
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0
    else
        fac_ed = f
        fac_ex = 0
        fac_hd = f
        fac_hx = 0
        fac_1 = 0
        fac_2 = 0 
        fac_3 = 0
        fac_4 = 0
        fac_5 = 0
        fac_6 = 0
        fac_7= 0
        fac_8 = 0
        fac_9 = -f
        fac_10 = 0
        fac_11 = 0
        fac_12 = -f
        fac_13 = -f
        fac_14 = 0
        fac_15 = 0
        fac_16 = -f
    end if
  endif
  if (keyword==5) then !H13
     if (imatrix==1) then
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = w/(r3t*2)
        fac_2 = -w/(r3t*2) 
        fac_3 = -w/(r3t*2)
        fac_4 = w/(r3t*2)
        fac_5 = w/(r3t*2)
        fac_6 = -w/(r3t*2)
        fac_7=  -w/(r3t*2)
        fac_8 =  w/(r3t*2)
        fac_9 =  -w/(r3t*2)
        fac_10 = w/(r3t*2)
        fac_11 = w/(r3t*2)
        fac_12 = -w/(r3t*2)
        fac_13 = -w/(r3t*2)
        fac_14 = w/(r3t*2)
        fac_15 = w/(r3t*2)
        fac_16 = -w/(r3t*2)
     elseif (imatrix == 4) then 
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = -r3t*f/2
        fac_2 = 0
        fac_3 = r3t*f/2
        fac_4 = 0
        fac_5 = 0
        fac_6 = r3t*f/2
        fac_7=  0
        fac_8 = -r3t*f/2
        fac_9 = 0
        fac_10 = -r3t*f/2
        fac_11 = 0
        fac_12 = r3t*f/2
        fac_13 = r3t*f/2
        fac_14 = 0
        fac_15 = -r3t*f/2
        fac_16 = 0
     else
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = f/(r3t*2)
        fac_2 = -f/(r3t*2) 
        fac_3 = -f/(r3t*2)
        fac_4 = f/(r3t*2)
        fac_5 = f/(r3t*2)
        fac_6 = -f/(r3t*2)
        fac_7=  -f/(r3t*2)
        fac_8 =  f/(r3t*2)
        fac_9 = -f/(r3t*2)
        fac_10 = f/(r3t*2)
        fac_11 = f/(r3t*2)
        fac_12 = -f/(r3t*2)
        fac_13 = -f/(r3t*2)
        fac_14 = f/(r3t*2)
        fac_15 = f/(r3t*2)
        fac_16 = -f/(r3t*2)
     end if
   end if
   if (keyword==6) then !H23
     if (imatrix==1) then
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = -w/2
        fac_2 = -w/2 
        fac_3 =  -w/2 
        fac_4 = -w/2
        fac_5 = -w/2
        fac_6 =  -w/2 
        fac_7=   -w/2 
        fac_8 =  -w/2
        fac_9 =   -w/2 
        fac_10 = -w/2
        fac_11 = -w/2
        fac_12 =  -w/2 
        fac_13 =  -w/2 
        fac_14 = -w/2
        fac_15 = -w/2
        fac_16 =  -w/2 
     elseif (imatrix == 4) then 
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = f/2
        fac_2 = 0
        fac_3 = f/2
        fac_4 = 0
        fac_5 = 0
        fac_6 = f/2
        fac_7=  0
        fac_8 = f/2
        fac_9 = 0
        fac_10 = f/2
        fac_11 = 0
        fac_12 = f/2
        fac_13 = f/2
        fac_14 = 0
        fac_15 = f/2
        fac_16 = 0
     else
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = -f/2
        fac_2 = -f/2 
        fac_3 = -f/2
        fac_4 = -f/2
        fac_5 = -f/2
        fac_6 = -f/2
        fac_7=  -f/2
        fac_8 = -f/2
        fac_9 = -f/2
        fac_10 =-f/2
        fac_11 = -f/2
        fac_12 = -f/2
        fac_13 = -f/2
        fac_14 = -f/2
        fac_15 = -f/2
        fac_16 = -f/2 
     end if
   end if
   if (keyword==7) then !singlet c=c
     if (imatrix==1) then
        fac_ed = w
        fac_ex = 0
        fac_hd = w
        fac_hx = w
        fac_1 = -w 
        fac_2 = 0
        fac_3 = 0
        fac_4 = -w
        fac_5 = -w
        fac_6 = 0
        fac_7= 0
        fac_8 = -w
        fac_9 = -w
        fac_10 = 0
        fac_11 = 0
        fac_12 = -w
        fac_13 = -w
        fac_14 = 0
        fac_15 = 0
        fac_16 = -w
     elseif (imatrix == 4) then 
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = 0.5*f
        fac_2 = 0
        fac_3 = 0
        fac_4 = 0.5*f
        fac_5 = 0.5*f
        fac_6 = 0 
        fac_7= 0
        fac_8 = 0.5*f
        fac_9 = 0.5*f
        fac_10 = 0
        fac_11 = 0
        fac_12 = 0.5*f
        fac_13 = 0.5*f
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0.5*f
     else
        fac_ed = f
        fac_ex = 0
        fac_hd = f
        fac_hx = f
        fac_1 = -f
        fac_2 = 0 
        fac_3 = 0
        fac_4 = -f
        fac_5 = -f
        fac_6 = 0
        fac_7= 0
        fac_8 = -f
        fac_9 = -f
        fac_10 = 0
        fac_11 = 0
        fac_12 = -f
        fac_13 = -f
        fac_14 = 0
        fac_15 = 0
        fac_16 = -f
    end if
  endif
  if (keyword==8) then !singlet14
    if (imatrix==1) then
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = 0
        fac_2 = 0 
        fac_3 = 0
        fac_4 = 0
        fac_5 = 0
        fac_6 = 0
        fac_7= 0
        fac_8 = 0
        fac_9 = 0
        fac_10 = 0
        fac_11 = 0
        fac_12 = 0
        fac_13 = 0
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0
     elseif (imatrix == 4) then 
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = -(r3t*f)/(r2t*2)
        fac_2 = -(r3t*f)/(r2t*2)
        fac_3 = (r3t*f)/(r2t*2)
        fac_4 = (r3t*f)/(r2t*2)
        fac_5 = (r3t*f)/(r2t*2)
        fac_6 = (r3t*f)/(r2t*2)
        fac_7= -(r3t*f)/(r2t*2)
        fac_8 =-(r3t*f)/(r2t*2)
        fac_9 = -(r3t*f)/(r2t*2)
        fac_10 =-(r3t*f)/(r2t*2)
        fac_11 = (r3t*f)/(r2t*2)
        fac_12 = (r3t*f)/(r2t*2)
        fac_13 = (r3t*f)/(r2t*2)
        fac_14 = (r3t*f)/(r2t*2)
        fac_15 = -(r3t*f)/(r2t*2)
        fac_16 = -(r3t*f)/(r2t*2)
     else
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = 0
        fac_2 = 0 
        fac_3 = 0
        fac_4 = 0
        fac_5 = 0
        fac_6 = 0
        fac_7= 0
        fac_8 = 0
        fac_9 = 0
        fac_10 = 0
        fac_11 = 0
        fac_12 = 0
        fac_13 = 0
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0

    end if
  endif
  if (keyword==9) then !singlet24
    if (imatrix==1) then
        fac_ed = w/(r2t)
        fac_ex = w/(r2t)
        fac_hd = 0
        fac_hx =0 
        fac_1 = -w/(r2t) 
        fac_2 = -w/(r2t)
        fac_3 = -w/(r2t)
        fac_4 = -w/(r2t)
        fac_5 = -w/(r2t)
        fac_6 = -w/(r2t)
        fac_7=  -w/(r2t)
        fac_8 = -w/(r2t)
        fac_9 = -w/(r2t)
        fac_10 = -w/(r2t)
        fac_11 = -w/(r2t)
        fac_12 = -w/(r2t)
        fac_13 = -w/(r2t)
        fac_14 = -w/(r2t)
        fac_15 = -w/(r2t)
        fac_16 = -w/(r2t)
     elseif (imatrix == 4) then 
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = f/(r2t*2)
        fac_2 = f/(r2t*2)
        fac_3 = f/(r2t*2)
        fac_4 = f/(r2t*2)
        fac_5 = f/(r2t*2)
        fac_6 = f/(r2t*2)
        fac_7=  f/(r2t*2)
        fac_8 = f/(r2t*2)
        fac_9 = f/(r2t*2)
        fac_10 = f/(r2t*2)
        fac_11 = f/(r2t*2)
        fac_12 = f/(r2t*2)
        fac_13 = f/(r2t*2)
        fac_14 = f/(r2t*2)
        fac_15 = f/(r2t*2)
        fac_16 = f/(r2t*2)
     else
        fac_ed = f/(r2t)
        fac_ex = f/(r2t)
        fac_hd = 0
        fac_hx = 0
        fac_1 = -f/(r2t)
        fac_2 = -f/(r2t)
        fac_3 = -f/(r2t)
        fac_4 = -f/(r2t)
        fac_5 = -f/(r2t)
        fac_6 = -f/(r2t)
        fac_7=  -f/(r2t)
        fac_8 = -f/(r2t)
        fac_9 = -f/(r2t)
        fac_10 = -f/(r2t)
        fac_11 = -f/(r2t)
        fac_12 = -f/(r2t)
        fac_13 = -f/(r2t)
        fac_14 = -f/(r2t)
        fac_15 = -f/(r2t)
        fac_16 = -f/(r2t)

    end if
  endif
  if (keyword==10) then !singlet34
    if (imatrix==1) then
        fac_ed = 0
        fac_ex = 0
        fac_hd = w/(r2t) 
        fac_hx = w/(r2t)  
        fac_1 = -w/(r2t) 
        fac_2 = 0
        fac_3 =  0
        fac_4 =  -w/(r2t)
        fac_5 =  -w/(r2t)
        fac_6 =  0
        fac_7=   0
        fac_8 =  -w/(r2t)
        fac_9 =  -w/(r2t)
        fac_10 = 0
        fac_11 = 0
        fac_12 =  -w/(r2t)
        fac_13 =  -w/(r2t)
        fac_14 =  0
        fac_15 =  0
        fac_16 =  -w/(r2t)
     elseif (imatrix == 4) then 
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = 0
        fac_2 = 0
        fac_3 = 0
        fac_4 = f/(r2t)
        fac_5 = f/(r2t)
        fac_6 = 0
        fac_7=  0
        fac_8 = 0
        fac_9 = 0
        fac_10 = 0
        fac_11 = 0
        fac_12 = f/(r2t)
        fac_13 = f/(r2t)
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0
     else
        fac_ed = 0
        fac_ex = 0
        fac_hd = f/(r2t)
        fac_hx = f/(r2t)
        fac_1 = -f/(r2t)
        fac_2 =0
        fac_3 =0
        fac_4 = -f/(r2t)
        fac_5 = -f/(r2t)
        fac_6 =0
        fac_7=  0
        fac_8 = -f/(r2t)
        fac_9 = -f/(r2t)
        fac_10 = 0
        fac_11 = 0
        fac_12 = -f/(r2t)
        fac_13 = -f/(r2t)
        fac_14 = 0
        fac_15 = 0
        fac_16 = -f/(r2t)

    end if
  endif
  if (keyword==11) then !singlet v=v
     if (imatrix==1) then
        fac_ed = w
        fac_ex = w
        fac_hd = w
        fac_hx = 0
        fac_1 = 0 
        fac_2 = 0
        fac_3 = 0
        fac_4 = 0
        fac_5 = 0
        fac_6 = 0
        fac_7= 0
        fac_8 = 0
        fac_9 = -w
        fac_10 = -w
        fac_11 = -w
        fac_12 = -w
        fac_13 = -w
        fac_14 = -w
        fac_15 = -w
        fac_16 = -w
     elseif (imatrix == 4) then 
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = 0
        fac_2 = 0
        fac_3 = 0
        fac_4 = 0
        fac_5 = 0
        fac_6 = 0 
        fac_7= 0
        fac_8 = 0
        fac_9 = 0.5*f
        fac_10 = 0.5*f
        fac_11 = 0.5*f
        fac_12 = 0.5*f
        fac_13 = 0.5*f
        fac_14 = 0.5*f
        fac_15 = 0.5*f
        fac_16 = 0.5*f
     else
        fac_ed = f
        fac_ex = f
        fac_hd = f
        fac_hx = 0
        fac_1 = 0
        fac_2 = 0 
        fac_3 = 0
        fac_4 = 0
        fac_5 = 0
        fac_6 = 0
        fac_7= 0
        fac_8 = 0
        fac_9 = -f
        fac_10 = -f
        fac_11 = -f
        fac_12 = -f
        fac_13 = -f
        fac_14 = -f
        fac_15 = -f
        fac_16 = -f
    end if
  endif
  if (keyword==12) then !H15
    if (imatrix==1) then
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = 0
        fac_2 = 0 
        fac_3 = 0
        fac_4 = 0
        fac_5 = 0
        fac_6 = 0
        fac_7= 0
        fac_8 = 0
        fac_9 = 0
        fac_10 = 0
        fac_11 = 0
        fac_12 = 0
        fac_13 = 0
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0
     elseif (imatrix == 4) then 
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = -(r3t*f)/(r2t*2)
        fac_2 = -(r3t*f)/(r2t*2)
        fac_3 = (r3t*f)/(r2t*2)
        fac_4 = (r3t*f)/(r2t*2)
        fac_5 = (r3t*f)/(r2t*2)
        fac_6 = (r3t*f)/(r2t*2)
        fac_7= -(r3t*f)/(r2t*2)
        fac_8 =-(r3t*f)/(r2t*2)
        fac_9 = -(r3t*f)/(r2t*2)
        fac_10 =-(r3t*f)/(r2t*2)
        fac_11 = (r3t*f)/(r2t*2)
        fac_12 = (r3t*f)/(r2t*2)
        fac_13 = (r3t*f)/(r2t*2)
        fac_14 = (r3t*f)/(r2t*2)
        fac_15 = -(r3t*f)/(r2t*2)
        fac_16 = -(r3t*f)/(r2t*2)
     else
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = 0
        fac_2 = 0 
        fac_3 = 0
        fac_4 = 0
        fac_5 = 0
        fac_6 = 0
        fac_7= 0
        fac_8 = 0
        fac_9 = 0
        fac_10 = 0
        fac_11 = 0
        fac_12 = 0
        fac_13 = 0
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0

    end if
  endif
  if (keyword==13) then !singlet25
    if (imatrix==1) then
        fac_ed = 0
        fac_ex = 0
        fac_hd = w/(r2t)
        fac_hx = w/(r2t) 
        fac_1 = -w/(r2t) 
        fac_2 = -w/(r2t)
        fac_3 = -w/(r2t)
        fac_4 = -w/(r2t)
        fac_5 = -w/(r2t)
        fac_6 = -w/(r2t)
        fac_7=  -w/(r2t)
        fac_8 = -w/(r2t)
        fac_9 =  -w/(r2t)
        fac_10 = -w/(r2t)
        fac_11 = -w/(r2t)
        fac_12 = -w/(r2t)
        fac_13 = - w/(r2t)
        fac_14 = -w/(r2t)
        fac_15 =  -w/(r2t)
        fac_16 =  -w/(r2t)
     elseif (imatrix == 4) then 
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = f/(r2t*2)
        fac_2 = f/(r2t*2)
        fac_3 = f/(r2t*2)
        fac_4 = f/(r2t*2)
        fac_5 = f/(r2t*2)
        fac_6 = f/(r2t*2)
        fac_7=  f/(r2t*2)
        fac_8 = f/(r2t*2)
        fac_9 = f/(r2t*2)
        fac_10 = f/(r2t*2)
        fac_11 = f/(r2t*2)
        fac_12 = f/(r2t*2)
        fac_13 = f/(r2t*2)
        fac_14 = f/(r2t*2)
        fac_15 = f/(r2t*2)
        fac_16 = f/(r2t*2)
     else
        fac_ed = 0
        fac_ex = 0
        fac_hd = f/(r2t)
        fac_hx = f/(r2t)
        fac_1 = -f/(r2t)
        fac_2 = -f/(r2t)
        fac_3 = -f/(r2t)
        fac_4 = -f/(r2t)
        fac_5 = -f/(r2t)
        fac_6 = -f/(r2t)
        fac_7=  -f/(r2t)
        fac_8 = -f/(r2t)
        fac_9 = -f/(r2t)
        fac_10 = -f/(r2t)
        fac_11 = -f/(r2t)
        fac_12 = -f/(r2t)
        fac_13 = -f/(r2t)
        fac_14 = -f/(r2t)
        fac_15 = -f/(r2t)
        fac_16 = -f/(r2t)

    end if
  endif
  if (keyword==14) then !singlet35
    if (imatrix==1) then
        fac_ed =w/(r2t) 
        fac_ex =w/(r2t)
        fac_hd = 0
        fac_hx = 0  
        fac_1 = 0 
        fac_2 = 0 
        fac_3 =  0 
        fac_4 =  0 
        fac_5 =  0 
        fac_6 =  0 
        fac_7=   0 
        fac_8 =  0 
        fac_9 = -w/(r2t)
        fac_10 = -w/(r2t)
        fac_11 = -w/(r2t)
        fac_12 = -w/(r2t)
        fac_13 = -w/(r2t)
        fac_14 = -w/(r2t)
        fac_15 = -w/(r2t)
        fac_16 = -w/(r2t)
     elseif (imatrix == 4) then 
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = 0
        fac_2 = 0
        fac_3 = 0
        fac_4 = 0
        fac_5 = 0
        fac_6 = 0
        fac_7=  0
        fac_8 = 0
        fac_9 = 0
        fac_10 = 0
        fac_11 = f/(r2t)
        fac_12 =  f/(r2t)
        fac_13 =  f/(r2t)
        fac_14 =  f/(r2t)
        fac_15 = 0
        fac_16 = 0
     else
        fac_ed = f/(r2t)
        fac_ex = f/(r2t)
        fac_hd = 0
        fac_hx = 0
        fac_1 =0
        fac_2 =0
        fac_3 =0
        fac_4 =0
        fac_5 =0
        fac_6 =0
        fac_7= 0
        fac_8 =0
        fac_9 = -f/(r2t)
        fac_10 = -f/(r2t)
        fac_11 = -f/(r2t)
        fac_12 = -f/(r2t)
        fac_13 = -f/(r2t)
        fac_14 = -f/(r2t)
        fac_15 = -f/(r2t)
        fac_16 = -f/(r2t)

    end if
  endif
  if (keyword==15) then !singlet45
    if (imatrix==1) then
        fac_ed = w/2
        fac_ex = w/2
        fac_hd = w/2
        fac_hx = w/2  
        fac_1 = -w/2
        fac_2 = -w/2
        fac_3 =  -w/2
        fac_4 =  -w/2
        fac_5 =  -w/2
        fac_6 =  -w/2
        fac_7=   -w/2
        fac_8 =  -w/2
        fac_9 =  -w/2
        fac_10 = -w/2
        fac_11 = -w/2
        fac_12 =  -w/2
        fac_13 =  -w/2
        fac_14 =  -w/2
        fac_15 =  -w/2
        fac_16 =  -w/2
     elseif (imatrix == 4) then 
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = f/2
        fac_2 = f/2
        fac_3 = 0
        fac_4 = 0
        fac_5 = 0
        fac_6 = 0
        fac_7=  f/2
        fac_8 = f/2
        fac_9 = 0
        fac_10 = 0
        fac_11 = f/2
        fac_12 = f/2
        fac_13 = f/2
        fac_14 = f/2
        fac_15 = 0
        fac_16 = 0
     else
        fac_ed = f/2
        fac_ex = f/2
        fac_hd = f/2
        fac_hx = f/2
        fac_1 = -f/2
        fac_2 =-f/2
        fac_3 =-f/2
        fac_4 = -f/2
        fac_5 = -f/2
        fac_6 =-f/2
        fac_7=  -f/2
        fac_8 = -f/2
        fac_9 = -f/2
        fac_10 = -f/2
        fac_11 = -f/2
        fac_12 = -f/2
        fac_13 = -f/2
        fac_14 = -f/2
        fac_15 = -f/2
        fac_16 = -f/2

    end if
  endif
  if (keyword==16) then !singlet
    if (imatrix==1) then
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = -w
        fac_2 = 0
        fac_3 = 0
        fac_4 = 0
        fac_5 = 0
        fac_6 = 0
        fac_7= 0
        fac_8 = 0
        fac_9 = 0
        fac_10 = 0
        fac_11 = 0
        fac_12 = 0
        fac_13 = 0
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0
    elseif (imatrix == 4) then
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = 2*f
        fac_2 = 0
        fac_3 = 0
        fac_4 = 0
        fac_5 = 0
        fac_6 = 0
        fac_7=  0
        fac_8 = 0
        fac_9 = 0
        fac_10 = 0
        fac_11 = 0
        fac_12 = 0
        fac_13 = 0
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0
    else
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = -f
        fac_2 = 0
        fac_3 = 0
        fac_4 = 0
        fac_5 = 0
        fac_6 = 0
        fac_7= 0
        fac_8 = 0
        fac_9 = 0
        fac_10 = 0
        fac_11 = 0
        fac_12 = 0
        fac_13 = 0
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0
    end if
  endif
  if (keyword==17) then !singlet16
    if (imatrix==1) then
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0 
        fac_hx = 0  
        fac_1 = (r3t/r2t)*w
        fac_2 =-(r3t/r2t)*w
        fac_3 = -(r3t/r2t)*w
        fac_4 =(r3t/r2t)*w
        fac_5 =-(r3t/r2t)*w
        fac_6 = (r3t/r2t)*w
        fac_7= (r3t/r2t)*w
        fac_8 = -(r3t/r2t)*w
        fac_9 =  0
        fac_10 = 0
        fac_11 = 0
        fac_12 = 0
        fac_13 = 0
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0
     elseif (imatrix == 4) then 
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = 0
        fac_2 = 0
        fac_3 = 0
        fac_4 = 0
        fac_5 = 0
        fac_6 = 0
        fac_7=  0
        fac_8 = 0
        fac_9 = 0
        fac_10 = 0
        fac_11 = 0
        fac_12 = 0
        fac_13 = 0
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0
     else
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0 
        fac_hx = 0  
        fac_1 = (r3t/r2t)*f
        fac_2 =-(r3t/r2t)*f
        fac_3 = -(r3t/r2t)*f
        fac_4 =(r3t/r2t)*f
        fac_5 =-(r3t/r2t)*f
        fac_6 = (r3t/r2t)*f
        fac_7= (r3t/r2t)*f
        fac_8 = -(r3t/r2t)*f
        fac_9 =  0
        fac_10 = 0
        fac_11 = 0
        fac_12 = 0
        fac_13 = 0
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0
        print*,"root3/2",r3t/r2t
    end if
  endif
  if (keyword==18) then !singlet26
    if (imatrix==1) then
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0 
        fac_hx = 0  
        fac_1 = w/r2t
        fac_2 = w/r2t
        fac_3 = -w/r2t
        fac_4 = -w/r2t
        fac_5 = w/r2t
        fac_6 = w/r2t
        fac_7=  -w/r2t
        fac_8 = -w/r2t
        fac_9 =  0
        fac_10 = 0
        fac_11 = 0
        fac_12 = 0
        fac_13 = 0
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0
     elseif (imatrix == 4) then 
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = 0 
        fac_2 = 0
        fac_3 = 0
        fac_4 = 0
        fac_5 = 0
        fac_6 = 0
        fac_7=  0
        fac_8 = 0
        fac_9 = 0
        fac_10 = 0
        fac_11 = 0
        fac_12 = 0
        fac_13 = 0
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0
     else
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0 
        fac_hx = 0  
        fac_1 = f/r2t
        fac_2 = f/r2t
        fac_3 = -f/r2t
        fac_4 = -f/r2t
        fac_5 = f/r2t
        fac_6 = f/r2t
        fac_7=  -f/r2t
        fac_8 = -f/r2t
        fac_9 =  0
        fac_10 = 0
        fac_11 = 0
        fac_12 = 0
        fac_13 = 0
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0
     end if
  endif
  if (keyword==19) then !singlet36
    if (imatrix==1) then
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0 
        fac_hx = 0  
        fac_1 = w/r2t
        fac_2 = 0
        fac_3 = -w/r2t
        fac_4 = 0
        fac_5 = 0
        fac_6 = w/r2t
        fac_7=  0
        fac_8 = -w/r2t
        fac_9 =  0
        fac_10 = 0
        fac_11 = 0
        fac_12 = 0
        fac_13 = 0
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0
     elseif (imatrix == 4) then 
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = 0 
        fac_2 = 0
        fac_3 = 0
        fac_4 = 0
        fac_5 = 0
        fac_6 = 0
        fac_7=  0
        fac_8 = 0
        fac_9 = 0
        fac_10 = 0
        fac_11 = 0
        fac_12 = 0
        fac_13 = 0
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0
     else
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0 
        fac_hx = 0  
        fac_1 = f/r2t
        fac_2 = 0
        fac_3 = -f/r2t
        fac_4 = 0
        fac_5 = 0
        fac_6 = f/r2t
        fac_7= 0
        fac_8 = -f/r2t
        fac_9 =  0
        fac_10 = 0
        fac_11 = 0
        fac_12 = 0
        fac_13 = 0
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0
     end if
  endif
  if (keyword==20) then !singlet46
    if (imatrix==1) then
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0 
        fac_hx = 0  
        fac_1 = w/2
        fac_2 = w/2
        fac_3 = -w/2
        fac_4 = -w/2
        fac_5 = w/2
        fac_6 = w/2
        fac_7=  -w/2
        fac_8 = -w/2
        fac_9 =  0
        fac_10 = 0
        fac_11 = 0
        fac_12 = 0
        fac_13 = 0
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0
     elseif (imatrix == 4) then 
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = 0
        fac_2 = 0
        fac_3 = 0
        fac_4 = 0
        fac_5 = 0
        fac_6 = 0
        fac_7=  0
        fac_8 = 0
        fac_9 = 0
        fac_10 = 0
        fac_11 = 0
        fac_12 = 0
        fac_13 = 0
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0
     else
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0 
        fac_hx = 0  
        fac_1 = f/2
        fac_2 = f/2
        fac_3 = -f/2
        fac_4 = -f/2
        fac_5 = f/2
        fac_6 = f/2
        fac_7=  -f/2
        fac_8 = -f/2
        fac_9 =  0
        fac_10 = 0
        fac_11 = 0
        fac_12 = 0
        fac_13 = 0
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0
     end if
  endif
  if (keyword==21) then !singlet56
    if (imatrix==1) then
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0 
        fac_hx = 0  
        fac_1 = w/2
        fac_2 = w/2
        fac_3 = -w/2
        fac_4 = -w/2
        fac_5 = w/2
        fac_6 = w/2
        fac_7=  -w/2
        fac_8 = -w/2
        fac_9 =  0
        fac_10 = 0
        fac_11 = 0
        fac_12 = 0
        fac_13 = 0
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0
     elseif (imatrix == 4) then 
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = 0
        fac_2 = 0
        fac_3 = 0
        fac_4 = 0
        fac_5 = 0
        fac_6 = 0
        fac_7=  0
        fac_8 = 0
        fac_9 = 0
        fac_10 = 0
        fac_11 = 0
        fac_12 = 0
        fac_13 = 0
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0
     else
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0 
        fac_hx = 0  
        fac_1 = f/2
        fac_2 = f/2
        fac_3 = -f/2
        fac_4 = -f/2
        fac_5 = f/2
        fac_6 = f/2
        fac_7=  -f/2
        fac_8 = -f/2
        fac_9 =  0
        fac_10 = 0
        fac_11 = 0
        fac_12 = 0
        fac_13 = 0
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0
     end if
  endif
  if (keyword==23) then !singlet
    if (imatrix==1) then
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = r3t*w
        fac_2 = -r3t*w
        fac_3 = 0
        fac_4 = 0
        fac_5 = 0
        fac_6 = 0
        fac_7= 0
        fac_8 = 0
        fac_9 = 0
        fac_10 = 0
        fac_11 = 0
        fac_12 = 0
        fac_13 = 0
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0
    elseif (imatrix == 4) then
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = 0
        fac_2 = 0
        fac_3 = 0
        fac_4 = 0
        fac_5 = 0
        fac_6 = 0
        fac_7=  0
        fac_8 = 0
        fac_9 = 0
        fac_10 = 0
        fac_11 = 0
        fac_12 = 0
        fac_13 = 0
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0
    else
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = r3t*f
        fac_2 = -r3t*f
        fac_3 = 0
        fac_4 = 0
        fac_5 = 0
        fac_6 = 0
        fac_7= 0
        fac_8 = 0
        fac_9 = 0
        fac_10 = 0
        fac_11 = 0
        fac_12 = 0
        fac_13 = 0
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0
    end if
  endif
  if (keyword==24) then !singlet
    if (imatrix==1) then
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = w
        fac_2 = w
        fac_3 = 0
        fac_4 = 0
        fac_5 = 0
        fac_6 = 0
        fac_7= 0
        fac_8 = 0
        fac_9 = 0
        fac_10 = 0
        fac_11 = 0
        fac_12 = 0
        fac_13 = 0
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0
    elseif (imatrix == 4) then
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = 0
        fac_2 = 0
        fac_3 = 0
        fac_4 = 0
        fac_5 = 0
        fac_6 = 0
        fac_7=  0
        fac_8 = 0
        fac_9 = 0
        fac_10 = 0
        fac_11 = 0
        fac_12 = 0
        fac_13 = 0
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0
    else
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = f
        fac_2 = f
        fac_3 = 0
        fac_4 = 0
        fac_5 = 0
        fac_6 = 0
        fac_7= 0
        fac_8 = 0
        fac_9 = 0
        fac_10 = 0
        fac_11 = 0
        fac_12 = 0
        fac_13 = 0
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0
    end if
  endif
  if (keyword==25) then !singlet
    if (imatrix==1) then
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = w
        fac_2 = 0
        fac_3 = 0
        fac_4 = 0
        fac_5 = 0
        fac_6 = 0
        fac_7= 0
        fac_8 = 0
        fac_9 = 0
        fac_10 = 0
        fac_11 = 0
        fac_12 = 0
        fac_13 = 0
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0
    elseif (imatrix == 4) then
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = 0
        fac_2 = 0
        fac_3 = 0
        fac_4 = 0
        fac_5 = 0
        fac_6 = 0
        fac_7=  0
        fac_8 = 0
        fac_9 = 0
        fac_10 = 0
        fac_11 = 0
        fac_12 = 0
        fac_13 = 0
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0
    else
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = f
        fac_2 = 0
        fac_3 = 0
        fac_4 = 0
        fac_5 = 0
        fac_6 = 0
        fac_7= 0
        fac_8 = 0
        fac_9 = 0
        fac_10 = 0
        fac_11 = 0
        fac_12 = 0
        fac_13 = 0
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0
    end if
  endif
  if (keyword==26) then !singlet
    if (imatrix==1) then
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = w/r2t
        fac_2 = w/r2t
        fac_3 = 0
        fac_4 = 0
        fac_5 = 0
        fac_6 = 0
        fac_7= 0
        fac_8 = 0
        fac_9 = 0
        fac_10 = 0
        fac_11 = 0
        fac_12 = 0
        fac_13 = 0
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0
    elseif (imatrix == 4) then
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = 0
        fac_2 = 0
        fac_3 = 0
        fac_4 = 0
        fac_5 = 0
        fac_6 = 0
        fac_7=  0
        fac_8 = 0
        fac_9 = 0
        fac_10 = 0
        fac_11 = 0
        fac_12 = 0
        fac_13 = 0
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0
    else
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = f/r2t
        fac_2 = f/r2t
        fac_3 = 0
        fac_4 = 0
        fac_5 = 0
        fac_6 = 0
        fac_7= 0
        fac_8 = 0
        fac_9 = 0
        fac_10 = 0
        fac_11 = 0
        fac_12 = 0
        fac_13 = 0
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0
    end if
  endif
  if (keyword==27) then !singlet
    if (imatrix==1) then
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = w/r2t
        fac_2 = w/r2t
        fac_3 = 0
        fac_4 = 0
        fac_5 = 0
        fac_6 = 0
        fac_7= 0
        fac_8 = 0
        fac_9 = 0
        fac_10 = 0
        fac_11 = 0
        fac_12 = 0
        fac_13 = 0
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0
    elseif (imatrix == 4) then
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = 0
        fac_2 = 0
        fac_3 = 0
        fac_4 = 0
        fac_5 = 0
        fac_6 = 0
        fac_7=  0
        fac_8 = 0
        fac_9 = 0
        fac_10 = 0
        fac_11 = 0
        fac_12 = 0
        fac_13 = 0
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0
    else
        fac_ed = 0
        fac_ex = 0
        fac_hd = 0
        fac_hx = 0
        fac_1 = f/r2t
        fac_2 = f/r2t
        fac_3 = 0
        fac_4 = 0
        fac_5 = 0
        fac_6 = 0
        fac_7= 0
        fac_8 = 0
        fac_9 = 0
        fac_10 = 0
        fac_11 = 0
        fac_12 = 0
        fac_13 = 0
        fac_14 = 0
        fac_15 = 0
        fac_16 = 0
    end if
  endif
end subroutine fac_values

end module hdf5_read_m
