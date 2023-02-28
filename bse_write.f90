 program hdf5_write
     
    use hdf5_read_m
    use create_H
    use hdf5
    implicit none

    double precision :: wcoul0,vol,ellow,elup,abstol,t1,t2
    
    complex(kind=8), dimension(:), allocatable :: mf_v
    complex(kind=8), dimension(:), allocatable :: mf_c

    integer :: nb,info,ilow,iup,error,nv,nc,r,s,t,x,y,m,q,n,o,p,l,ml,z,c

    
     
     !!!!! read the QP-energeis and calculate the size of matrix !!!!!!!!!
     call read_inp(nv,nc,vol)
     allocate(mf_v(nv))
     allocate(mf_c(nc))
     call load_mf(nv,nc,mf_v,mf_c)
     call load_weff(wcoul0)
     nb = nv + nc
     s= nc*(nc-1)*nv*(nv-1)/4     !size of 1 + 2
     p = nc*(nc-1)*nv/2
     q = nv*(nv-1)*nc/2
     n = (2*s) + (nv*nc)        !size of 1+2+3
     o = (2*s) + (nv*nc) + q    !size of 1+2+3+4
     m = (2*s) + (nv*nc) + q + p !size of 1+2+3+4+5
     l = (2*s) + (2*(nv*nc)) + q + p !size of 1+2+3+4+5+6
     ml = (nv*nc) + q + p !size of 3+4+5
     z = l+1
     print*,l
     c = nc !variable that stores chunk size
     
     call MPI_INIT(error)
     call cpu_time(t1)
     !! write the hamiltonian within group1 elements!!!!!!!!!
     call create_H1(nv,nc,mf_v,mf_c,wcoul0,s,vol,c)
     print*,"2"
      call create_H2(nv,nc,mf_v,mf_c,wcoul0,s,vol,c)
      print*,"2"

     call create_H12(nv,nc,mf_v,mf_c,wcoul0,s,vol,c)
      print*,"2"

     !call create_H3(nv,nc,mf_v,mf_c,wcoul0,s,vol)
      print*,"2"

     call create_H13(nv,nc,mf_v,mf_c,wcoul0,s,vol,c)
      print*,"2"

     call create_H23(nv,nc,mf_v,mf_c,wcoul0,s,vol,c)
      print*,"2"

     !call create_H4(nv,nc,mf_v,mf_c,wcoul0,s,vol,q,n)
      print*,"H4"

     call create_H14(nv,nc,mf_v,mf_c,wcoul0,s,vol,q,n,c)
      print*,"H14"

     call create_H24(nv,nc,mf_v,mf_c,wcoul0,s,vol,q,n,c)
      print*,"H24"

     !call create_H34(nv,nc,mf_v,mf_c,wcoul0,s,vol,q,n)
      print*,"H34"

     call create_H5(nv,nc,mf_v,mf_c,wcoul0,s,vol,p,o,c)

      print*,"2"

     call create_H15(nv,nc,mf_v,mf_c,wcoul0,s,vol,p,o,c)
      print*,"2"

     call create_H25(nv,nc,mf_v,mf_c,wcoul0,s,vol,p,o,c)
      print*,"2"

     !call create_H35(nv,nc,mf_v,mf_c,wcoul0,s,vol,p,o)
      print*,"2"

     !call create_H45(nv,nc,mf_v,mf_c,wcoul0,s,vol,q,p,n,o)
      print*,"2"

     !call create_H6(nv,nc,mf_v,mf_c,wcoul0,s,vol,m)
     call create_H16(nv,nc,mf_v,mf_c,wcoul0,s,vol,m,c)
     call create_H26(nv,nc,mf_v,mf_c,wcoul0,s,vol,m,c)
    ! call create_H36(nv,nc,mf_v,mf_c,wcoul0,s,vol,m)
     !call create_H46(nv,nc,mf_v,mf_c,wcoul0,s,vol,q,n,m)
     
     call create_H56(nv,nc,mf_v,mf_c,wcoul0,s,vol,p,o,m,c)

     !call cpu_time(t2)
     print*,l,t1-t2

    deallocate(mf_v)
    deallocate(mf_c)
    call MPI_FINALIZE(error)



 end program hdf5_write
