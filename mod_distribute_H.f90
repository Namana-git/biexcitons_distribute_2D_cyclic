module distribute_H


      public :: calculate_chunk_size

      contains

      subroutine calculate_chunk_size(chunk_size)

               integer , intent(inout)  :: chunk_size 
               integer         :: tot_length,c,r 


               tot_length =sys_var%nc*(sys_var%nc-1)/2
               
               c = tot_length/grid%nprow
               r =  mod(tot_length,grid%nprow)

               if (grid%myprow < r) then
                   
                  chunk_size = c + 1 
               else
                  chunk_size = c
               
               endif

        



      end subroutine calculate_chunk_size






































       



end module distribute_H
