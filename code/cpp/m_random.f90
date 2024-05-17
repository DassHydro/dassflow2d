MODULE m_random
   USE m_common
   implicit none
CONTAINS
   SUBROUTINE init_random_seed
      implicit none
      integer(ip) :: clock !> random value generation (based on time ?)
      integer(ip) :: iseed !> id of the vector to modify
      integer(ip) :: nseed !> size of the random vector
      integer(ip), dimension(:), allocatable :: seed
      call random_seed( size = nseed )
      allocate( seed( nseed ) )
      call system_clock( count = clock )
      seed = clock + 37 * (/ ( iseed - 1 , iseed = 1 , nseed ) /)
      call random_seed( put = seed )
      deallocate( seed )
      return
   END SUBROUTINE init_random_seed
END MODULE m_random
