MODULE m_linear_algebra
   USE m_common
   implicit none
    !> 2D Vector type
    !!
    !! \details This type is used for the definition of vector of size 2 ( Coordinates \f$(x,y)\f$ for example).
   TYPE vec2d
      real(rp) :: x
      real(rp) :: y
   END TYPE vec2d
    !> 3D Vector type
    !!
    !! \details This type is used for the definition of vector of size 3 ( Coordinates \f$(x,y,z)\f$ for example).
   TYPE vec3d
      real(rp) :: x , y , z
   END TYPE vec3d
    !> 2D Tensor Structure
    !!
    !! \details This type is used for the definition of Tensor of size 2*2 ( ... for example).
   TYPE tens2d
      real(rp) :: xx , xy , yx , yy
   END TYPE tens2d
   type(tens2d) :: Id2d
    !> 3D Tensor Structure
    !!
    !! \details This type is used for the definition of Tensor of size 3*3 ( ... for example).
   TYPE tens3d
      real(rp) :: xx , xy , xz , yx , yy , yz , zx , zy , zz
   END TYPE tens3d
   type(tens3d) :: Id3d
    !> Linear System Structure ( Full Matrix )
   TYPE sys_lin_full
        !> size ???
      integer(ip) :: n
        !> Data array ??? (COLUMN MAJOR : column index, row index)
      real(rp), dimension(:,:), allocatable :: a
        !> ???
      real(rp), dimension(:), allocatable :: rhs
        !> determinant ???
      real(rp) :: det
   END TYPE sys_lin_full
    !> Linear System Structure ( Sparse Matrix )
   TYPE sys_lin_sparse
      integer(ip) :: n , nz
      integer(ip), dimension(:), allocatable :: ia , ja , swap
      real(rp), dimension(:), allocatable :: a , x , x0 , rhs
   END TYPE sys_lin_sparse
CONTAINS
END MODULE m_linear_algebra
