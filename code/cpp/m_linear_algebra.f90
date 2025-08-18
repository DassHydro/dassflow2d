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
   interface operator (+)
      module procedure add_vec2d , add_tens2d
   end interface
   interface operator (-)
      module procedure sub_vec2d , sub_tens2d
   end interface
   interface operator (*)
      module procedure scalar_vec2d , vec2d_scalar , scalar_tens2d , tens2d_scalar , tens2d_vec2d
   end interface
   interface operator (/)
      module procedure vec2d_div_scalar , tens2d_div_scalar
   end interface
   interface operator (.tensprod.)
      module procedure tens_prod2d
   end interface
   interface operator (.norm.)
      module procedure norm_vec2d
   end interface
   interface operator (.dotprod.)
      module procedure dot_product_vec2d
   end interface
   interface operator (.t.)
      module procedure transpose_tens2d
   end interface
CONTAINS
    !> sum of 2D vectors
   FUNCTION add_vec2d( vec1 , vec2 ) RESULT( res )
      implicit none
      type(vec2d), intent(in) :: vec1 , vec2
      type(vec2d) :: res
      res%x = vec1%x + vec2%x
      res%y = vec1%y + vec2%y
   END FUNCTION add_vec2d
    !> sum of 2D tensors
   FUNCTION add_tens2d( tens1 , tens2 ) RESULT( res )
      implicit none
      type(tens2d), intent(in) :: tens1 , tens2
      type(tens2d) :: res
      res%xx = tens1%xx + tens2%xx
      res%xy = tens1%xy + tens2%xy
      res%yx = tens1%yx + tens2%yx
      res%yy = tens1%yy + tens2%yy
   END FUNCTION add_tens2d
    !> substraction of 2D vectors
   FUNCTION sub_vec2d( vec1 , vec2 ) RESULT( res )
      implicit none
      type(vec2d), intent(in) :: vec1 , vec2
      type(vec2d) :: res
      res%x = vec1%x - vec2%x
      res%y = vec1%y - vec2%y
   END FUNCTION sub_vec2d
    !> substraction of 2D tensors
   FUNCTION sub_tens2d( tens1 , tens2 ) RESULT( res )
      implicit none
      type(tens2d), intent(in) :: tens1 , tens2
      type(tens2d) :: res
      res%xx = tens1%xx - tens2%xx
      res%xy = tens1%xy - tens2%xy
      res%yx = tens1%yx - tens2%yx
      res%yy = tens1%yy - tens2%yy
   END FUNCTION sub_tens2d
    !> multiplication of 2D vectors with a scalar (r)
   FUNCTION scalar_vec2d( r , vec ) RESULT( res )
      implicit none
      type(vec2d), intent(in) :: vec
      real(rp), intent(in) :: r
      type(vec2d) :: res
      res%x = r * vec%x
      res%y = r * vec%y
   END FUNCTION scalar_vec2d
    !> multiplication of 2D vectors with a scalar (r) == scalar_vec2d ?
   FUNCTION vec2d_scalar( vec , r ) RESULT( res )
      implicit none
      type(vec2d), intent(in) :: vec
      real(rp), intent(in) :: r
      type(vec2d) :: res
      res%x = r * vec%x
      res%y = r * vec%y
   END FUNCTION vec2d_scalar
    !> multiplication of 2D tensors with a scalar (r)
   FUNCTION scalar_tens2d( r , vec ) RESULT( res )
      implicit none
      type(tens2d), intent(in) :: vec
      real(rp), intent(in) :: r
      type(tens2d) :: res
      res%xx = r * vec%xx
      res%xy = r * vec%xy
      res%yx = r * vec%yx
      res%yy = r * vec%yy
   END FUNCTION scalar_tens2d
    !> multiplication of 2D tensors with a scalar (r) == scalar_tens2d?
   FUNCTION tens2d_scalar( vec , r ) RESULT( res )
      implicit none
      type(tens2d), intent(in) :: vec
      real(rp), intent(in) :: r
      type(tens2d) :: res
      res%xx = r * vec%xx
      res%xy = r * vec%xy
      res%yx = r * vec%yx
      res%yy = r * vec%yy
   END FUNCTION tens2d_scalar
    !> multiply a tensor with two vectors
   FUNCTION tens2d_vec2d( tens , vec ) RESULT( res )
      implicit none
      type(tens2d), intent(in) :: tens
      type(vec2d), intent(in) :: vec
      type(vec2d) :: res
      res%x = tens%xx * vec%x + tens%xy * vec%y
      res%y = tens%yx * vec%x + tens%yy * vec%y
   END FUNCTION tens2d_vec2d
    !> division of 2D vectors with a scalar (r)
   FUNCTION vec2d_div_scalar( vec , r ) RESULT( res )
      implicit none
      type(vec2d), intent(in) :: vec
      real(rp), intent(in) :: r
      type(vec2d) :: res
      res%x = vec%x / r
      res%y = vec%y / r
   END FUNCTION vec2d_div_scalar
    !> division of 2D tensor with a scalar (r)
   FUNCTION tens2d_div_scalar( tens , r ) RESULT( res )
      implicit none
      type(tens2d), intent(in) :: tens
      real(rp), intent(in) :: r
      type(tens2d) :: res
      res%xx = tens%xx / r
      res%xy = tens%xy / r
      res%yx = tens%yx / r
      res%yy = tens%yy / r
   END FUNCTION tens2d_div_scalar
    !> calculate (euclidian) norm of a 2D vector
   FUNCTION norm_vec2d( vec ) RESULT( res )
      implicit none
      type(vec2d), intent(in) :: vec
      real(rp) :: res
      res = sqrt( vec%x**2 + vec%y**2 )
   END FUNCTION norm_vec2d
    !> dot_product of two 2Dvectors
   FUNCTION dot_product_vec2d( vec1 , vec2 ) RESULT( res )
      implicit none
      type(vec2d), intent(in) :: vec1 , vec2
      real(rp) :: res
      res = vec1%x * vec2%x + &
              vec1%y * vec2%y
   END FUNCTION dot_product_vec2d
    !> tensorial product of two 2Dvectors
   FUNCTION tens_prod2d( vec1 , vec2 ) RESULT( res )
      implicit none
      type(vec2d), intent(in) :: vec1 , vec2
      type(tens2d) :: res
      res%xx = vec1%x * vec2%x
      res%xy = vec1%x * vec2%y
      res%yx = vec1%y * vec2%x
      res%yy = vec1%y * vec2%y
   END FUNCTION tens_prod2d
    !> transpose a tensor 2D
   FUNCTION transpose_tens2d( tens ) RESULT( res )
      implicit none
      type(tens2d), intent(in) :: tens
      type(tens2d) :: res
      res%xx = tens%xx
      res%xy = tens%yx
      res%yx = tens%xy
      res%yy = tens%yy
   END FUNCTION transpose_tens2d
END MODULE m_linear_algebra
