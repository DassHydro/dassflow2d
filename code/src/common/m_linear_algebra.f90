!======================================================================================================================!
!
!                    DassFlow Version 3.0
!
!======================================================================================================================!
!
!  Copyright University of Toulouse-INSA, Univ. of Strasbourg, INRAE & CNRS (France)
!
!  This file is part of the DassFlow software (Data Assimilation for Free Surface Flows).
!  DassFlow is a computational software whose purpose is to simulate geophysical free surface flows,
!  designed for variational sensitivities and data assimilation (4D-var). Inverse capabilities are
!  based on the adjoint code generation by a source-to-source algorithmic differentiation (Tapenade software used).
!
!  DassFlow software includes few mostly independent "modules" with common architectures and structures.
!  Please consult the DassFlow webpage for more details: http://www-gmm.insa-toulouse.fr/~monnier/DassFlow/.
!
!  Many people have contributed to the DassFlow development from the initial version to the latest ones.
!  Current contributions:
!               L. Pujol (PhD Unistra)
!               L. Villenave (PhD student)
!               P.-A. Garambois (INRAE Aix-en-Provence)
!               J. Monnier (INSA & Mathematics Institute of Toulouse IMT).
!               K. Larnier (CS group - IMT-INSA).
!  Former scientific or programming contributions of:
!               F. Couderc (CNRS & Mathematics Institute of Toulouse IMT).
!               J.-P. Vila (INSA & Mathematics Institute of Toulouse IMT).
!               R. Madec   (Mathematics Institute of Toulouse IMT).
!  plus less recent other developers (M. Honnorat and J. Marin).
!
!  Contact : see the DassFlow webpage 
!
!  This software is governed by the CeCILL license under French law and abiding by the rules of distribution
!  of free software. You can use, modify and/or redistribute the software under the terms of the CeCILL license
!  as circulated by CEA, CNRS and INRIA at the following URL: "http://www.cecill.info".
!
!  As a counterpart to the access to the source code and rights to copy, modify and redistribute granted by the
!  license, users are provided only with a limited warranty and the software's author, the holder of the economic
!  rights, and the successive licensors have only limited liability.
!
!  In this respect, the user's attention is drawn to the risks associated with loading, using, modifying and/or
!  developing or reproducing the software by the user in light of its specific status of free software, that may
!  mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and
!  experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the
!  software's suitability as regards their requirements in conditions enabling the security of their systems and/or
!  data to be ensured and, more generally, to use and operate it in the same conditions as regards security.
!
!  The fact that you are presently reading this means that you have had knowledge of the CeCILL license and that you
!  accept its terms.
!
!======================================================================================================================!
!> @file m_linear_algebra.f90
!> @brief This file includes the m_linear_algebra module.



!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Module m_linear_algebra.f90
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

!> @brief Module m_linear_algebra.
!> @details TO ADD : CHOLESKY, etc... for var_chg
MODULE m_linear_algebra

   USE m_common

   implicit none

!======================================================================================================================!
!  2D Vector Structure
!======================================================================================================================!
    !> 2D Vector type
    !!
    !! \details This type is used for the definition of vector of size 2 ( Coordinates \f$(x,y)\f$ for example).
   TYPE vec2d

      real(rp)  ::  x 
      real(rp) :: y

   END TYPE vec2d

!======================================================================================================================!
!  3D Vector Structure
!======================================================================================================================!
    !> 3D Vector type
    !!
    !! \details This type is used for the definition of vector of size 3 ( Coordinates \f$(x,y,z)\f$ for example).
   TYPE vec3d

      real(rp)  ::  x , y , z

   END TYPE vec3d

!======================================================================================================================!
!  2D Tensor Structure
!======================================================================================================================!
    !> 2D Tensor Structure
    !!
    !! \details This type is used for the definition of Tensor of size 2*2  ( ... for example).
   TYPE tens2d

      real(rp)  ::  xx , xy , yx , yy

   END TYPE tens2d

   type(tens2d)  :: Id2d

!======================================================================================================================!
!  3D Tensor Structure
!======================================================================================================================!
    !> 3D Tensor Structure
    !!
    !! \details This type is used for the definition of Tensor of size 3*3  ( ... for example).
   TYPE tens3d

      real(rp)  ::  xx , xy , xz , yx , yy , yz , zx , zy , zz

   END TYPE tens3d

   type(tens3d)  :: Id3d

!======================================================================================================================!
!  Linear System Structure ( Full Matrix )
!======================================================================================================================!
    !> Linear System Structure ( Full Matrix )
   TYPE sys_lin_full
        !> size ??? 
      integer(ip)  ::  n
        !> Data array  ??? (COLUMN MAJOR : column index, row index)
      real(rp), dimension(:,:), allocatable  ::  a
        !>   ??? 
      real(rp), dimension(:), allocatable  ::  rhs
        !>   determinant ???
      real(rp)  ::  det

   END TYPE sys_lin_full

!======================================================================================================================!
!  Linear System Structure ( Sparse Matrix )
!======================================================================================================================!
    !>  Linear System Structure ( Sparse Matrix )
   TYPE sys_lin_sparse

      integer(ip)  ::  n , nz

      integer(ip), dimension(:), allocatable  ::  ia , ja , swap

      real(rp), dimension(:), allocatable  ::  a , x , x0 , rhs

   END TYPE sys_lin_sparse

!======================================================================================================================!
!  Operator overloading interfaces
!======================================================================================================================!

#ifndef CPP_PASS
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
#endif

CONTAINS
#ifndef CPP_PASS

!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Linear Algebra Basic Operations
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


!======================================================================================================================!
!  Operator (+)
!======================================================================================================================!

    !>   sum of 2D vectors
   FUNCTION add_vec2d( vec1 , vec2 ) RESULT( res )

      implicit none

      type(vec2d), intent(in)  ::  vec1 , vec2

      type(vec2d)  ::  res

      res%x  =  vec1%x  +  vec2%x
      res%y  =  vec1%y  +  vec2%y

   END FUNCTION add_vec2d

    !>   sum of 2D tensors
   FUNCTION add_tens2d( tens1 , tens2 ) RESULT( res )

      implicit none

      type(tens2d), intent(in)  ::  tens1 , tens2

      type(tens2d)  ::  res

      res%xx  =  tens1%xx  +  tens2%xx
      res%xy  =  tens1%xy  +  tens2%xy
      res%yx  =  tens1%yx  +  tens2%yx
      res%yy  =  tens1%yy  +  tens2%yy

   END FUNCTION add_tens2d


!======================================================================================================================!
!  Operator (-)
!======================================================================================================================!

    !>   substraction of 2D vectors
   FUNCTION sub_vec2d( vec1 , vec2 ) RESULT( res )

      implicit none

      type(vec2d), intent(in)  ::  vec1 , vec2

      type(vec2d)  ::  res

      res%x  =  vec1%x  -  vec2%x
      res%y  =  vec1%y  -  vec2%y

   END FUNCTION sub_vec2d

    !>   substraction of 2D tensors
   FUNCTION sub_tens2d( tens1 , tens2 ) RESULT( res )

      implicit none

      type(tens2d), intent(in)  ::  tens1 , tens2

      type(tens2d)  ::  res

      res%xx  =  tens1%xx  -  tens2%xx
      res%xy  =  tens1%xy  -  tens2%xy
      res%yx  =  tens1%yx  -  tens2%yx
      res%yy  =  tens1%yy  -  tens2%yy

   END FUNCTION sub_tens2d


!======================================================================================================================!
!  Operator (*)
!======================================================================================================================!

    !>    multiplication of 2D vectors with a scalar (r)
   FUNCTION scalar_vec2d( r , vec ) RESULT( res )

      implicit none

      type(vec2d), intent(in)  ::  vec

      real(rp), intent(in)  ::  r

      type(vec2d)  ::  res

      res%x  =  r * vec%x
      res%y  =  r * vec%y

   END FUNCTION scalar_vec2d

    !>   multiplication of 2D vectors with a scalar (r)  == scalar_vec2d ?
   FUNCTION vec2d_scalar( vec , r ) RESULT( res )

      implicit none

      type(vec2d), intent(in)  ::  vec

      real(rp), intent(in)  ::  r

      type(vec2d)  ::  res

      res%x  =  r * vec%x
      res%y  =  r * vec%y

   END FUNCTION vec2d_scalar

    !>   multiplication of 2D tensors with a scalar (r)
   FUNCTION scalar_tens2d( r , vec ) RESULT( res )

      implicit none

      type(tens2d), intent(in)  ::  vec

      real(rp), intent(in)  ::  r

      type(tens2d)  ::  res

      res%xx  =  r * vec%xx
      res%xy  =  r * vec%xy
      res%yx  =  r * vec%yx
      res%yy  =  r * vec%yy

   END FUNCTION scalar_tens2d

    !>   multiplication of 2D tensors with a scalar (r) == scalar_tens2d?
   FUNCTION tens2d_scalar( vec , r ) RESULT( res )

      implicit none

      type(tens2d), intent(in)  ::  vec

      real(rp), intent(in)  ::  r

      type(tens2d)  ::  res

      res%xx  =  r * vec%xx
      res%xy  =  r * vec%xy
      res%yx  =  r * vec%yx
      res%yy  =  r * vec%yy

   END FUNCTION tens2d_scalar

    !>   multiply a tensor  with two vectors 
   FUNCTION tens2d_vec2d( tens , vec ) RESULT( res )

      implicit none

      type(tens2d), intent(in)  ::  tens

      type(vec2d), intent(in)  ::  vec

      type(vec2d)  ::  res

      res%x  =  tens%xx * vec%x + tens%xy * vec%y
      res%y  =  tens%yx * vec%x + tens%yy * vec%y

   END FUNCTION tens2d_vec2d


!======================================================================================================================!
!  Operator (/)
!======================================================================================================================!

    !>   division of 2D vectors with a scalar (r) 
   FUNCTION vec2d_div_scalar( vec , r ) RESULT( res )

      implicit none

      type(vec2d), intent(in)  ::  vec

      real(rp), intent(in)  ::  r

      type(vec2d)  ::  res

      res%x  =  vec%x / r
      res%y  =  vec%y / r

   END FUNCTION vec2d_div_scalar

    !>   division of 2D tensor with a scalar (r) 
   FUNCTION tens2d_div_scalar( tens , r ) RESULT( res )

      implicit none

      type(tens2d), intent(in)  ::  tens

      real(rp), intent(in)  ::  r

      type(tens2d)  ::  res

      res%xx  =  tens%xx / r
      res%xy  =  tens%xy / r
      res%yx  =  tens%yx / r
      res%yy  =  tens%yy / r

   END FUNCTION tens2d_div_scalar


!======================================================================================================================!
!  Operator (.norm.)
!======================================================================================================================!

    !>   calculate (euclidian) norm of a 2D vector
   FUNCTION norm_vec2d( vec ) RESULT( res )

      implicit none

      type(vec2d), intent(in)  ::  vec

      real(rp)  ::  res

      res  =  sqrt( vec%x**2 + vec%y**2 )

   END FUNCTION norm_vec2d


!======================================================================================================================!
!  Operator (.dotprod.)
!======================================================================================================================!

    !>   dot_product of two 2Dvectors
   FUNCTION dot_product_vec2d( vec1 , vec2 ) RESULT( res )

      implicit none

      type(vec2d), intent(in)  ::  vec1 , vec2

      real(rp)  ::  res

      res  =  vec1%x * vec2%x + &
              vec1%y * vec2%y

   END FUNCTION dot_product_vec2d


!======================================================================================================================!
!  Operator (.tensprod.)
!======================================================================================================================!
    !>   tensorial  product of two 2Dvectors
   FUNCTION tens_prod2d( vec1 , vec2 ) RESULT( res )

      implicit none

      type(vec2d), intent(in)  ::  vec1 , vec2

      type(tens2d)  ::  res

      res%xx  =  vec1%x * vec2%x
      res%xy  =  vec1%x * vec2%y

      res%yx  =  vec1%y * vec2%x
      res%yy  =  vec1%y * vec2%y

   END FUNCTION tens_prod2d


!======================================================================================================================!
!  Operator (.transpose.)
!======================================================================================================================!

    !>   transpose a tensor 2D
   FUNCTION transpose_tens2d( tens ) RESULT( res )

      implicit none

      type(tens2d), intent(in)  ::  tens

      type(tens2d)  ::  res

      res%xx  =  tens%xx
      res%xy  =  tens%yx
      res%yx  =  tens%xy
      res%yy  =  tens%yy

   END FUNCTION transpose_tens2d
#endif


END MODULE m_linear_algebra
