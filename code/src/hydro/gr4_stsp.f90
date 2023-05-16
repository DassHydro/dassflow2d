! Model core subroutine
! Aim: Running the model on a single time-step by solving the differential equations
! Author: Leonard Santos, Irstea-HBAN, France (leonard.santos@irstea.fr)
! Date: 2017/09/28

      SUBROUTINE GR4_STSP(E,P,ll)
      
      USE m_common
      USE m_model
      USE m_gr4
      
      Implicit none
      
 ! Variables declaration      
      ! Loop index
      integer(ip), intent(in) ::   ll
      integer(ip)  ::  ll_ststp
   
      ! Reservoir count
      integer(ip) ::   nres = 11
      
      ! Fixed parameters of the model
      integer(ip) ::   al
      integer(ip) ::   be
      integer(ip) ::   ga
      real(rp) ::   om
      real(rp) ::   nu
      real(rp) ::   Phi
      real(rp), dimension(4) ::   Coef
      real(rp), dimension(4) ::   Coef_st

      
      ! Inputs
      real(rp) ::   E
      real(rp) ::   EN
      real(rp) ::   P
      real(rp) ::   PN
      real(rp), dimension(4) ::   Inpu
      real(rp), dimension(4) ::   Inpu_st
      
      ! Time-step and sub-time-step
      real(rp) ::   tstp
      real(rp) ::   ststp
      integer(ip) ::   nststp
      
      ! Temporary state vector and values
      real(rp), dimension(4)  ::  Param

      real(rp), dimension(13) ::   State_st
      real(rp), dimension(13) ::   State_st1
      real(rp), dimension(13) ::   State_st_rec
      real(rp) ::   Temp
      real(rp) ::   S0
      real(rp) ::   S1
      real(rp) ::   fS1
      real(rp) ::   fS0
      
      ! Substepping parameters
      real(rp) ::   StopVal
      real(rp) ::   Err
      real(rp) ::   Tol
      
      ! Internal fluxes variables
      real(rp) ::   ETR
      real(rp) ::   QD
      real(rp) ::   QR
      real(rp) ::   PERC
      real(rp) ::   PS
      real(rp) ::   PR
      real(rp) ::   Q9
      real(rp) ::   Q1
      real(rp) ::   ECH
      real(rp) ::   ECT
      
      ! Output of the model
      real(rp) ::    Q     
 

 
 
 ! Time step value
      !if (Tstep.eq."H") then
        tstp=0.0417
!       else
!         tstp=1.0000
!       endif
 

 ! Fixed parameters of the model
      al=2         ! alpha :: Net rainfall exponent
      be=5         ! beta  :: Percolation exponent
      ga=5         ! gamma :: Routing outflow exponent
      om=7./2.     ! omega :: Exchange exponent
      nu=4./9.     ! nu    :: Percolation coefficient
      Phi=0.9      ! Phi   :: Flow partition coefficient
 

 ! Calculation of outflow coefficient (taking time-step into account)
   ! Stored in the same array
   ! See Tab. 1 in Santos et al. (2017) for full equations definitions
      Coef(1)=tstp*(nu/bc%gr4(catchnb)%params(1))**(be-1)/(be-1)  ! kperc :: percolation coefficient (Perc = kperc * S**beta)
      Coef(2)=tstp*(nres-1)/bc%gr4(catchnb)%params(4)             ! kcasc :: Nash cascade coefficient (Qsh,k = kcasc * Sh,k)
      Coef(3)=tstp/(bc%gr4(catchnb)%params(3)**(ga-1)*(ga-1))     ! krout :: routing store outflow coefficient (QR = krout * R  gamma)
      Coef(4)=tstp*bc%gr4(catchnb)%params(2)/bc%gr4(catchnb)%params(3)**om         ! kech  :: routing store exchange coefficient (F = kech * R  omega)
 
 ! Interception step (neutralization between rainfall and evap)
      ! Create input data for the state-space model
      E=E*(1./tstp)
      P=P*(1./tstp)
      EN=max(E-P,0.)
      PN=max(P-E,0.)
 

 ! Calculation of input coefficient values after interception (taking time-step into account)
   ! Stored in the same array
      Inpu(1)=tstp*PN                 ! Net rainfall after interception step
      Inpu(2)=tstp*PN/bc%gr4(catchnb)%params(1)**al    ! Quadratic input in the Nash cascade coefficient
      Inpu(3)=tstp*EN/bc%gr4(catchnb)%params(1)**al    ! Positive component of the evap from production store
      Inpu(4)=tstp*2.*EN/bc%gr4(catchnb)%params(1)     ! Negative component of the evap from production store
 
 ! Initialization and fixed values for adaptaive time-stepping
      StopVal=0.01_rp        ! Stop value for the secant method
      Tol=0.1_rp         ! Fixed error tolerance at 0.1
      Err=12._rp               ! Error initialization (greater than tol)
      State_st_rec=bc%gr4(catchnb)%state    ! Initialize the reference to calculate the substepping error (Err)
      nststp=1_ip              ! Initialize the value and the number of substep
      ststp=1._rp/float(nststp)       
 

 ! Adaptative subtime-stepping (Press et al., 1992)
      do while (Err.gt.Tol)
       
        State_st=bc%gr4(catchnb)%state   ! Initialization of the state array for subtime step
         
   ! Subtime coefficients
        Coef_st=ststp*Coef
        Inpu_st=ststp*Inpu
 

   ! Initialization of internal fluxes values calculation (see MISC at the end for the explations)
        ETR=min(P,E)*tstp
        QD=0._rp
        QR=0._rp
        PERC=0._rp
        PS=0._rp
        PR=0._rp
        Q9=0._rp
        Q1=0._rp
        ECH=0._rp
        ECT=0._rp
 

   ! Loop on all the sub-time-steps     

        do ll_ststp=1,nststp
         
   ! Production store
     ! Secant method to resolve the Euler implicit equation on the store level

       ! Initial solution guess of the equation
          S0=State_st(1)      ! Former state value
          S1=S0-Coef_st(1)*S0**be+Inpu_st(1)-Inpu_st(2)*S0**al+&
          Inpu_st(3)*S0**al-Inpu_st(4)*S0  ! substep explicit Euler for a good start point
       
       ! Initial result of the equation
          call fzer_prod(S0,State_st(1),Coef_st,Inpu_st,fS0)
          call fzer_prod(S1,State_st(1),Coef_st,Inpu_st,fS1)
           
     ! Secant iterations
          do while (abs(fS1).gt.StopVal)
            Temp=S1-fS1*(S1-S0)/(fS1-fS0)  ! Sequant iteration calculation
            
            S0=S1                          ! Keep the second limit
            fS0=fS1
             
            S1=Temp                        ! Calculate the solution for the new limit
            call fzer_prod(S1,State_st(1),Coef_st,Inpu_st,fS1)
          enddo
 
          State_st1(1)=S1                  ! record the resulting production store level
 
 

   ! Nash cascade store 1
       ! analytical resolution of the implict Euler equation
          State_st1(2)=(State_st(2)+Inpu_st(2)*State_st1(1)**al+&
          Coef_st(1)*State_st1(1)**be)/(1+Coef_st(2))
 

    ! All other Nash cascade stores
       ! analytical resolution of the implict Euler equation
          do k=3,nres+1
            State_st1(k)=(State_st(k)+Coef_st(2)*State_st1(k-1))/&
            (1+Coef_st(2))
          enddo
 

    ! Routing store
     ! Secant method to solve the implicit Euler equation on the store level

       ! Initial solution guess of the equation
          S0=State_st(nres+2)      ! Former state value
          S1=S0+Phi*Coef_st(2)*State_st1(nres+1)+Coef_st(4)*S0**om-&
          Coef_st(3)*S0**ga                 ! substep explicit Euler for a good start point
       
       ! Initial result of the equation
          call fzer_rout(S0,State_st(nres+2),State_st1(nres+1),Coef_st,&
          Phi,fS0)
          call fzer_rout(S1,State_st(nres+2),State_st1(nres+1),Coef_st,&
          Phi,fS1)
           
     ! Secant iterations
          do while (abs(fS1).gt.StopVal)
            Temp=S1-fS1*(S1-S0)/(fS1-fS0)  ! Sequant iteration calculation
             
            S0=S1                          ! Keep the second limit
            fS0=fS1
             
            S1=Temp                        ! Calculate the solution for the new limit
            call fzer_rout(S1,State_st(nres+2),State_st1(nres+1),Coef_st&
            ,Phi,fS1)
          enddo
 
          State_st1(nres+2)=S1             ! record the resulting routing store level
 
 

     ! Internal fluxes calculation by summing all substep fluxes values
          QD=QD+(1-Phi)*Coef_st(2)*State_st1(nres+1)+Coef_st(4)*State_st1(nres+2)**om
          QR=QR+Coef_st(3)*State_st1(nres+2)**ga

          ETR=ETR+Inpu_st(4)*State_st1(1)-Inpu_st(3)*State_st(1)**al
          PERC=PERC+Coef_st(1)*State_st1(1)**be
          PS=PS+Inpu_st(1)-Inpu_st(2)*State_st1(1)**al
          PR=PR+Inpu_st(2)*State_st1(1)**al+Coef_st(1)*State_st1(1)**be
          Q9=Q9+Phi*Coef_st(2)*State_st1(nres+1)
          Q1=Q1+(1-Phi)*Coef_st(2)*State_st1(nres+1)
          ECH=ECH+Coef_st(4)*State_st1(nres+2)**om
          if (Coef_st(4).gt.0.) then
            ECT=ECT+2*Coef_st(4)*State_st1(nres+2)**om
          else
            ECT=ECT+Coef_st(4)*State_st1(nres+2)**om-min(abs((1-Phi)*&
            Coef_st(2)*State_st1(1+nres)),abs(Coef_st(4)*State_st1(nres+2)**om))
          endif
 
  
     ! Update of the state variables for the following subtime-step
          State_st=State_st1
 
           
        enddo    ! End of loop on all the substeps

    ! Procedure for updating the new subtime-step value      
          
         if (nststp.gt.100) exit ! finish the loop if the substep is over 100
          
         ! Calculation of error as the maximum difference between the new and the previous states
         Err=maxval(abs(State_st-State_st_rec))  
         State_st_rec=State_st  ! Recording of the new state calculation
          
         ! Update the subtime step value (Press et al., 1992)
         ststp=0.9_rp*ststp*(Tol/Err)
         nststp=int(1./ststp)  ! get the number of substeps
         ststp=1./float(nststp)  ! get an exact correspondance between the number and the value of substeps
          
         ! Maximum number of substeps : 150
         if (nststp.gt.150_ip) then
           nststp=150_ip
          ststp=1._rp/float(nststp)
         endif

 
         
      enddo    ! End of while loop for adaptative step
 

 ! Flow calcultaion at the outlet for the time-step
      Q=QR+max(0.,QD)

 ! Update of the state variables for the following time-step
      bc%gr4(catchnb)%state=State_st
 

 ! Array of output variables
      MISC( 1)=E*tstp                  ! PE     ! potential evapotranspiration  [mm/t]
      MISC( 2)=P*tstp                  ! Precip ! total precipitation  [mm/t]
      MISC( 3)=bc%gr4(catchnb)%state(1)                ! Prod   ! production store level [mm]
      MISC( 4)=bc%gr4(catchnb)%state(1)/bc%gr4(catchnb)%params(1)       ! RatPro ! production store ratio level [-]
      MISC( 5)=ETR                     ! ETR    ! actual evapotranspiration [mm/t]
      MISC( 6)=PERC                    ! Perc   ! percolation from production store [mm/t]
      MISC( 7)=PS                      ! PS     ! rain in production store [mm/t]
      MISC( 8)=PR                      ! PR     ! precipitation inflow of the Nash cascade [mm/t]
      MISC( 9)=Q9                      ! Q9     ! inflow in the main routing branche [mm/t]
      MISC(10)=Q1                      ! Q1     ! inflow in the secondary routing branche [mm/t]
      MISC(11)=bc%gr4(catchnb)%state(2+nres)           ! Rout   ! routing store level [mm]
      MISC(12)=bc%gr4(catchnb)%state(2+nres)/bc%gr4(catchnb)%params(3) ! RatRou ! routing store level ratio [-]
      MISC(13)=ECH                     ! Exch   ! potential semi-exchange between catchments [mm/t]
      MISC(14)=ECT                     ! AExch  ! actual total exchange between catchments [mm/t]
      MISC(15)=QR                      ! QR     ! outflow from routing store [mm/t]
      MISC(16)=max(0.,QD)              ! QD     ! outflow from secondary branch after exchange [mm/t]
      MISC(17)=Q
      bc%gr4(catchnb)%Q(ll)=Q * bc%gr4(catchnb)%surf / 3.6_rp                       ! Qsim   ! outflow at catchment outlet [mm/h] => [m3/s]
      MISC(18)=bc%gr4(catchnb)%state(2)
      MISC(19)=bc%gr4(catchnb)%state(1+nres)           ! Casc   ! Nash cascade stores levels [mm]
      
      !bc%gr4(catchnb)%state(:) = State(:)
 
      end SUBROUTINE GR4_STSP


! Euler Implicit equation for production store
! Aim: Implement the equation to be solved by sequant method for the production store
! Author: Leonard Santos, Irstea-HBAN, France (leonard.santos@irstea.fr)

      SUBROUTINE fzer_prod(x,y,Coef,Inp,out)

      USE m_common
      
      implicit none

      real(rp), dimension(4) ::   Coef
      real(rp), dimension(4) ::   Inp
      real(rp) ::   y,x
      real(rp) ::   out

 ! Function to bring to zero (x is the variable, y the previous state)
      out=y-Coef(1)*x*x*x*x*x+Inp(1)-Inp(2)*x*x+Inp(3)*x*x-Inp(4)*x-x
      
      end SUBROUTINE fzer_prod
      

! Euler Implicit equation for routing store
! Aim: Implement the equation to be solved by sequant method for the routing store
! Author: Leonard Santos, Irstea-HBAN, France (leonard.santos@irstea.fr)

      SUBROUTINE fzer_rout(x,y,Casc,Coef,Phi,out)
      
      USE m_common

      implicit none

      real(rp), dimension(4) ::   Coef
      real(rp) ::   Phi
      real(rp) ::   Casc
      real(rp) ::   y,x
      real(rp) ::   out
 ! Function to bring to zero (x is the variable, y the previous state)
      out=y+Phi*Coef(2)*Casc-Coef(3)*x*x*x*x*x+Coef(4)*x*x*x*sqrt(x)-x
      
      end SUBROUTINE fzer_rout
