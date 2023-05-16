! Module for the whole program
! Aim: declaring the  variables that are common between the subroutines
! Author: Leonard Santos, Irstea-HBAN, France (leonard.santos@irstea.fr) 
! Date: 2017/09/28

   MODULE m_gr4
      
   USE m_common
   
   implicit none

 ! Modelling time-step information
      character :: Tstep
     
  ! Tstep="H" for hourly time-step
  ! Tstep="D" for daily time-step
 
           
 ! Number of data values at the different time-steps (total number and warm-up)
      integer(ip)   ::  Ndata_Hour ! Number of hourly time step to compute
      integer(ip)   ::  Nwup_Hour ! Number of hourly warmup time step to compute    !Keep for DF??



 ! Number of stores of the Nash cascade
      !integer(ip)   ::  nres  = 11
      integer(ip)   ::  catchnb
 

!  ! Parameters and state vectors of the model
    !  real(rp), dimension(13)  ::  State
      
     
  ! Param(1): x1 ! Max capacity of the production store [mm]
  ! Param(2): x2 ! Inter-catchment exchange coefficient [mm/t]
  ! Param(3): x3 ! Max capacity of the routing store [mm]
  ! Param(4): x4 ! Base time of the unit hydrograph [t]
  
  ! State(1):          S ! Production store level [mm]
  ! State(2:nres+1): Sh ! Nash cascade stores levels [mm]
  ! State(nres+2):    R ! Routing store level [mm]
 

 ! Input data
      integer, dimension(:), allocatable :: Year,Month,Day,Hour
      real(rp), dimension(:), allocatable :: Precip,Pot_Ev
      real(rp), dimension(:), allocatable :: Obs_Fl

  ! Year(:): year at each time-step
  ! Month(:): month number at each time-step
  ! Day(:): day number at each time-step
  ! Hour(:): hour at each time-step
  
  ! Precip(:): Rainfall amount at each time-step [mm]
  ! Pot_Ev(:): Potential evapotranspiration at each time-step [mm]
  
  ! Obs_Fl(:): Observed flow at each time-step [mm]
 

 ! Output data
      real(rp),  dimension(19) ::  MISC
      real(rp),  dimension(:,:), allocatable :: MISCO
     
  ! MISC( 1): PE     ! potential evapotranspiration [mm/t]
  ! MISC( 2): Precip ! total precipitation [mm/t]
  ! MISC( 3): Prod   ! production store level [mm]
  ! MISC( 4): RatPro ! production store ratio level [-]
  ! MISC( 5): ETR    ! actual evapotranspiration [mm/t]
  ! MISC( 6): Perc   ! percolation from production store [mm/t]
  ! MISC( 7): PS     ! rain in production store [mm/t]
  ! MISC( 8): PR     ! precipitation inflow of the Nash cascade [mm/t]
  ! MISC( 9): Q9     ! inflow in the main routing branche [mm/t]
  ! MISC(10): Q1     ! inflow in the secondary routing branche [mm/t]
  ! MISC(11): Rout   ! routing store level [mm]
  ! MISC(12): RatRou ! routing store level ratio [-]
  ! MISC(13): Exch   ! potential semi-exchange between catchments [mm/t]
  ! MISC(14): AExch  ! actual total exchange between catchments [mm/t]
  ! MISC(15): QR     ! outflow from routing store [mm/t]
  ! MISC(16): QD     ! outflow from secondary branch after exchange [mm/t]
  ! MISC(17): Qsim   ! outflow at catchment outlet [mm/t]
  ! MISC(18)
  ! MISC(19): Casc   ! Nash cascade stores levels [mm]
  
  ! MISCO(:,k): MISC variable number k at each time-step
 

 ! Catchment informations
      real(rp) ::  Utils(2)

  ! Utils(1): Catchment area [km2]
  ! Utils(2): Men observed flow [mm]
 

 ! Performances array
      real(rp) ::  Perf(10)
      
 

      end module m_gr4
