! Model launch subroutine
! Aim: Temporal loop to launch the model at each time-step and get output data
! Author: Leonard Santos, Irstea-HBAN, France (leonard.santos@irstea.fr)
! Date: 2017/09/28

      Subroutine MOD_LAUNCH
      
      USE m_common
      USE m_model
      USE m_gr4
      
      Implicit none
      
 ! Variables declaration
      integer(ip) ::   ll
      integer(ip) ::   ll_tot
      
      real(rp) ::   Rain_gr4
      real(rp) ::   Evap
 
       character(50)  ::  filename
       
       
 ! Choice of the total number of time-steps
   ! and allocate output array
   Ndata_Hour = SIZE(bc%gr4(catchnb)%t)

!       if (Tstep.eq."H") then
        ll_tot=Ndata_Hour
        allocate(MISCO(Ndata_Hour,19))
!       else
!         ll_tot=Ndata_Day
!         allocate(MISCO(Ndata_Day,19))
!       endif
 

 ! Time loop
 write(filename,'("res/GR4outputs_",1I4.4,".txt")') catchnb

 
      open(10,file=filename,status='replace',form='formatted')
      write(10,*) ' Hour E(mm/h) P(mm/h) Q(mm/h) States(13, in mm)'
      do ll=1,ll_tot
        ! Extract input data
        Rain_gr4=bc%gr4(catchnb)%P(ll)
        Evap=bc%gr4(catchnb)%E(ll)

        ! Launch the model
        Call GR4_STSP(Evap,Rain_gr4,ll)
        write(10,*) ll,Evap, Rain_gr4, MISC(17), bc%gr4(catchnb)%state
        ! Get output data
        MISCO(ll,:)=MISC(:)


        
      enddo
      close(10)

      end subroutine MOD_LAUNCH
      
      
      
      Subroutine MOD_WARMUP
      
      USE m_common
      USE m_model
      USE m_gr4
      
      Implicit none
      
 ! Variables declaration
      integer(ip) ::   ll
      integer(ip) ::   ll_tot
      
      real(rp) ::   Rain_gr4
      real(rp) ::   Evap
 
       character(50)  ::  filename
       
       
 ! Choice of the total number of time-steps
   ! and allocate output array
   Ndata_Hour = size(bc%gr4(catchnb)%P0)

        ll_tot=Ndata_Hour
        allocate(MISCO(Ndata_Hour,19))


      ! Time loop
      do ll=1,ll_tot
        ! Extract input data
        Rain_gr4=bc%gr4(catchnb)%P0(ll)
        Evap=bc%gr4(catchnb)%E0(ll)

        ! Launch the model
        Call GR4_STSP(Evap,Rain_gr4,ll)
        
        ! Get output data
        ! MISCO(ll,:)=MISC(:)


      enddo

      end subroutine MOD_WARMUP
