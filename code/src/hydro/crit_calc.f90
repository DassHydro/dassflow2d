! Perfomance calculation subroutine
! Aim: calculating the performances of the model using different flow transformations and different criterions
! Author: Leonard Santos, Irstea-HBAN, France (leonard.santos@irstea.fr)
! Date: 2017/09/28

      SUBROUTINE CRIT_CALC(SIM,OBS,Trans,Critnam,Per)
      
      USE m_common
      USE m_model
      USE m_gr4
      
      Implicit none
      
! Variables declaration

      Character*2 Trans
      Character*3 Critnam
      
      integer(ip)  ::  ll
      
      integer(ip)  ::  ll_tot
      integer(ip)  ::  ll_wup
      integer(ip)  ::  IO1
      
      real(rp) ::   QQ,CC,eps
      real(rp) ::   SQ1,SQ2,SC1,SC2,SCQ,SDIFF
      real(rp) ::   MEAQ,MEAC,SIGQ,SIGC,CORQC,COVQC
      
      real(rp) ::   SIM(*),OBS(*)
      
      real(rp) ::   Per
      
      

 ! Initialization of criterion components

      IO1=0      ! Number of observations
      SQ1=0.     ! Sum of observed flows (transformed)
      SQ2=0.     ! Sum of squared observed flows (transformed)
      SC1=0.     ! Sum of simulated flows (transformed)
      SC2=0.     ! Sum of squared simulated flows (transformed)
      SCQ=0.     ! Sum of products between observed and simulated flows (transformed)
      SDIFF=0.   ! Sum of errors between observed and simulated flows (transformed)

 

 ! Choice of the total and warm-up number of time step
      if (Tstep.eq."H") then
      
        ll_tot=SIZE(bc%gr4(catchnb)%t)
        ll_wup=Nwup_Hour
      
      else
      
        ll_tot=0!Ndata_Day
        ll_wup=0!Nwup_Day
      
      endif
 
      
 ! Time loop
      do ll=(ll_wup+1),ll_tot
        
        QQ=OBS(ll)  ! observed flow at the time step
        CC=SIM(ll)  ! simulated flow at the time step
        
        if ( QQ.ge.0.) then  ! avoid the gaps in observed measurment
        
   ! Transformation of the flows (Pushpalatha et al., 2012)
          ! "HQ" ! High flows: no transformation
          ! "LQ" ! Low flows: logarithmic transformation
          ! "IQ" ! Low flows: inverted transformation
          ! "QI" ! Transitional flows: square root transformation
          ! "PQ" ! Peak flows: square transformation
           
          ! epsilon to avoid null flows in log and inverted transformantion
          ! calculated as a function of mean flow and catchment area (Pushpalatha et al., 2012)
          eps=0.01*(Utils(2)*Utils(1)/86.4)
          
          if (Trans.eq."LQ") then
            QQ=log(QQ+eps)
            CC=log(CC+eps)
          elseif (Trans.eq."QI") then
            QQ=sqrt(QQ)
            CC=sqrt(CC)
          elseif (Trans.eq."PQ") then
            QQ=QQ*QQ
            CC=CC*CC
          elseif (Trans.eq."IQ") then
            QQ=1./(QQ+eps)
            CC=1./(CC+eps)
          endif

 

   ! Components calculation

          IO1=IO1+1
          
          SQ1=SQ1+QQ
          SQ2=SQ2+QQ*QQ
        
          SC1=SC1+CC
          SC2=SC2+CC*CC
          
          SCQ=SCQ+QQ*CC
          SDIFF=SDIFF+(QQ-CC)*(QQ-CC)

 
        
        endif
        
      enddo
 

 ! Statistic values calculated using components
      ! Means
      MEAQ=SQ1/float(IO1)
      MEAC=SC1/float(IO1)
      ! Standard deviations
      SIGQ=sqrt(SQ2/float(IO1)-MEAQ*MEAQ)
      SIGC=sqrt(SC2/float(IO1)-MEAC*MEAC)
      ! Covariance
      COVQC=SCQ/float(IO1)-MEAQ*MEAC
      ! Pearson corelation coefficient
      CORQC=COVQC/(SIGQ*SIGC)
 

 ! Performance calculation

      ! KGE'
      if (Critnam.eq.'KGE') then

        Per=1-sqrt((CORQC-1)**2+((SIGC*MEAQ)/(SIGQ*MEAC)-1)**2+ &
        (MEAC/MEAQ-1)**2)
     
      ! NSE
      elseif (Critnam.eq.'NSE') then
      
        Per=1-SDIFF/(SQ2-SQ1*SQ1/float(IO1))
      
      endif

 

      end SUBROUTINE CRIT_CALC


