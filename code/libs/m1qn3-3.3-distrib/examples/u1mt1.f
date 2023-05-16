c-----------------------------------------------------------------------
c
c     Test problem: u1mt1 (meteorology), version 1 in double precision
c
c     The difference with version 0 is in the computation of the
c     gradient g of the function f, which is done in the subroutine
c     "u1mt1t" (only this subroutine has been changed). 
c     Here the gradient is computed by an "automatique differentiation"
c     technique, the reverse mode of differentiation implemented in an
c     adjoint code (see the review paper "La differentiation automatique
c     de fonctions representees par des programmes" by J.Ch. Gilbert,
c     G. Le Vey and J. Masse (1991) INRIA Research Report).
c     This technique is about twice faster than the one in the original
c     code.
c
c-----------------------------------------------------------------------
c
      implicit none
      integer i, j, ilec, irep, isav, n2, n3, nv, np, itlag, itmax,
     &  nemax, izs(1)
      real rzs(1)
      double precision cola, sp, epsu, epsv, epsp, au, av, ap, alx, aly,
     &  alpha, beta, gamma, udx, vdx, pdx, x(1875), ff, g(1875), epsf,
     &  zero, ec, rr, dxz(1875), dzs(1875)
c    &  e1u(625), e1v(625), e1p(625)
ctmp
      double precision ytab, stab
      common /ystab/ ytab(1875, 10), stab(1875, 10)
ctmp
c
c --- common
c
      integer n
      double precision dx, dy, fm, fc, u0, v0, ph0, x0, eu, ev, ep, ru,
     &  rv, rp
      common /u1mt1c/ dx, dy, fm(625), fc(625), u0(625), v0(625),
     &  ph0(625), x0(1875), eu(625), ev(625), ep(625), ru(625), rv(625),
     &  rp(625), n
c
c --- external
c
      external u1mt1s, euclid, ctonbe, ctcabe
c
c --- working zone for the optimizer
c
      character*3 normtype
      integer impres, io, iter, nsim, imode(3), omode, ntrav, iz(5),
     &  nmtrav, indic, reverse
      parameter (nmtrav=85000)
      double precision dxmin, df1, epsrel, trav(nmtrav)
c     integer nfich
c     parameter (nfich=9)
c     character fich*80
c
c --- fichiers
c
c     open (6,form='print')
      ilec=14
      irep=15
      isav=16
c
c --- initialisation
c
      n=25
      n2=n*n
      nv=n2+1
      np=2*n2+1
      n3=3*n2
      cola=0.1d+0
      epsu=1.d-5
      epsv=1.d-5
      epsp=1.d-5
      alx=2500000.d+0
      aly=2500000.d+0
      dx=alx/dble(n-1)
      dy=aly/dble(n-1)
c     gaa=10.d+0
c     bmult=2.d+0
c
c --- lecture du fichier de donnees
c
      open (ilec,file="u1mt1.data",
     /      status='old',access='sequential',form='formatted')
      call u1mt1l (ilec)
      read (ilec,*) alpha,beta,gamma
      read (ilec,*) udx,vdx,pdx
c
c     point initial et multiplicateurs
c
      do i=1,n3
          x(i)=x0(i)
          dzs(i)=0.d+0
      enddo
c
      do i=1,n2
          ru(i)=alpha
          rv(i)=beta
          rp(i)=gamma
          u0(i)=x0(i)
          v0(i)=x0(i+n2)
          ph0(i)=x0(i+2*n2)
          dxz(i)=udx
          dxz(i+n2)=vdx
          dxz(i+2*n2)=pdx
c         e1u(i)=0.d+0
c         e1v(i)=0.d+0
c         e1p(i)=0.d+0
          eu(i)=0.d+0
          ev(i)=0.d+0
          ep(i)=0.d+0
      enddo
c
c --- iteration de lagrangien augmente
c
      itlag=0
      itmax=1
  100 itlag=itlag+1
      if (itlag.gt.itmax) goto 9999
c
c --- nouveau point de depart
c     prend la valeur dans le fichier `fich' lorsqu'on fait
c         u1mt1 fich
c     fich doit contenir 1875 valeurs formatees
c
c     call getarg (1,fich)
c     open (nfich,file=fich,status='old',access='sequential',
c    &      form='formatted')
c     read (nfich,*) (x(i),i=1,1875)
c
c --- appel simulateur
c
      indic=4
      call u1mt1s (indic,n3,x,ff,g,izs,rzs,dzs)
c
c --- appel optimiseur
c     normtype can be 'sup' for sup-norm, 'two' for 2-norm, 'dfn' for the norm defined by prosca
c
      dxmin=dxz(1)
      df1=ff
      epsrel=1.d-10
      epsf=epsrel*df1
      zero=1.d-6
      nemax=21
      normtype = 'dfn'
      impres=5
      io=6
c
      iter=15
      nsim=10000
      imode(1)=0
      imode(2)=0
      imode(3)=0
c
      ntrav=4*n3+10*(2*n3+1)
      if (ntrav.gt.nmtrav) then
          write (6,900)
  900     format ("u1mt1: ntrav too large")
          stop
      endif
      reverse=0
      call m1qn3 (u1mt1s, euclid, ctonbe, ctcabe, n3, x, ff, g, dxmin,
     &            df1, epsrel, normtype, impres, io, imode, omode, iter,
     &            nsim, iz, trav, ntrav, reverse, indic, izs, rzs, dzs)
c
c --- rappel optimiseur
c
c     epsrel=1.d-4
c     mode=2
c     iter=5
c     nsim=3000
c     do 888 i=1876,ntrav
c        trav(i)=0.d+0
c 888 continue
c     call m1qn3 (u1mt1s,euclid,ctonbe,ctcabe,n3,x,ff,g,dxmin,df1,
c    /            epsrel,impres,io,mode,iter,nsim,iz,trav,ntrav,izs,rzs,
c    /            dzs)
c
c --- rappel optimiseur
c
c     epsrel=1.d-4
c     mode=2
c     iter=30
c     nsim=3000
c     do 889 i=1876,ntrav
c        trav(i)=0.d+0
c 889 continue
c     call m1qn3 (u1mt1s,euclid,ctonbe,ctcabe,n3,x,ff,g,dxmin,df1,
c    /            epsrel,impres,io,mode,iter,nsim,iz,trav,ntrav,izs,rzs,
c    /            dzs)
c     if (mode.ne.2) then
c         write (6,9998) (x(i),i=1,n3)
c         write (6,9997)
c         write (6,9998) (g(i),i=1,n3)
c9997     format (x)
c9998     format (3e22.15)
c     endif
c     ntrav=26260
c     call conmin (u1mt1s,euclid,n3,x,ff,g,dxmin,df1,epsrel,
c    /            impres,io,mode,iter,nsim,trav,ntrav,izs,rzs,dzs)
c
      indic=4
      call u1mt1s (indic,n3,x,ff,g,izs,rzs,dzs)
      call euclid (n2,eu,eu,sp,izs,rzs,dzs)
      au=sqrt(sp)/dble(n2)
      call euclid (n2,ev,ev,sp,izs,rzs,dzs)
      av=sqrt(sp)/dble(n2)
      call euclid (n2,ep,ep,sp,izs,rzs,dzs)
      ap=sqrt(sp)/dble(n2)
      ec=sqrt(au+av+ap)/dble(n3)
c     if(au.le.epsu.and.av.le.epsv.and.ap.le.epsp)goto 9999
c
c     update les multiplicateurs de lagrange
c
      do 150 i=1,n2
          rr=2.*cola*ru(i)*eu(i)
          dzs(i)=dzs(i)+rr
  150 continue
c
      do 180 i=1,n2
          j=i+n2
          rr=2.*cola*rv(i)*ev(i)
          dzs(j)=dzs(j)+rr
  180 continue
c
      do 200 i=1,n2
          j=i+2*n2
          rr=2.*cola*rp(i)*ep(i)
          dzs(j)=dzs(j)+rr
  200 continue
c
c     if (epsmult.lt.0.01) goto 9999
      if (omode.eq.1 .or. omode.eq.4 .or. omode.eq.5) go to 100
c
c --- sauvegarde
c
9999  continue
c
      stop
      end
c
c-----------------------------------------------------------------------
c
      subroutine u1mt1s (indic,n,x,ff,g,izs,rzs,dzs)
c
      implicit none
      integer indic, n, izs(*)
      real rzs(*)
      double precision x(n), ff, g(n), dzs(*)
c
c --- local variables
c
      integer n1, n2
c
      if (indic.eq.1) then
          write (6,'(/a/)') "Free call to u1mt1s"
          return
      endif
c
      n1=n/3+1
      n2=2*n1-1
      call u1mt1t(indic,x(1),x(n1),x(n2),ff,g,izs,dzs(1),
     1 dzs(n1),dzs(n2))
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine u1mt1t (indic,u,v,phi,fff,g,izs,ul,vl,pl)
c
      implicit none
      integer indic, izs(*)
      double precision u(*), v(*), phi(*), fff, g(*), ul(*), vl(*),
     &  pl(*)
c
c --- common
c
c     fc est le facteur de coriolis, fm celui d echelle
c
      integer n
      double precision dx, dy, fm, fc, u0, v0, ph0, x0, eu, ev, ep, ru,
     &  rv, rp
      common /u1mt1c/ dx, dy, fm(625), fc(625), u0(625), v0(625),
     &  ph0(625), x0(1875), eu(625), ev(625), ep(625), ru(625), rv(625),
     &  rp(625), n
c
c --- local variables
c
      integer i, j, k, m, iiz, ia, ib, nm1, nm2
      double precision poids, d1x, d2x, d1y, d2y, dxi2, dxi4, dyi2,
     &  dyi4, uu, vv, pp, uxp, uxm, uyp, uym, vxp, vxm, vyp, vym, pxp,
     &  pxm, pyp, pym, euu, evv, epp, zxc

c   calcul des ecarts au geostrophisme

      poids=1.d-2
      m=n*n
      nm2=n-2
      nm1=n-1
      d1x=0.5/dx
      d2x=d1x/2.d+0
      d1y=0.5/dy
      d2y=d1y/2.d+0
      dxi2=d1x
      dyi2=d1y
      dxi4=d2x
      dyi4=d2y
      fff=0.d+0
      do 1201 iiz=2,n-1
         uu=u(iiz)-u0(iiz)
         vv=v(iiz)-v0(iiz)
         pp=phi(iiz)-ph0(iiz)
         fff=fff+uu*uu+vv*vv+pp*pp
         uu=u(m+1-iiz)-u0(m+1-iiz)
         vv=v(m+1-iiz)-v0(m+1-iiz)
         pp=phi(m+1-iiz)-ph0(m+1-iiz)
         fff=fff+uu*uu+vv*vv+pp*pp*poids
 1201 continue
      do 100 j=1,nm2

c     calcul des contraintes sur une ligne

         ia=j*n+1
         ib=n*(j+1)
         uu=u(ia)-u0(ia)
         vv=v(ia)-v0(ia)
         pp=phi(ia)-ph0(ia)
         fff=fff+uu*uu+vv*vv+pp*pp
         uu=u(ib)-u0(ib)
         vv=v(ib)-v0(ib)
         pp=phi(ib)-ph0(ib)
         fff=fff+uu*uu+vv*vv+pp*pp*poids
         do 10 i=2,nm1
            k=j*n+i
            uxp=u(k+1)+u(k)
            uxm=u(k)+u(k-1)
            uyp=u(k+n)+u(k)
            uym=u(k)+u(k-n)
            vxp=v(k+1)+v(k)
            vxm=v(k)+v(k-1)
            vyp=v(k+n)+v(k)
            vym=v(k)+v(k-n)
            pxp=phi(k+1)+phi(k)
            pxm=phi(k)+phi(k-1)
            pyp=phi(k+n)+phi(k)
            pym=phi(k-n)+phi(k)
            euu=(uxp*uxp-uxm*uxm)*d2x+(uyp*vyp-uym*vym)*d2y
     1      -u(k)*((uxp-uxm)*d1x+(vyp-vym)*d1y)+(pxp-pxm)*d1x
            euu=euu*fm(k)-fc(k)*(vxp+vxm+vyp+vym)/8.d+0
            evv=(uxp*vxp-uxm*vxm)*d2x+(vyp*vyp-vym*vym)*d2y
     1      -v(k)*((uxp-uxm)*d1x+(vyp-vym)*d1y)+(pyp-pym)*d1y
            evv=evv*fm(k)+fc(k)*(uxp+uxm+uyp+uym)/8.d+0
            epp=fm(k)*((pxp*uxp-pxm*uxm)*d2x+(pyp*vyp-pym*vym)*
     1      d2y)
            uu=u(k)-u0(k)
            vv=v(k)-v0(k)
            pp=phi(k)-ph0(k)
            zxc=ru(k)*euu*euu
            zxc=zxc+rv(k)*evv*evv
            zxc=zxc+rp(k)*epp*epp
            fff=fff+uu*uu+vv*vv+pp*pp*poids+zxc
     1       +ul(k)*euu+vl(k)*evv+pl(k)*epp
            eu(k)=euu
            ev(k)=evv
            ep(k)=epp
   10    continue
  100 continue

c     fin du calcul du critere

c     calcul du gradient

      call u1mt1m (indic,u,v,phi,g(1),g(626),g(1251),izs,ul,vl,pl)

c     fin du calcul du gradient

      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine u1mt1m (indic,u,v,phi,u_d,v_d,phi_d,izs,ul,vl,pl)
c
      implicit none
      integer indic, izs(1)
      double precision u(625), v(625), phi(625), u_d(1), v_d(1),
     &  phi_d(1), ul(625), vl(625), pl(625)
c
c---
c
c     Code adjoint de u1mt1f.f
c
c     ATTENTION: ce code ne marche que si les variables
c                euu, evv et epp du common /u1mt1c/ ont bien ete
c                calculees lors du calcul de la fonction. Il faut pour
c                cela que les trois instructions suivantes
c                   eu(k)=euu
c                   ev(k)=evv
c                   ep(k)=epp
c                soient presentes dans u1mt1f.f qui doit etre execute
c                avant u1mt1m.
c
c     u_d, v_d, phi_d: sont les valeurs duales de u, v et phi (c est a
c        dire les derivees partielles de fff par rapport a u, v et phi)
c
c     ANALYSE DU CODE PRIMAL (u1mt1f)
c
c     constantes
c        u0,v0,p0
c        ru,rv,rp
c        ul,vl,pl
c        fm,fc
c     variables independantes
c        u,v,phi
c     variable dependante
c        fff
c     variables intervenant NON LINEAIREMENT dans le calcul
c        uu,vv,pp: recalcule
c        uxp,uxm,uyp,uym: recalcule
c        vxp,vxm,vyp,vym: recalcule
c        pxp,pxm,pyp,pym: recalcule
c        euu,evv,epp: memorise dans /u1mt1c/
c
c     La valeur duale d une variable "v" est notee "v_d"
c     les variables duales sont notees avec _d
c
c---
c
c --- common
c
c     fc est le facteur de coriolis, fm celui d echelle
c
      integer n
      double precision dx, dy, fm, fc, u0, v0, ph0, x0, euu, evv, epp,
     &  ru, rv, rp
      common /u1mt1c/ dx, dy, fm(625), fc(625), u0(625), v0(625),
     &  ph0(625), x0(1875), euu(625), evv(625), epp(625), ru(625),
     &  rv(625), rp(625), n

c --- variables locales

      integer m, nm1, nm2, i, j, k, ia, ib, iiz, no
      double precision r, rr
      double precision poids, d1x, d1y, d2x, d2y, r1x, r1y, r2x, r2y
      double precision uu, vv, pp
      double precision uxp, uxm, uyp, uym
      double precision vxp, vxm, vyp, vym
      double precision pxp, pxm, pyp, pym

c --- valeurs duales

      double precision fff_d,  euu_d,  evv_d,  epp_d,  uu_d,  vv_d,
     &  pp_d
      double precision zxc_d
      double precision uxp_d,  vxp_d,  pxp_d,  uyp_d,  vyp_d,  pyp_d
      double precision uxm_d,  vxm_d,  pxm_d,  uym_d,  vym_d,  pym_d

c --- initialisation des constantes

      poids=1.d-2
      m=n*n
      nm2=n-2
      nm1=n-1
      d1x=0.5/dx
      d2x=d1x/2.d+0
      d1y=0.5/dy
      d2y=d1y/2.d+0
      no=27

c     --- initialisation des valeurs duales

      fff_d = 1.d+0
      euu_d = 0.d+0
      evv_d = 0.d+0
      epp_d = 0.d+0
      uu_d = 0.d+0
      vv_d = 0.d+0
      pp_d = 0.d+0

      zxc_d = 0.d+0

      pxm_d = 0.d+0
      pxp_d = 0.d+0
      pym_d = 0.d+0
      pyp_d = 0.d+0
      uxm_d = 0.d+0
      uxp_d = 0.d+0
      uym_d = 0.d+0
      uyp_d = 0.d+0
      vxm_d = 0.d+0
      vxp_d = 0.d+0
      vym_d = 0.d+0
      vyp_d = 0.d+0

      do 1 i=1,m
         u_d(i) = 0.d+0
         v_d(i) = 0.d+0
         phi_d(i) = 0.d+0
    1 continue

c --- c'est parti

      do 100 j=nm2,1,-1
         do 10 i=nm1,2,-1
            k=j*n+i

c        re-calcul de certaines quantites

            uu=u(k)-u0(k)
            vv=v(k)-v0(k)
            pp=phi(k)-ph0(k)

c        code adjoint

            uu_d = uu_d + 2.*uu*fff_d
            vv_d = vv_d + 2.*vv*fff_d
            pp_d = pp_d + 2.*pp*poids*fff_d
            zxc_d = zxc_d + fff_d
            euu_d = euu_d + ul(k)*fff_d
            evv_d = evv_d + vl(k)*fff_d
            epp_d = epp_d + pl(k)*fff_d

            epp_d = epp_d + 2.*epp(k)*rp(k)*zxc_d
            evv_d = evv_d + 2.*evv(k)*rv(k)*zxc_d
            euu_d = euu_d + 2.*euu(k)*ru(k)*zxc_d
            zxc_d = 0.d+0

            phi_d(k) = phi_d(k) + pp_d
            pp_d=0.d+0
            v_d(k) = v_d(k) + vv_d
            vv_d=0.d+0
            u_d(k) = u_d(k) + uu_d
            uu_d=0.d+0
            
c        re-calcul de certaines quantites

            uxp=u(k+1)+u(k)
            uxm=u(k)+u(k-1)
            uyp=u(k+n)+u(k)
            uym=u(k)+u(k-n)
            vxp=v(k+1)+v(k)
            vxm=v(k)+v(k-1)
            vyp=v(k+n)+v(k)
            vym=v(k)+v(k-n)
            pxp=phi(k+1)+phi(k)
            pxm=phi(k)+phi(k-1)
            pyp=phi(k+n)+phi(k)
            pym=phi(k-n)+phi(k)

            rr = (uxp-uxm)*d1x+(vyp-vym)*d1y

c        code adjoint

            r = fm(k)*d2x*epp_d
            pxp_d = pxp_d + uxp*r
            uxp_d = uxp_d + pxp*r
            pxm_d = pxm_d - uxm*r
            uxm_d = uxm_d - pxm*r
            r = fm(k)*d2y*epp_d
            pyp_d = pyp_d + vyp*r
            vyp_d = vyp_d + pyp*r
            pym_d = pym_d - vym*r
            vym_d = vym_d - pym*r
            epp_d=0.d+0

            r = fc(k)*evv_d/8.d+0
            uxp_d = uxp_d + r
            uxm_d = uxm_d + r
            uyp_d = uyp_d + r
            uym_d = uym_d + r
            evv_d = fm(k)*evv_d

            r1x = v(k)*d1x*evv_d
            r2x = d2x*evv_d
            r1y = v(k)*d1y*evv_d
            r2y = 2.*d2y*evv_d
            uxp_d = uxp_d + vxp*r2x - r1x
            vxp_d = vxp_d + uxp*r2x
            uxm_d = uxm_d - vxm*r2x + r1x
            vxm_d = vxm_d - uxm*r2x
            vyp_d = vyp_d + vyp*r2y - r1y
            vym_d = vym_d - vym*r2y + r1y
            r = d1y*evv_d
            pyp_d = pyp_d + r
            pym_d = pym_d - r
            v_d(k) = v_d(k) - rr*evv_d
            evv_d = 0.d+0

            r = fc(k)*euu_d/8.d+0
            vxp_d = vxp_d - r
            vxm_d = vxm_d - r
            vyp_d = vyp_d - r
            vym_d = vym_d - r
            euu_d = fm(k)*euu_d

            r1y = u(k)*d1y*euu_d
            r2y = d2y*euu_d
            uyp_d = uyp_d + vyp*r2y
            vyp_d = vyp_d + uyp*r2y - r1y
            uym_d = uym_d - vym*r2y
            vym_d = vym_d - uym*r2y + r1y
            r1x = u(k)*d1x*euu_d
            r2x = 2.*d2x*euu_d
            uxp_d = uxp_d + uxp*r2x - r1x
            uxm_d = uxm_d - uxm*r2x + r1x
            r = d1x*euu_d
            pxp_d = pxp_d + r
            pxm_d = pxm_d - r
            u_d(k) = u_d(k) - rr*euu_d
            euu_d = 0.d+0

            phi_d(k-n) = phi_d(k-n) + pym_d
            phi_d(k) = phi_d(k) + pym_d
            pym_d = 0.d+0
            phi_d(k+n) = phi_d(k+n) + pyp_d
            phi_d(k) = phi_d(k) + pyp_d
            pyp_d = 0.d+0
            phi_d(k-1) = phi_d(k-1) + pxm_d
            phi_d(k) = phi_d(k) + pxm_d
            pxm_d = 0.d+0
            phi_d(k+1) = phi_d(k+1) + pxp_d
            phi_d(k) = phi_d(k) + pxp_d
            pxp_d = 0.d+0

            v_d(k-n) = v_d(k-n) + vym_d
            v_d(k) = v_d(k) + vym_d
            vym_d = 0.d+0
            v_d(k+n) = v_d(k+n) + vyp_d
            v_d(k) = v_d(k) + vyp_d
            vyp_d = 0.d+0
            v_d(k-1) = v_d(k-1) + vxm_d
            v_d(k) = v_d(k) + vxm_d
            vxm_d = 0.d+0
            v_d(k+1) = v_d(k+1) + vxp_d
            v_d(k) = v_d(k) + vxp_d
            vxp_d = 0.d+0

            u_d(k-n) = u_d(k-n) + uym_d
            u_d(k) = u_d(k) + uym_d
            uym_d = 0.d+0
            u_d(k+n) = u_d(k+n) + uyp_d
            u_d(k) = u_d(k) + uyp_d
            uyp_d = 0.d+0
            u_d(k-1) = u_d(k-1) + uxm_d
            u_d(k) = u_d(k) + uxm_d
            uxm_d = 0.d+0
            u_d(k+1) = u_d(k+1) + uxp_d
            u_d(k) = u_d(k) + uxp_d
            uxp_d = 0.d+0

   10    continue

         ia=j*n+1
         ib=n*(j+1)

c     re-calcul de certaines quantites

         uu=u(ib)-u0(ib)
         vv=v(ib)-v0(ib)
         pp=phi(ib)-ph0(ib)

c     code adjoint

         uu_d = uu_d + 2.*uu*fff_d
         vv_d = vv_d + 2.*vv*fff_d
         pp_d = pp_d + 2.*pp*poids*fff_d

         phi_d(ib) = phi_d(ib) + pp_d
         pp_d = 0.d+0
         v_d(ib) = v_d(ib) + vv_d
         vv_d = 0.d+0
         u_d(ib) = u_d(ib) + uu_d
         uu_d = 0.d+0

c     re-calcul de certaines quantites

         uu=u(ia)-u0(ia)
         vv=v(ia)-v0(ia)
         pp=phi(ia)-ph0(ia)

c     code adjoint

         uu_d = uu_d + 2.*uu*fff_d
         vv_d = vv_d + 2.*vv*fff_d
         pp_d = pp_d + 2.*pp*fff_d

         phi_d(ia) = phi_d(ia) + pp_d
         pp_d = 0.d+0
         v_d(ia) = v_d(ia) + vv_d
         vv_d = 0.d+0
         u_d(ia) = u_d(ia) + uu_d
         uu_d = 0.d+0

  100 continue

      do 1201 iiz=n-1,2,-1

c     re-calcul de certaines quantites

         uu=u(m+1-iiz)-u0(m+1-iiz)
         vv=v(m+1-iiz)-v0(m+1-iiz)
         pp=phi(m+1-iiz)-ph0(m+1-iiz)

c     code adjoint

         uu_d = uu_d + 2.*uu*fff_d
         vv_d = vv_d + 2.*vv*fff_d
         pp_d = pp_d + 2.*pp*poids*fff_d

         phi_d(m+1-iiz) = phi_d(m+1-iiz) + pp_d
         pp_d = 0.d+0
         v_d(m+1-iiz) = v_d(m+1-iiz) + vv_d
         vv_d = 0.d+0
         u_d(m+1-iiz) = u_d(m+1-iiz) + uu_d
         uu_d = 0.d+0

c     re-calcul de certaines quantites

         uu=u(iiz)-u0(iiz)
         vv=v(iiz)-v0(iiz)
         pp=phi(iiz)-ph0(iiz)

c     code adjoint

         uu_d = uu_d + 2.*uu*fff_d
         vv_d = vv_d + 2.*vv*fff_d
         pp_d = pp_d + 2.*pp*fff_d

         phi_d(iiz) = phi_d(iiz) + pp_d
         pp_d = 0.d+0
         v_d(iiz) = v_d(iiz) + vv_d
         vv_d = 0.d+0
         u_d(iiz) = u_d(iiz) + uu_d
         uu_d = 0.d+0

 1201 continue

      fff_d = 0.d+0

c --- Controle du bon codage: Toutes les valeurs duales
c     (autres que u_d, v_d et phi_d) doivent etre nulles.

c     if (fff_d.ne.0.) print *,"fff_d non nul"
c     if (euu_d.ne.0.) print *,"euu_d non nul"
c     if (evv_d.ne.0.) print *,"evv_d non nul"
c     if (epp_d.ne.0.) print *,"epp_d non nul"
c     if (uu_d.ne.0.) print *,"uu_d non nul"
c     if (vv_d.ne.0.) print *,"vv_d non nul"
c     if (pp_d.ne.0.) print *,"pp_d non nul"

c     if (zxc_d.ne.0.) print *,"zxc_d non nul"

c     if (pxm_d.ne.0.) print *,"pxm_d non nul"
c     if (pxp_d.ne.0.) print *,"pxp_d non nul"
c     if (pym_d.ne.0.) print *,"pym_d non nul"
c     if (pyp_d.ne.0.) print *,"pyp_d non nul"
c     if (uxm_d.ne.0.) print *,"uxm_d non nul"
c     if (uxp_d.ne.0.) print *,"uxp_d non nul"
c     if (uym_d.ne.0.) print *,"uym_d non nul"
c     if (uyp_d.ne.0.) print *,"uyp_d non nul"
c     if (vxm_d.ne.0.) print *,"vxm_d non nul"
c     if (vxp_d.ne.0.) print *,"vxp_d non nul"
c     if (vym_d.ne.0.) print *,"vym_d non nul"
c     if (vyp_d.ne.0.) print *,"vyp_d non nul"

      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine u1mt1l (ilec)
c
      implicit none
      integer ilec
c
c --- common
c
      integer n
      double precision dx, dy, fm, fc, u0, v0, ph0, x0, eu, ev, ep, ru,
     &  rv, rp
      common /u1mt1c/ dx, dy, fm(625), fc(625), u0(625), v0(625),
     &  ph0(625), x0(1875), eu(625), ev(625), ep(625), ru(625), rv(625),
     &  rp(625), n
c     common /u1mt1c/ dx,dy,fm(625),fc(625),u0(625),v0(625),ph0(625),
c    1 x0(1875),eu(625),ev(625),ep(625),ru(625),rv(625),rp(625),n
c
c --- local variables
c
      integer i, j, k, ii, num, in, n2, nc, ns
      double precision fco, pi, a, ph1, ph2, alpha, beta, ps1, ps2, ps4,
     &  psi, an, uu(55), vv(55), geo(55), gis, vit, xx(55), yy(55),
     &  lon, lat, ro, ro1, zz, ud, x, y, u7, v7, ge, pond, xxx, yyy,
     &  xcv, dist, vav, vau, vap, um, vm, pm, ua, va, pa, ecu, ecv, ecp
c
      fco=1854.d+0/3600.d+0
      pi=3.141593
      a=6356152.d+0
      ph1=pi/3.d+0
      ph2=pi/6.d+0
c ph1 et ph2 definissent le cone de projection
      beta=5.d+0*pi/180.d+0
      alpha=9.d+0*beta
c   alpha et beta coordonnees du centre
c
c
c calcul des coordonnees
c
      ps1=pi/2.-ph1
      ps2=pi/2.-ph2
      an=(log(cos(ph2)/cos(ph1)))/log(tan(ps2/2.)/tan(ps1/2.))
c    ns:nombre de stations
      ns=55
      do 10 i=1,ns
          read (ilec,*) num,lat,lon,in,geo(i),gis,vit
          if(in.eq.0)lon=-lon
          lat=lat*pi/180.d+0
          lon=lon*pi/180.d+0
          psi=pi/2.-lat
          ro1=a*sin(ps1)/an
          ro=ro1*(abs(tan(psi/2.)/tan(ps1/2.)))**an
          ps4=an*(lon-beta)
          xx(i)=ro1*sin(ps4)
          yy(i)=ro1*(abs(tan(pi/4.-alpha/2.)/tan(ps1/2.)))**an
     /          -ro*cos(ps4)
          zz=pi/2.-an*(lon-beta)+gis*pi/180.d+0
          uu(i)=-vit*cos(zz)*fco
          vv(i)=-vit*sin(zz)*fco
 10   continue
      ud=.1
c    unite de distance
      n2=n*n
      nc=(n+1)/2
      do 50i=0,n-1
      do60j=1,n
      ii=i*n+j
      x=(j-nc)*dx
      y=(i+1-nc)*dy
      u7=0
      v7=0
        fc(ii)=1.d-4
        fm(ii)=1.d+0
      ge=0.d+0
      pond=0.d+0
      do 160 k=1,ns
       xxx=(x-xx(k))/1.e5
       yyy=(y-yy(k))/1.e5
       xcv=xxx*xxx+yyy*yyy
      dist=xcv
      if(dist.le.ud)goto80
      pond=pond+1./dist
      u7=u7+uu(k)/dist
      v7=v7+vv(k)/dist
      ge=ge+geo(k)/dist
      goto 160
  80  u7=uu(k)
      v7=vv(k)
      ge=geo(k)
      pond=1.d+0
      goto 200
  160 continue
 200   x0(ii)=u7/pond
      x0(ii+n2)=v7/pond
      x0(ii+2*n2)=ge/pond
  60  continue
  50  continue
c   calcul des moyennes et variances
       vav=0.d+0
       vau=0.d+0
       vap=0.d+0
      um=0.d+0
      vm=0.d+0
      pm=0.d+0
      do 617 i=1,n2
      ua=x0(i)
      va=x0(i+n2)
      pa=x0(i+2*n2)
      um=um+ua
      vm=vm+va
      pm=pm+pa
      vau=vau+ua*ua
      vav=vav+va*va
      vap=vap+pa*pa
 617  continue
      vau=(vau-um*um/dble(n2))/dble(n2-1)
      vav=(vav-vm*vm/dble(n2))/dble(n2-1)
      vap=(vap-pm*pm/dble(n2))/dble(n2-1)
      um=um/dble(n2)
      vm=vm/dble(n2)
      pm=pm/dble(n2)
      ecu=sqrt(vau)
      ecv=sqrt(vav)
      ecp=sqrt(vap)
 9999 return
      end
c
c-----------------------------------------------------------------------
c
c     subroutine mupdts (sscale,inmemo,n,m,nrz)
c
c         arguments
c
c     logical sscale,inmemo
c     integer n,m,nrz
c----
c
c     This routine has to return:
c       m:      the numer of updates to form the approximate Hessien H,
c       inmemo = .true.  if the vectors y and s used to form H are
c                        stored in core memory,
c                .false. otherwise,
c     When inmemo=.false., the routine `ystbl', which stores and
c     restores (y,s) pairs, has to be rewritten.
c
c----
c
c     m=10
c     inmemo=.false.
c     return
c     end
c
c-----------------------------------------------------------------------
c
c     subroutine ystbl (store,ybar,sbar,n,j)
c----
c
c     This subroutine should store (if store = .true.) or restore
c     (if store = .false.) a pair (ybar,sbar) at or from position
c     j in memory. Be sure to have 1 <= j <= m, where m in the number
c     of updates specified by subroutine mupdts.
c
c     The subroutine is used only when the (y,s) pairs are not
c     stored in core memory in the arrays ybar(.,.) and sbar(.,.).
c     In this case, the subroutine has to be written by the user.
c
c----
c
c         arguments
c
c     logical store
c     integer n,j
c     real ybar(n),sbar(n)
c
c     real ytab,stab
c     common /ystab/ytab(1875,10),stab(1875,10)
c
c     if (store) then
c         do 10 i=1,n
c             ytab(i,j)=ybar(i)
c             stab(i,j)=sbar(i)
c  10     continue
c     else
c         do 20 i=1,n
c             ybar(i)=ytab(i,j)
c             sbar(i)=stab(i,j)
c  20     continue
c     endif
c
c     return
c     end
