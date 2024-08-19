c    SCAN correlation functional 
c           META GGA
C         utilizes ingredients:
c                              rho   -  density
c                              delrho - gradient of density
c                              tau (tauN)- K.S kinetic energy density

      Subroutine xc_cdscan(tol_rho, cfac, rho, delrho, Amat, Cmat,
     &                     nq, ipol, Ec, qwght, ldew, func,
     &                     tau, Mmat, laprho, Lmat, whichk)     
     

c
c$Id: xc_cm06.F 26990 2015-04-15 22:00:47Z d3y133 $
c
c

      implicit none
c
#include "errquit.fh"
#include "dft2drv.fh"
c
c     Input and other parameters
c
      character(*) whichk
      integer lnq,ispin
      integer ipol, nq
      double precision dummy(1)
      double precision cfac(*)
      logical ldew
      logical lfac, nlfac,lfacl,nlfacl
      double precision func(*)
      double precision fac,facl
      double precision lqwght
      double precision tol_rho
c
c     Threshold parameters
c
      double precision TOLL
      double precision EXPTOL
      double precision EPS
      double precision eps_argexp
      double precision thr1,thr2
      parameter (TOLL = 1.0D-40)
      parameter (EXPTOL = 40.0D0)
      parameter (EPS = 1.0D-8)
      parameter (eps_argexp = 1.0D-14)
      parameter (thr1=0.996d0,thr2=1.004d0)
c
c     Correlation energy
c
      double precision Ec
c
c     Charge Density 
c
      double precision rho(nq,ipol*(ipol+1)/2)
c
c     Charge Density Gradient
c
      double precision delrho(nq,3,ipol), gammaval, gam12
c
c     Kinetic Energy Density
c
      double precision tau(nq,ipol), tauN
c
c     Charge Density Laplacian
c
      double precision laprho(nq,ipol*(ipol+1)/2)
c
c     Quadrature Weights
c
      double precision qwght(nq)
c
c     Sampling Matrices for the XC Potential
c
      double precision Amat(nq,ipol), Cmat(nq,*)
      double precision Mmat(nq,*), Lmat(nq,*)
c
c     Intermediate derivatives, results, etc.
c
      integer n
      double precision ntot,n13,n53,n83,tautot
      double precision dn2,p,rs,rs12,drsdn,dpdg,dpdn
      double precision epsc,depscdrs,depscdzeta
      double precision zeta,omz,opz
      double precision phi,phi2,phi3,dphidzeta
      double precision t2,dt2dp,dt2drs,dt2dzeta
      double precision BETA,dBETAdrs
      double precision A,dAdrs,dAdzeta
      double precision At2,Gaux,Gat2
      double precision w1fac,expw1,w1,dw1drs,dw1dzeta
      double precision arg1,darg1dp,darg1drs,darg1dzeta
      double precision H1,dH1dp,dH1drs,dH1dzeta
      double precision Ec1,dEc1drs,dEc1dp,dEc1dzeta,dEc1dn,dEc1dg
      double precision epsc0,depsc0drs,depsc0dn
      double precision dx,ddxdzeta
      double precision gc,dgcdzeta
      double precision ginf,dginfdp
      double precision w0fac,expw0,w0,dw0drs
      double precision arg0,darg0dp,darg0drs
      double precision H0,dH0dp,dH0drs
      double precision Ec0,dEc0dp,dEc0drs,dEc0dzeta,dEc0dn,dEc0dg
      double precision ds,ddsdzeta
      double precision tueg,tvw
      double precision alpha,dalphadzeta,dalphadn,dalphadg,dalphadt
      double precision oma,oma2
      double precision fca,dfcada,dfcadg,dfcadzeta,dfcadn
      double precision exp5,exp6
      double precision vcpol,vcn

      double precision ge2,ge4,mge4,ft,zpc,fpc,bfpc,q,lapval
      double precision exp8,exp9,factor,ftw,mge4den
      double precision dge2dp,dge2dq,dge4dp,dge4dq,dmge4dp,dmge4dq
      double precision dtauudr,dfpcdz,dfpc,dftdz,dftdp,dftdq
      double precision dqdl,dqdn,dtdn(ipol),dtdg(ipol),dtdl(ipol)
      double precision dadp,dadq,dadl
      double precision dalphadl,dtuegdn,dfcadt

      double precision GAMMA,BETAzero,pi,p14a,p14b,p14f,ckf,ckf2
      parameter (GAMMA = 0.03109069086965489503494086371273d0)
      parameter (BETAzero = 0.06672455060314922d0)
      parameter (p14a=0.1d0,p14b=0.1778d0)

      double precision b1c,b2c,b3c,c1c,c2c,dc,dxc,xi
      parameter (b1c=0.0285764d0,b2c=0.0889d0,b3c=0.125541d0)
      parameter (c1c=0.64d0,c2c=1.5d0,dc=0.7d0,dxc=2.3631d0)
      parameter (xi=0.12802585262625815d0)

      double precision F4,F13,F23,F43,F53,F83,F14,F8,F5,F16
      parameter (F4=4d0,F13=1d0/3d0,F23=F13+F13,F43=F13+1d0)
      parameter (F53=F23+1d0,F83=F53+1d0,F14=0.25d0)
      parameter (F8=8d0,F5=5d0,F16=1d0/6d0)
c      parameter (amgga=1.86962d0, bmgga=0.228503d0)
c      parameter (amgga=1.0860d0, bmgga=0.45d0)
c      parameter (amgga=1.2880d0, bmgga=0.340d0)
c      parameter (amgga=1.7990d0, bmgga=0.2550d0)
c      parameter (amgga=3.9300d0, bmgga=0.080d0)
c      parameter (amgga=1.3770d0, bmgga=0.300d0)
c      parameter (amgga=0.5389d0, bmgga=3.0d0)
c      parameter (amgga=1.203d0, bmgga=0.47d0)

c
c     ======> BOTH SPIN-RESTRICETED AND UNRESTRICTED <======
c
      Pi = dacos(-1d0)
      p14f = (3d0/(4d0*Pi))**F13
      ckf = (3d0*Pi*Pi)**F13
      ckf2 = ckf*ckf

      do 20 n = 1, nq
c
         ntot=rho(n,1)

         if (ntot.lt.tol_rho) goto 20

         n13=ntot**F13
         n53=ntot**F53
         n83=ntot**F83
c
         if (ipol.eq.1) then
           zeta = 0d0
         else
           zeta = (rho(n,2) - rho(n,3))/ntot      
           if (zeta.lt.-1d0) zeta=-1d0
           if (zeta.gt. 1d0) zeta= 1d0
         endif
c
         if (ipol.eq.1) then
           dn2=delrho(n,1,1)**2 + delrho(n,2,1)**2 + delrho(n,3,1)**2
         else
           dn2=(delrho(n,1,1)+delrho(n,1,2))**2 +
     &         (delrho(n,2,1)+delrho(n,2,2))**2 +
     &         (delrho(n,3,1)+delrho(n,3,2))**2
         end if

         if (ipol.eq.1) then
           tautot=tau(n,1)
         else
           tautot=tau(n,1)+tau(n,2)
         endif

         p=dn2/(F4*ckf2*n83)
         dpdg = 1d0/(F4*ckf2*n83)
         dpdn = -F83*p/ntot

         rs=p14f/n13
         drsdn=-F13*rs/ntot
         rs12=dsqrt(rs)
c
         call lsdac(tol_rho,rs,zeta,epsc,depscdrs,depscdzeta,dummy,
     &              dummy,dummy)
c
         opz = 1d0 + zeta
         omz = 1d0 - zeta
         phi = 0.5d0*(opz**F23 + omz**F23)
         phi2 = phi*phi
         phi3 = phi2*phi
c
         BETA = BETAzero*(1d0 + p14a*rs)/(1d0 + p14b*rs)
         dBETAdrs = BETAzero*(p14a - p14b)/(1d0 + p14b*rs)**2
c
         w1fac = epsc/(GAMMA*phi3)
         expw1 = dexp(-w1fac)
         w1 = expw1 - 1d0
         dw1drs = -expw1*depscdrs/(GAMMA*phi3)

         A = BETA/GAMMA/w1
         dAdrs = dBETAdrs/GAMMA/w1 - A*dw1drs/w1

         t2 = ckf2*p/(4d0*phi2*rs*2d0**F23)
         dt2dp = ckf2/(4d0*phi2*rs*2d0**F23)
         dt2drs = -t2/rs

         At2 = A*t2
         Gaux = 1d0 + 4d0*At2
         GAt2 = 1d0/dsqrt(dsqrt(Gaux))

         arg1 = 1d0 + w1*(1d0 - GAt2)
         darg1drs = dw1drs*(1d0 - GAt2) +
     &              Gat2*w1*(t2*dAdrs + A*dt2drs)/Gaux
         darg1dp = Gat2*w1*A*dt2dp/Gaux


         H1 = GAMMA*phi3*dlog(arg1)
         dH1dp = GAMMA*phi3*darg1dp/arg1
         dH1drs = GAMMA*phi3*darg1drs/arg1

         Ec1 = epsc + H1
         dEc1drs = depscdrs + dH1drs
         dEc1dp = dH1dp
c
c        --------------------------------------------------------------
c
         epsc0 = -b1c/(1d0 + b2c*rs12 + b3c*rs)
         depsc0drs = b1c*(b3c + 0.5d0*b2c/rs12)/
     &               (1d0 + b2c*rs12 + b3c*rs)**2
 
         dx = 0.5d0*(opz**F43 + omz**F43)
         gc = (1d0 - dxc*(dx - 1d0))*(1d0 - zeta**12)

         w0fac = epsc0/b1c
         expw0 = dexp(-w0fac)
         w0 = expw0 - 1d0
         dw0drs = -depsc0drs*expw0/b1c

         ginf = (1d0/(1d0 + 4d0*xi*p))**F14
         dginfdp = -xi*ginf/(1d0 + 4d0*xi*p)

         arg0 = 1d0 + w0*(1d0 - ginf)
         darg0drs = dw0drs*(1d0 - ginf)
         darg0dp = -w0*dginfdp

         H0 = b1c*dlog(arg0)
         dH0dp = b1c*darg0dp/arg0
         dH0drs = b1c*darg0drs/arg0

         Ec0 = (epsc0 + H0)*gc
         dEc0drs = gc*(depsc0drs + dH0drs)
         dEc0dp = gc*dH0dp
c
c        --------------------------------------------------------------
c
c        Kinetic energy density
c
         ds = 0.5d0*(opz**F53 + omz**F53)
         tueg = 0.3d0*ckf*ckf*ds*n53
         tvw = 0.125d0*dn2/ntot

         tautot = 0d0
         if (ipol.eq.1) then
             if (whichk(1:2).eq.'pc') then
               call ts_mgga(tol_rho,rho(n,1),delrho(n,1:3,1),
     &                      laprho(n,1),dtdn,dtdg,
     &                      dtdl,tautot,1d0)
             else if (whichk(1:2).eq.'cr') then
               call ts_cr(tol_rho,rho(n,1),delrho(n,1:3,1),
     &                    laprho(n,1),dtdn,dtdg,
     &                    dtdl,tautot,1d0)
             else if (whichk(1:3).eq.'tfl') then
               call ts_tfl(tol_rho,rho(n,1),delrho(n,1:3,1),
     &                     laprho(n,1),dtdn,dtdg,
     &                     dtdl,tautot,1d0)
             else if (whichk(1:4).eq.'tanh') then
               call ts_tanh(tol_rho,rho(n,1),delrho(n,1:3,1),
     &                      laprho(n,1),dtdn,dtdg,
     &                      dtdl,tautot,1d0)
             endif
         else if (ipol.eq.2) then
           do ispin=1,2
             if (whichk(1:2).eq.'pc') then
               call ts_mgga(tol_rho,rho(n,ispin+1),delrho(n,1:3,ispin),
     &                      laprho(n,ispin+1),dtdn(ispin),dtdg(ispin),
     &                      dtdl(ispin),tautot,2d0)
             else if (whichk(1:2).eq.'cr') then
               call ts_cr(tol_rho,rho(n,ispin+1),delrho(n,1:3,ispin),
     &                    laprho(n,ispin+1),dtdn(ispin),dtdg(ispin),
     &                    dtdl(ispin),tautot,2d0)
             else if (whichk(1:3).eq.'tfl') then
               call ts_tfl(tol_rho,rho(n,ispin+1),delrho(n,1:3,ispin),
     &                     laprho(n,ispin),dtdn(ispin),dtdg(ispin),
     &                     dtdl(ispin),tautot,2d0)
             else if (whichk(1:4).eq.'tanh') then
               call ts_tanh(tol_rho,rho(n,ispin+1),delrho(n,1:3,ispin),
     &                      laprho(n,ispin+1),dtdn(ispin),dtdg(ispin),
     &                      dtdl(ispin),tautot,2d0)
             endif
           enddo
         endif

         dtuegdn = F53*tueg/ntot

         alpha = (tautot - tvw)/tueg
         if (alpha.lt.0d0) alpha=0d0
         oma = 1d0 - alpha
         oma2 = oma*oma

         if (alpha.ge.thr1) then
           exp5 = 0d0
         else
           exp5 = dexp(-c1c*alpha/oma)
         end if
         if (alpha.le.thr2) then
           exp6 = 0d0
         else
           exp6 = dexp(c2c/oma)
         end if
         fca = exp5 - dc*exp6

         if (alpha.ge.thr1.and.alpha.le.thr2) then
           dfcada = 0d0
         else
           dfcada = -(c1c*exp5 + dc*exp6*c2c)/oma2
         end if
c
c        --------------------------------------------------------------
c
         Ec = Ec + ntot*(Ec1 + fca*(Ec0-Ec1))*qwght(n)
         if (ldew) func(n) = func(n) + ntot*(Ec1 + fca*(Ec0-Ec1))
c
         dalphadn = F13*(F8*tvw - F5*tautot)/(tueg*ntot)
         dalphadg = -0.125d0/(tueg*ntot)
         dalphadt = 1d0/tueg

         dEc1dn = dEc1drs*drsdn + dEc1dp*dpdn
         dEc0dn = dEc0drs*drsdn + dEc0dp*dpdn
         dfcadn = dfcada*dalphadn
         dfcadg = dfcada*dalphadg
         dfcadt = dfcada*dalphadt

         vcn = Ec1 + fca*(Ec0-Ec1) + 
     &         ntot*(dEc1dn + fca*(dEc0dn-dEc1dn) + dfcadn*(Ec0-Ec1))

         Amat(n,1) = Amat(n,1) + vcn + ntot*(dfcadt*(Ec0-Ec1)*dtdn(1))

         dEc1dg = dEc1dp*dpdg
         dEc0dg = dEc0dp*dpdg

         Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + ntot*(dEc1dg + 
     &   fca*(dEc0dg-dEc1dg) + dfcadg*(Ec0-Ec1) +
     &   dfcadt*(Ec0-Ec1)*dtdg(1))
         Cmat(n,D1_GAB) = Cmat(n,D1_GAB) + ntot*(dEc1dg + 
     &   fca*(dEc0dg-dEc1dg) + dfcadg*(Ec0-Ec1))*2d0

         Lmat(n,1) = Lmat(n,1) + ntot*dfcadt*dtdl(1)*(Ec0-Ec1)*qwght(n)

         if (ipol.eq.2) then
           Amat(n,2) = Amat(n,2) + vcn + ntot*dfcadt*(Ec0-Ec1)*dtdn(2)

           Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + ntot*(dEc1dg +
     &     fca*(dEc0dg-dEc1dg) + dfcadt*(Ec0-Ec1)*dtdg(2) +
     &     dfcadg*(Ec0-Ec1))
           
           Lmat(n,2) = Lmat(n,2) + ntot*dfcadt*dtdl(2)*(Ec0-Ec1)*
     &     qwght(n)

           if (omz.lt.tol_rho) then
             dphidzeta = 0.5d0*F23*(opz**F23/opz)
           else if (opz.lt.tol_rho) then
             dphidzeta = -0.5d0*F23*(omz**F23/omz)
           else
             dphidzeta = 0.5d0*F23*(opz**F23/opz - omz**F23/omz)
           end if

           dt2dzeta = -2d0*t2*dphidzeta/phi
           dw1dzeta = (3d0*w1fac*dphidzeta/phi - 
     &                 depscdzeta/(GAMMA*phi3))*expw1
           dAdzeta = -A*dw1dzeta/w1
           darg1dzeta = dw1dzeta*(1d0 - Gat2) + 
     &                  Gat2*w1*(t2*dAdzeta + A*dt2dzeta)/Gaux
           dH1dzeta = 3d0*H1*dphidzeta/phi + GAMMA*phi3*darg1dzeta/arg1
           dEc1dzeta = depscdzeta + dH1dzeta

           ddxdzeta = 0.5d0*F43*(opz**F13 - omz**F13)
           dgcdzeta = -dxc*ddxdzeta*(1d0 - zeta**12) -
     &                12d0*zeta**11*(1d0 - dxc*(dx - 1d0))
           dEc0dzeta = dgcdzeta*(epsc0 + H0)

           ddsdzeta = 0.5d0*F53*(opz**F23 - omz**F23)
           dalphadzeta = -alpha*ddsdzeta/ds
           dfcadzeta = dfcada*dalphadzeta

           vcpol = dEc1dzeta + dfcadzeta*(Ec0-Ec1) + 
     &             fca*(dEc0dzeta-dEc1dzeta)

           Amat(n,1) = Amat(n,1) + (1d0-zeta)*vcpol
           Amat(n,2) = Amat(n,2) - (1d0+zeta)*vcpol
         
         end if

20    continue
      end

      Subroutine xc_cdscan_d2()
      implicit none
      call errquit(' xc_cdscan: d2 not coded ',0,0)
      return
      end

