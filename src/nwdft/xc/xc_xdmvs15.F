c   
c$Id: xc_xmvs15.F 21740 2012-01-11 00:25:15Z edo $
c
#include "dft2drv.fh"
c     Made Very Simple (MVS) functional (Exchange part only)
c           META GGA
C         utilizes ingredients:
c                              rho   -  density
c                              delrho - gradient of density
c                              tau - K.S kinetic energy density
c                              tauW - von Weiszacker kinetic energy density
c                              tauU - uniform-gas KE density
c     References:
c     J. Sun, J.P. Perdew, A. Ruzsinszky
c     PNAS 2015, 112, 685-689
c     DOI: 10.1073/pnas.1423145112 

      Subroutine xc_xdmvs15(tol_rho, fac,  rho, delrho, 
     &                     Amat, Cmat, nq, ipol, Ex, 
     &                     qwght, ldew, func, tau,Mmat,
     &                     laprho,Lmat,whichk)
      implicit none
c      
      double precision fac, Ex
      integer nq, ipol
      logical ldew
      double precision func(*)  ! value of the functional [output]
c
c     Charge Density & Its Cube Root
c
      double precision rho(nq,ipol*(ipol+1)/2)
c
c     Charge Density Gradient
c
      double precision delrho(nq,3,ipol)
c
c     Charge Density Laplacian
c
      double precision laprho(nq,ipol*(ipol+1)/2)
c
c     Quadrature Weights
c
      double precision qwght(nq)
c
c     Sampling Matrices for the XC Potential & Energy
c
      double precision Amat(nq,ipol), Cmat(nq,*), Lmat(nq,ipol)
c
c     kinetic energy density   or  tau
c
      double precision tau(nq,ipol), Mmat(nq,*)
      double precision tol_rho
      character*(*) whichk
c
      integer ispin,cmatpos
c
      if (ipol.eq.1 )then
c     
c     SPIN-RESTRICTED
c     Ex = Ex[n]
c
         call xc_xdmvs15_cs(tol_rho, fac,  rho, delrho, 
     &                     Amat, Cmat, nq, Ex, 1d0,
     &                     qwght, ldew, func, tau,Mmat,
     &                     laprho, Lmat, whichk)
      else
c     
c     SPIN-UNRESTRICTED
c     Ex = Ex[2n_up]/2 + Ex[2n_down]/2

         do ispin=1,2
            if (ispin.eq.1) cmatpos=D1_GAA
            if (ispin.eq.2) cmatpos=D1_GBB
            call xc_xdmvs15_cs(tol_rho, fac,  
     R           rho(1,ispin+1), delrho(1,1,ispin), 
     &           Amat(1,ispin), Cmat(1,cmatpos), 
     &           nq, Ex, 2d0,
     &           qwght, ldew, func, 
     T           tau(1,ispin),Mmat(1,ispin),
     &           laprho(1,ispin+1),Lmat(1,ispin),whichk)
         enddo

      endif
      return
      end
      Subroutine xc_xdmvs15_cs(tol_rho, fac,  rho, delrho, 
     &                     Amat, Cmat, nq, Ex, facttwo,
     &                     qwght, ldew, func, tau,Mmat,
     &                     laprho, Lmat, whichk)
      implicit none
c      
      double precision fac, Ex
      integer nq
      logical ldew
      double precision func(*)  ! value of the functional [output]
c
c     Charge Density & Its Cube Root
c
      double precision rho(*)
c
c     Charge Density Gradient
c
      double precision delrho(nq,3)
c
c     Charge Density Laplacian
c
      double precision laprho(*)
c
c     Quadrature Weights
c
      double precision qwght(nq)
c
c     Sampling Matrices for the XC Potential & Energy
c
      double precision Amat(nq), Cmat(nq), Lmat(nq)
c
c     kinetic energy density   or  tau
c
      double precision tau(nq,*), Mmat(nq)
c
      double precision facttwo ! 2 for o.s. 1 for c.s.
c
      character*(*) whichk
c
      integer n
      double precision tol_rho, pi
      double precision rrho, rho43, rho13, rho23, rho83
      double precision tauN, tauW, tauU, gamma
      double precision p, a, z, rz, g2, g, uge
      double precision F83, F23, F53, F43, F13, F18, F5, F8
      double precision afact2, Ax, rhoval, Pconst
      double precision rK0, rB, rE1, rC1

c     functional derivatives below FFFFFFFFFFFF

      double precision derivr, derivg, derivl
      double precision FaNum, BFaDen, FaDen, Fa
      double precision dFaNumda, dFaDenda, dBFaDenda, dFada
      double precision BFxDen, FxDen, FxNum
      double precision dFxNumda, dFxDendp
      double precision dadp, dadz, dpdg, dpdr, dzdr, dzdg, dzdt, dadt
      double precision Fx, dFxdp, dFxdz, dFxdr, dFxdg, dFxdl, dFxda
      double precision dtdr, dtdg, dtdl
      double precision dadr, dadg, dadl
      
c     functional derivatives above FFFFFFFFFFFF
      
      parameter (F5=5d0, F8=8d0)
      parameter (F43=4.d0/3.d0, F13=1.d0/3.d0)
      parameter (F83=8.d0/3.d0, F23=2.d0/3.d0)
      parameter (F18=1.d0/8.d0, F53=5.d0/3.d0)
      parameter (rK0=0.174d0, rB=0.0233d0,
     &           rE1=-1.66651d0, uge=10.d0/81.d0)
c
      pi=acos(-1d0)
      Ax = (-0.75d0)*(3d0/pi)**F13
      Pconst = (3.d0*pi**2)**F23
      afact2=1d0/facttwo
c
      rC1=(20.d0*rK0/(27.d0*uge))**(4d0)-(1.d0+rE1)*(1.d0+rE1) ! gfortran error if this is defined as a parameter
c
      do n = 1, nq
         if (rho(n).ge.tol_rho) then

           tauN=0d0
           if (whichk(1:2).eq.'pc') then
             call ts_mgga(tol_rho, rho(n), delrho(n,1:3),
     &                    laprho(n),dtdr,dtdg,dtdl,tauN,facttwo)
           else if (whichk(1:2).eq.'cr') then
             call ts_cr(tol_rho, rho(n), delrho(n,1:3),
     &                  laprho(n),dtdr,dtdg,dtdl,tauN,facttwo)
           else if (whichk(1:3).eq.'tfl') then
             call ts_tfl(tol_rho, rho(n), delrho(n,1:3),
     &                   laprho(n),dtdr,dtdg,dtdl,tauN,facttwo)
           else if (whichk(1:4).eq.'tanh') then
             call ts_tanh(tol_rho, rho(n), delrho(n,1:3),
     &                    laprho(n),dtdr,dtdg,dtdl,tauN,facttwo)
           endif


            rhoval=rho(n)*facttwo
            rho43 = rhoval**F43  ! rho^4/3
            rrho  = 1d0/rhoval   ! reciprocal of rho
            rho13 = rho43*rrho 
            rho23 = rhoval**F23
            rho83 = rhoval**F83
            
            g2 = delrho(n,1)*delrho(n,1) +
     &           delrho(n,2)*delrho(n,2) +
     &           delrho(n,3)*delrho(n,3)

            if (dsqrt(g2).gt.tol_rho)then
               g2 = g2 *facttwo*facttwo
               g = dsqrt(g2)
            else
               g  = tol_rho
               g2 = tol_rho*tol_rho
            endif

            tauN = tauN*facttwo
            tauW = F18*g2*rrho
            tauU = 0.3d0*Pconst*rhoval**F53

c     
c     Evaluate the Fx
c     
            p=g2/(4d0*Pconst*rho83)
c
c            a=(tauN-tauW)/tauU
c
            a = (tauN-tauW)/tauU
            if(a.lt.0d0) a=0d0
 
            FaNum = (1d0 - a)
            BFaDen = (1d0 + rE1*a*a)**2d0 + rC1*a**4d0
            FaDen = BFaDen**0.25d0
            Fa = Fanum / FaDen

            BFxDen = 1d0 + rB*p*p
            FxDen = BFxDen**F18
            FxNum = 1.d0 + rK0*Fa
            Fx = FxNum / FxDen

            Ex = Ex + Fx*Ax*rho43*qwght(n)*fac*afact2
            if (ldew)  func(n)= func(n) + Fx*Ax*rho43*fac*afact2

c     functional derivatives FFFFFFFFFFFFFFFFFFFFFFFFFFFF
         
            dpdr = -F83*p*rrho
            dpdg = 1d0/(4d0*Pconst*rho83)

            dadt = 1d0/tauU
            dadg = -0.125d0/(tauU*rhoval)
            dadr = F13*(F8*tauW - F5*tauN)/(tauU*rhoval)

            dadr = dadr + dadt*dtdr
            dadg = dadg + 0.5d0*dadt*dtdg
            dadl = dadt*dtdl

            dFaNumda = -1d0
            dBFaDenda = 2d0*(1d0 + rE1*a*a)*2d0*rE1*a + 4d0*rC1*a**3d0
            dFaDenda = 0.25d0*(FaDen/BFaDen)*dBFaDenda
            dFada = (dFaNumda * FaDen - FaNum*dFaDenda)/(FaDen*FaDen)

            dFxda = rk0*dFada/FxDen
            dFxdp = -0.25d0*FxNum*rB*p/(FxDen*BFxDen)

            dFxdr = dFxda * dadr + dFxdp * dpdr
            dFxdg = dFxda * dadg + dFxdp * dpdg
            dFxdl = dFxda * dadl

            derivr = F43*Ax*rho13*Fx + Ax*rho43*dFxdr
            derivg = Ax*rho43*dFxdg
            derivl = Ax*rho43*dFxdl

            Amat(n) = Amat(n) + derivr*fac
c     
c     4x factor comes from gamma_aa = gamma_total/4
c     
            Cmat(n)=  Cmat(n) + 2.0d0*derivg*fac
            Lmat(n)=  Lmat(n) + derivl*fac*qwght(n)
         endif
      enddo
      return
      end

      Subroutine xc_xdmvs15_d2()
      call errquit(' xmvs15: d2 not coded ',0,0)
      return
      end



