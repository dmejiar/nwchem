c   
c$Id: xc_xscan.F 21740 2012-01-11 00:25:15Z dmr $
c
#include "dft2drv.fh"
c     Strongly constrained and appropriately normed (SCAN) 
c     functional (Exchange part only)
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

      Subroutine xc_xdscan(tol_rho, fac,  rho, delrho, 
     &                    Amat, Cmat, nq, ipol, Ex, 
     &                    qwght, ldew, func, tau,Mmat,
     &                    laprho,Lmat,whichk)
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
c     Kinetic energy density
c
      if (ipol.eq.1 )then
c     
c     SPIN-RESTRICTED
c     Ex = Ex[n]
c
         call xc_xdscan_cs(tol_rho, fac,  rho, delrho, 
     &                    Amat, Cmat, nq, Ex, 1d0,
     &                    qwght, ldew, func, tau,Mmat,
     &                    laprho,Lmat,whichk)
      else
c     
c     SPIN-UNRESTRICTED
c     Ex = Ex[2n_up]/2 + Ex[2n_down]/2

         do ispin=1,2
            if (ispin.eq.1) cmatpos=D1_GAA
            if (ispin.eq.2) cmatpos=D1_GBB
            call xc_xdscan_cs(tol_rho, fac,  
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
      Subroutine xc_xdscan_cs(tol_rho, fac,  rho, delrho, 
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
      character*(*) whichk
c
      double precision facttwo ! 2 for o.s. 1 for c.s.
c
      integer n
      double precision tol_rho, pi
      double precision rrho, rho43, rho13, rho23, rho83, rho53
      double precision tauN, tauW, tauU, gamma
      double precision p, p14, a, z, rz, g2, uge
      double precision F83, F23, F53, F43, F13, F18, F8, F5
      double precision afact2, Ax, rhoval, Pconst
      double precision rH0, rK1, rA1, rC1, rC2, rD, rMu
      double precision rB1, rB2, rB3, rB4
      double precision oma, oma2
      double precision exp1, exp2, exp3, exp4, exp5
      double precision x1, x2, x
      double precision H, Hden, Hnum
      double precision G
      double precision thr1,thr2
c
      double precision q,lapval

c     functional derivatives below FFFFFFFFFFFF

      double precision derivr, derivg, derivt
      double precision Fa,dFada
      double precision dGdp, dHdp, dHdx, dHda
      double precision dxdp, dx1dp, dx2dp, dxda
      double precision dadz, dpdg, dpdr, dzdr, dzdg, dzdt
      double precision Fx, dFxda, dFxdp, dFxdz, dFxdr, dFxdg, dFxdt
      double precision dadg,dadt,dadr

      double precision dtauudr,dfpcdz,dfpc,dftdz,dftdp,dftdq
      double precision dqdl,dqdr,dtdr,dtdg,dtdl,dadp,dadq,dadl
      double precision dFxdl,derivl
      
c     functional derivatives above FFFFFFFFFFFF
      
      parameter (F43=4.d0/3.d0, F13=1.d0/3.d0)
      parameter (F83=8.d0/3.d0, F23=2.d0/3.d0)
      parameter (F18=1.d0/8.d0, F53=5.d0/3.d0)
      parameter (F5=5.d0, F8=8.d0)
      parameter (thr1=0.996d0, thr2=1.004d0)
      parameter (rH0=1.174d0, rK1=0.065d0, rA1=4.9479d0)
      parameter (rC1=0.667d0, rC2=0.8d0, rD=1.24d0)
c      parameter (amgga=1.86962d0, bmgga=0.228503d0)
c      parameter (amgga=1.0860d0, bmgga=0.45d0)
c      parameter (amgga=1.2880d0, bmgga=0.340d0)
c      parameter (amgga=1.7990d0, bmgga=0.2550d0)
c      parameter (amgga=3.9300d0, bmgga=0.080d0)
c      parameter (amgga=1.3770d0, bmgga=0.300d0)
c      parameter (amgga=0.5389d0, bmgga=3.0d0)
c      parameter (amgga=1.203d0, bmgga=0.47d0)

      rMu=10.d0/81.d0
      rB2=dsqrt(5913.d0/405000.d0)
      rB1=(511.d0/13500.d0)/(2.d0*rB2)
      rB3=0.5d0
      rB4=rMu*rMu/rK1-1606.d0/18225.d0-rB1*rB1
c
      pi=acos(-1d0)
      Ax = (-0.75d0)*(3d0/pi)**F13
      Pconst = (3.d0*pi**2)**F23
      afact2=1d0/facttwo
c
      do n = 1, nq
         if (rho(n).ge.tol_rho) then

           tauN = 0d0
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
            rho53 = rho23*rhoval
            rho83 = rhoval**F83
      
      
            g2 = delrho(n,1)*delrho(n,1) +
     &           delrho(n,2)*delrho(n,2) +
     &           delrho(n,3)*delrho(n,3)

            g2 = g2*facttwo*facttwo

            lapval = laprho(n)*facttwo

            tauW = F18*g2*rrho
            tauU = 0.3d0*Pconst*rho53
            
            tauN=facttwo*tauN
            
            p = g2/(4d0*Pconst*rho83)
            q = lapval/(4d0*Pconst*rho53)
c     
c     Evaluate the Fx
c     
            p14 = dsqrt(dsqrt(p))
c
            a = (tauN-tauW)/tauU
            if(a.lt.0d0)  a=0d0

            oma = 1d0 - a
            oma2 = oma*oma
            
            exp1 = dexp(-rB4/rMu*p)
            exp2 = dexp(-rB3*oma2)
            x1 = rMu*p*(1d0 + rB4/rMu*p*exp1)
            x2 = rB1*p + rB2*oma*exp2

            x = x1 + x2*x2

            Hden = rK1 + x
            Hnum = hden + rK1*x
            H = Hnum/Hden

            if (p14.lt.0.001d0) then
              exp3 = 0d0
            else
              exp3 = dexp(-rA1/p14)
            endif
            G = 1d0 - exp3

            if (a.ge.thr1) then
              exp4 = 0d0
            else
              exp4 = dexp(-rC1*a/oma)
            end if

            if (a.le.thr2) then
              exp5 = 0d0
            else
              exp5 = dexp(rC2/oma)
            end if

            Fa = exp4 - rD*exp5

            Fx = G*(H + Fa*(rH0 - H))

            Ex = Ex + Fx*Ax*rho43*qwght(n)*fac*afact2
            if (ldew)  func(n)= func(n) + Fx*Ax*rho43*fac*afact2

c     functional derivatives FFFFFFFFFFFFFFFFFFFFFFFFFFFF
         
            
            dpdr = -F83*p*rrho
            dpdg = 1.d0/(4d0*Pconst*rho83)

            dqdr = -F53*q*rrho
            dqdl = 1.d0/(4d0*Pconst*rho53)

            dtauUdr = F53*tauU*rrho

            dadt = 1d0/tauU
            dadg = -0.125d0/(tauU*rhoval)
            dadr = F13*(F8*tauW - F5*tauN)/(tauU*rhoval)

            dadr = dadr + dadt*dtdr
            dadg = dadg + dadt*dtdg*0.5d0
            dadl = dadt*dtdl

            if (p14.lt.0.002d0) then
              dGdp = 0d0
            else
              dGdp = -0.25d0*rA1*exp3/(p*p14)
            end if

            dx1dp = rMu + rB4*p*exp1*(2d0 - p*rB4/rMu)
            dx2dp = rB1
            dxdp = dx1dp + 2d0*x2*dx2dp
            dxda = 2d0*rB2*exp2*x2*(2d0*rB3*oma2 - 1d0)

            dHdx = (rK1/Hden)**2
            dHdp = dHdx*dxdp
            dHda = dHdx*dxda

            if ((a.ge.thr1).and.(a.le.thr2)) then
              dFada = 0d0
            else
              dFada = -(rC1*exp4 + rD*exp5*rC2)/oma2
            end if

            dFxdp = dGdp*(H + Fa*(rH0 - H)) + G*dHdp*(1d0 - Fa)
            dFxda = G*(dHda + dFada*(rH0 - H) - Fa*dHda)

            dFxdr = dFxda*dadr + dFxdp*dpdr
            dFxdg = dFxda*dadg + dFxdp*dpdg
            dFxdl = dFxda*dadl

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

      Subroutine xc_xdscan_d2()
      call errquit(' xscan: d2 not coded ',0,0)
      return
      end



