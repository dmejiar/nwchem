      Subroutine ts_mgga(tol_rho, rho, delrho, laprho,
     &                   dtdn, dtdg, dtdl, tau, facttwo)
      implicit none
c
#include "xc_kedf.fh"
c      
      double precision tau
c
c     Charge Density & Its Cube Root
c
      double precision rho
c
c     Charge Density Gradient
c
      double precision delrho(3)
c
c     Charge Density Laplacian
c
      double precision laprho
c
c     Sampling Matrices for the XC Potential & Energy
c
      double precision dtdn, dtdg, dtdl
c
      double precision facttwo ! 2 for o.s. 1 for c.s.
c
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
      double precision ge2,ge4,mge4,ft,zpc,fpc,bfpc,q,lapval
      double precision exp8,exp9,factor,ftw,mge4den

c     functional derivatives below FFFFFFFFFFFF

      double precision derivr, derivg, derivt
      double precision Fa,dFada
      double precision dGdp, dHdp, dHdx, dHda
      double precision dxdp, dx1dp, dx2dp, dxda
      double precision dadz, dpdg, dpdr, dzdr, dzdg, dzdt
      double precision Fx, dFxda, dFxdp, dFxdz, dFxdr, dFxdg, dFxdt
      double precision dadg,dadt,dadr

      double precision dge2dp,dge2dq,dge4dp,dge4dq,dmge4dp,dmge4dq
      double precision dtauudr,dfpcdz,dfpc,dftdz,dftdp,dftdq
      double precision dqdl,dqdr,dadp,dadq,dadl
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
c
      pi=acos(-1d0)
      Ax = (-0.75d0)*(3d0/pi)**F13
      Pconst = (3.d0*pi**2)**F23
      afact2=1d0/facttwo
c
      if (rho.ge.tol_rho) then

            rhoval=rho*facttwo
            rho43 = rhoval**F43  ! rho^4/3
            rrho  = 1d0/rhoval   ! reciprocal of rho
            rho13 = rho43*rrho 
            rho23 = rhoval**F23
            rho53 = rho23*rhoval
            rho83 = rhoval**F83
      
      
            g2 = delrho(1)*delrho(1) +
     &           delrho(2)*delrho(2) +
     &           delrho(3)*delrho(3)

            g2 = g2*facttwo*facttwo

            lapval = laprho*facttwo

            tauW = F18*g2*rrho
            tauU = 0.3d0*Pconst*rho53
c
c     Kinetic Energy Density
c
            p = g2/(4d0*Pconst*rho83)
            q = lapval/(4d0*Pconst*rho53)

            ge2 = 5.d0/27.d0*p + 20.d0/9.d0*q
            dge2dp = 5.d0/27.d0
            dge2dq = 20.d0/9.d0

            ge4 = q*(8.d0/81.d0*q - 1.d0/9.d0*p) + 8.d0/243.d0*p*p
            dge4dp = 16.d0/243.d0*p - q/9.d0
            dge4dq = 16.d0/81.d0*q - p/9.d0

            ftw = F53*p
            factor = ge4/(1d0 + ftw)
            if (factor.gt.1.d6) then
              mge4 = (1d0/ge4 + ge2/ge4 + 1d0)*(1d0 + ftw)
              dmge4dp = (dge2dp+dge4dp)/factor -
     &                  mge4*(dge4dp/ge4 - F53/(1d0 + ftw))
              dmge4dq = (dge2dq+dge4dq)/factor - mge4*dge4dq/ge4
            else
              mge4den = dsqrt(1.d0 + factor**2)
              mge4 = (1d0 + ge2 + ge4)/mge4den
              dmge4dp = (dge2dp+dge4dp)/mge4den -
     &                  mge4*ge4/(mge4den*(1d0+ftw))**2*
     &                 (dge4dp - F53*factor)
              dmge4dq = (dge2dq+dge4dq)/mge4den - 
     &                  mge4*ge4*dge4dq/(mge4den*(1d0+ftw))**2
            endif
c
            zpc = mge4 - ftw
            if (zpc.le.0.001d0) then
              fpc=0d0
              dfpcdz=0d0
            else if (zpc.lt.amgga) then
              exp8 = dexp(-amgga/zpc)
              exp9 = dexp(-amgga/(amgga-zpc))
              bfpc = (exp8*exp9 + exp8)/(exp8 + exp9)
              fpc = bfpc**bmgga
              dfpc = amgga*bmgga*fpc/(exp8 + exp9)
              dfpcdz = dfpc*exp9/zpc**2  
              if (amgga-zpc.gt.0.0013d0) then
                dfpcdz = dfpcdz + 
     &          dfpc*exp9*(1d0-exp8)/((1d0+exp9)*(amgga-zpc)**2)
              end if
            else
              fpc=1d0
              dfpcdz=0d0
            endif
            ft=ftw+zpc*fpc
            tau=tau+tauU*ft*afact2

            dftdz = fpc + zpc*dfpcdz
            dftdp = F53 + dftdz*(dmge4dp - F53)
            dftdq = dftdz*dmge4dq
            
            dpdr = -F83*p*rrho
            dpdg = 1.d0/(4d0*Pconst*rho83)

            dqdr = -F53*q*rrho
            dqdl = 1.d0/(4d0*Pconst*rho53)

            dtauUdr = F53*tauU*rrho
            dtdn = (ft*dtauUdr + tauU*(dftdp*dpdr + dftdq*dqdr))
            dtdg = 2d0*tauU*(dftdp*dpdg)
            dtdl = tauU*(dftdq*dqdl)

      else
c
        dtdn = 0d0
        dtdg = 0d0
        dtdl = 0d0
c
      end if
c
      return
      end

