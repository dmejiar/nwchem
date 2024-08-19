      Subroutine ts_tanh(tol_rho, rho, delrho, laprho,
     &                 dtdn, dtdg, dtdl, tau, facttwo)
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
      double precision Pconst, afact2, dpdr, dpdg, rhoval
c
      double precision ft,zpc,fpc,bfpc,q,lapval
      double precision exp8,exp9,factor,ftw,mge4den

c     functional derivatives below FFFFFFFFFFFF

      double precision dzpcdp,dzpcdq
      double precision dtauudr,dfpcdz,dfpc,dftdz,dftdp,dftdq
      double precision dqdl,dqdr
      double precision rK,rMu
      double precision rz8,rz17
      integer irz
c     functional derivatives above FFFFFFFFFFFF
      
      parameter (rK = 0.804, rMu=0.23889)
      parameter (F43=4.d0/3.d0, F13=1.d0/3.d0)
      parameter (F83=8.d0/3.d0, F23=2.d0/3.d0)
      parameter (F18=1.d0/8.d0, F53=5.d0/3.d0)
      parameter (F5=5.d0, F8=8.d0)
c
      pi=acos(-1d0)
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

            ftw = F53*p
            zpc = adeorb*p + bdeorb*q - ftw
            dzpcdp = adeorb - F53
            dzpcdq = bdeorb

            if (zpc.gt.-0.55d0) then
              fpc=1d0
              dfpcdz=0d0
            else if (zpc.lt.-1d50) then
              fpc=-1d0/zpc
              dfpcdz=0d0
            else if (zpc.lt.-100d0) then
              rz=1d0/zpc
              fpc=-rz
              dfpcdz=fpc*fpc
            else if (zpc.lt.-30d0) then
              rz=1d0/zpc
              rz17=rz
              do irz=2,17
                rz17=rz17*rz
              enddo
              fpc = -rz + 1d0/24d0*rz17
              dfpcdz = rz*rz - 17d0/24d0*rz17*rz
            else
              rz=1d0/zpc
              rz8=rz
              do irz=2,8
                rz8=rz8*rz
              enddo
              fpc = dtanh(rz8)**0.125d0
              dfpcdz = -rz8/(zpc*dcosh(rz8)**2*dtanh(rz8)**0.875d0)
            endif

            ft=1d0+ftw+zpc*fpc
            tau=tau+tauU*ft*afact2

            dftdz = fpc + zpc*dfpcdz
            dftdp = F53 + dftdz*dzpcdp
            dftdq = dftdz*dzpcdq
            
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

