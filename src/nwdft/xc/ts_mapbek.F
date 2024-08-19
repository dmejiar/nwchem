      Subroutine ts_mapbek(tol_rho, rho, delrho, laprho,
     &                     dtdn, dtdg, dtdl, tau, facttwo)
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
      double precision Pconst, afact2, dpdg, rhoval, dpdr
c
      double precision ft,zpc,fpc,bfpc,q,lapval
      double precision exp8,exp9,factor,ftw,mge4den

c     functional derivatives below FFFFFFFFFFFF

      double precision dzpcdp,dzpcdq
      double precision dtauudr,dfpcdz,dfpc,dftdz,dftdp,dftdq
      double precision dqdl,dqdr
      double precision rK,rMu
      
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
            zpc = 1d0 + rMu*p/(1d0 + p*rMu/rK) + 20d0/9d0*q - ftw
            dzpcdp = rK*rK*rMu/(rK + rMu*p)**2 - F53
            dzpcdq = 20d0/9d0


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

