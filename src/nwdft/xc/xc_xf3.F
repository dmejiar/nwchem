c Copyright 2019 (C) Orbital-free DFT group at University of Florida
c Licensed under the Educational Community License, Version 2.0 
c (the "License"); you may not use this file except in compliance with 
c the License. You may obtain a copy of the License at
c
c    http://www.osedu.org/licenses/ECL-2.0
c
c Unless required by applicable law or agreed to in writing,
c software distributed under the License is distributed on an "AS IS"
c BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express
c or implied. See the License for the specific language governing
c permissions and limitations under the License.
c
#include "dft2drv.fh"
#include "dft3drv.fh"
c     Nearly correct asymptotic potential (NCAP) functional
c     (Exchange part only)
c           GGA
C         utilizes ingredients:
c                              rho   -  density
c                              delrho - gradient of density
c
c     Written by:
c     Daniel Mejia-Rodriguez
c     QTP, Department of Physics, University of Florida
c
c     References:
c     J. Carmona-Espindola, J.L. Gazquez, A. Vela, S.B. Trickey
c     JCTC 15, 303 (2019).
c     DOI: 10.1021/acs.jctc.8b00998
c
c
#if !defined SECOND_DERIV && !defined THIRD_DERIV
      Subroutine xc_xf3(tol_rho, fac, rho, delrho,
     &                    Amat, Cmat, nq, ipol, Ex, qwght,ldew,func,
     &                    version)
#elif defined(SECOND_DERIV) && !defined THIRD_DERIV
      Subroutine xc_xf3_d2(tol_rho, fac, rho, delrho,
     &                       Amat, Amat2, Cmat, Cmat2, nq, ipol, Ex,
     &                       qwght,ldew,func,
     &                       version)
#else
      Subroutine xc_xf3_d3(tol_rho, fac, rho, delrho,
     &                       Amat, Amat2, Amat3, Cmat, Cmat2, Cmat3, 
     &                       nq, ipol, Ex, qwght,ldew,func,
     &                       version)
#endif
c
      implicit none
c
      character*1 version
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
c     Quadrature Weights
c
      double precision qwght(nq)
c
c     Sampling Matrices for the XC Potential & Energy
c
      double precision amat(nq,ipol), cmat(nq,*)
c
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
      double precision Amat2(nq,NCOL_AMAT2), Cmat2(nq,NCOL_CMAT2)
#endif
#ifdef THIRD_DERIV
      double precision Amat3(nq,NCOL_AMAT3), Cmat3(nq,NCOL_CMAT3)
#endif
c
      double precision tol_rho, pi, mu, alpha, beta, zeta
      double precision ckf, Ax
      double precision F43, F13, F23, F49
c
      parameter(mu=0.2195149727645171d0)

      parameter (F43=4.d0/3.d0, F13=1.d0/3.d0, F23=2.0d0/3.0d0)
      parameter (F49=4d0/9d0)
c
      integer n
      double precision rrho, rho43, rho13, gamma, gam12, s, s2, d1s
      double precision arcsinh, darcsinh, d2arcsinh, d3arcsinh
      double precision Fx, dFx, d2Fx, d3Fx
      double precision tansin, dtansin
      double precision factor, dfactor, d2factor, d3factor
      double precision asymp, dasymp
      double precision denom, ddenom
      double precision small, dsmall, d2small, d3small
      double precision large, dlarge, d2large, d3large
      double precision temp, dtemp, d2temp, p, b1, b2
      double precision atanhs, sechs2, sech2s, tanh2s, tanhs

#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
      double precision d2denom, d2asymp, d2tansin, d2s, rhom23
#endif
#if defined(THIRD_DERIV)
      double precision d3denom, d3asymp, d3tansin, d3s, rhom53
#endif
      double precision ops, omz, logops

      double precision s2overgamma,s2overgamma_rho43
      double precision dFxd1s,dFxd2s,fact

c
      pi    = dacos(-1.d0)
      beta  = 0.00912345d0
      b2    = 0.34971d0
      b1    = 1.0d0 - 2*(1.0 - dlog(6.d0*pi)/3.d0)*b2

      ckf   = 2*(3d0*pi*pi)**F13
      alpha = 4d0*pi/3d0*beta/mu
      Ax    = -0.375d0*ckf/pi
c
      if (ipol.eq.1 )then
c
c        ======> SPIN-RESTRICTED <======
c
#ifdef IFCV81
CDEC$ NOSWP
#endif
         do 10 n = 1, nq
            if (rho(n,1).lt.tol_rho) goto 10
            rho43 = rho(n,1)**F43
            rrho = 1d0/rho(n,1)
            rho13 = rho43*rrho
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
            if (dsqrt(gamma).le.tol_rho**2) goto 10

            s = dsqrt(gamma)/(ckf*rho43)
            p = s*s

            if(p.lt.1d-8) then
              Fx  = 1.0d0 + mu*p + alpha*b1*mu*p*s - 
     $              (beta + F23)*p*p + 
     $              (alpha*b2 - alpha*b1*(beta + F23))*mu*p*p*s +
     $              (F23 + F43*beta + beta**2)*mu*p*p*p
              dFx  = mu + 1.5d0*alpha*b1*mu*s -
     $               2d0*(beta + F23)*p +
     $               2.5d0*(b2 - b1*(beta + F23))*alpha*mu*p*s +
     $               (2d0 + 4d0*beta + 3d0*beta**2)*mu*p*p
              d2Fx = 0.75d0*alpha*b1*mu/s - 2d0*(beta + F23) +
     $               3.75d0*(b2 - b1*(beta + F23))*alpha*mu*s +
     $               2d0*(2d0 + 4d0*beta + 3d0*beta**2)*mu*p
              d3Fx = -0.375d0*alpha*b1*mu/(p*s) + 
     $                1.875d0*(b2 - b1*(beta + F23))*alpha*mu/s +
     $                2d0*(2d0 + 4d0*beta + 3d0*beta**2)*mu
            else
              if(s.gt.19d0) then
                tanhs  = 1d0
                atanhs = 0.25d0*pi
                temp   = 0.5d0*atanhs/s
                dtemp  = -0.25d0*atanhs/(s*p)
                d2temp = 0.375d0*atanhs/(p*p*s)
              else
                tanhs  = dtanh(s)
                atanhs = datan(dtanh(s))
                sechs2 = 1.0d0 - dtanh(s)**2
                sech2s = sechs2/(1.0d0 + tanhs**2)
                tanh2s = dtanh(2*s)
                temp   = 0.5d0*(atanhs/s + sech2s)
                dtemp  = 0.25d0*sech2s/p - 0.25d0*atanhs/(s*p) -
     $                   0.5d0*sech2s*tanh2s/s
                d2temp = 0.125d0*(3d0*atanhs + 
     $              s*sech2s*(4d0*p - 8d0*p*sech2s**2 - 3d0))/(p*p*s)
              endif

              if(s.gt.1000d0) then
                fact = 0.25d0*pi*s*beta
                small = mu/beta*(1.0d0 - 1.0d0/fact + 1.0d0/fact**2 -
     $                           1.0d0/fact**3 + 1.0d0/fact**4 -
     $                           1.0d0/fact**5 + 1.0d0/fact**6)
                dsmall = mu/beta/p*(0.5d0/fact - 1.0d0/fact**2 +
     $                    1.5d0/fact**3 - 2/fact**4 + 2.5d0/fact**5 -
     $                    3.0d0/fact**6)
                d2small = mu/beta/p**2*(-0.75d0/fac + 2.0d0/fac**2 -
     $                    3.75d0/fact**3 + 6/fact**4 - 8.75d0/fact**5 +
     $                    12.0d0/fact**6)
                d3small = mu/beta/p**3*(1.875d0/fact - 6.0d0/fact**2 +
     $                    13.125d0/fact**3 - 24.0d0/fact**4 +
     $                    39.375d0/fact**5 - 60.0d0/fact**6)
              else
                small  = mu*s*atanhs/(1.0d0 + beta*s*atanhs)
                dsmall = mu*temp/(1.0d0 + beta*s*atanhs)**2
                d2small = mu*dtemp/(1.0d0 + beta*s*atanhs)**2 -
     $                    2*beta*temp*dsmall/(1d0 + beta*s*atanhs)
                d3small = mu*d2temp/(1d0 + beta*s*atanhs)**2 -
     $                  2*mu*dtemp/(1d0 + beta*s*atanhs)**3*beta*temp -
     $                  2*beta*dtemp*dsmall/(1d0 + beta*s*atanhs) -
     $                  2*beta*temp*d2small/(1d0 + beta*s*atanhs) +
     $             2*beta*temp*dsmall/(1d0 + beta*s*atanhs)**2*beta*temp
              endif

              if(s.gt.1000d0) then
                factor = b1 + b2*log(1.0d0+p)
                large = 1.0d0 + alpha*s*factor
                dlarge = alpha/s*(0.5d0*factor +
     $                   b2*(1d0 - 1d0/p + 1d0/(p*p) - 1d0/(p*p*p)))
                d2large = alpha/(s*p)*(-0.25d0*factor + 
     $                    b2*(1.0d0/p - 2d0/(p*p) + 3.0d0/(p*p*p)))
                d3large = alpha/(s*p*p)*(0.375d0*factor +
     $                    0.25d0*b2*(-1d0 - 9d0/p +27d0/(p*p)))
              else
                factor   = b1 + b2*log(1.0d0 + p)
                dfactor  = b2/(1.0d0 + p)
                d2factor = -dfactor/(1d0 + p)
                d3factor = -2d0*d2factor/(1d0 + p)

                large   = 1.0d0 + alpha*s*factor
                dlarge  = 0.5d0*alpha*factor/s + alpha*s*dfactor
                d2large = -0.25d0*alpha*factor/(s*p) + alpha*dfactor/s +
     $                    alpha*s*d2factor
                d3large = -0.75d0*alpha*dfactor/(s*p) +
     $                     0.375d0*alpha*factor/(s*p*p) +
     $                     1.5d0*alpha*d2factor/s +
     $                     alpha*s*d3factor
              endif
              
              Fx   = 1d0 + small*large
              dFx  = dsmall*large + small*dlarge
              d2Fx = d2small*large + 2d0*dsmall*dlarge + small*d2large
              d3Fx = d3small*large + 3*d2small*dlarge +
     $               3*dsmall*d2large + small*d3large
            endif

            s2overgamma = 1d0/(ckf*ckf*rho43*rho43)
            s2overgamma_rho43 = 1d0/(ckf*ckf*rho43)

            Ex = Ex + Ax*rho43*Fx*qwght(n)*fac
            if (ldew) func(n) = func(n) + Ax*rho43*Fx*fac

            Amat(n,1) = Amat(n,1) + F43*Ax*rho13*(Fx - 2d0*p*dFx)*fac
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + 
     $                          2d0*Ax*rho43*dFx*s2overgamma*fac

#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            rhom23 = rho13*rrho
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + 2d0*F49*Ax*rhom23*
     &                          (Fx + 6d0*p*dFx + 16d0*p*p*d2Fx)*fac

            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) -
     &            F43*Ax*rho13*s2overgamma*(dFx + 2d0*p*d2Fx)*4d0*fac

            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + Ax*
     &            s2overgamma_rho43*s2overgamma*d2Fx*8d0*fac
#endif
#ifdef THIRD_DERIV
            rhom53 = rhom23*rrho
            d3s = -1.5d0*d2s/gamma

            Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA) - 4d0*fac*
     &      F23*F49*Ax*rhom53*(Fx + 34d0*p*dFx + 
     $         168d0*p*p*d2Fx + 64d0*p*p*p*d3Fx)

            Cmat3(n,D3_RA_RA_GAA) = Cmat3(n,D3_RA_RA_GAA) +
     &      F49*Ax*rhom23*s2overgamma*(7d0*dFx + 38d0*p*d2Fx + 
     &                                 16d0*d3Fx*p*p)*8d0*fac

            Cmat3(n,D3_RA_GAA_GAA) = Cmat3(n,D3_RA_GAA_GAA) -
     &      F43*Ax*rho13*s2overgamma**2*(3d0*d2Fx + 2d0*p*d3Fx)*16d0*fac

            Cmat3(n,D3_GAA_GAA_GAA) = Cmat3(n,D3_GAA_GAA_GAA) +
     &                 Ax*s2overgamma_rho43*s2overgamma**2*d3Fx*32d0*fac

#endif
 10      continue
c
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
#ifdef IFCV81
CDEC$ NOSWP
#endif
         do 20 n = 1, nq
            if (rho(n,1).lt.tol_rho) goto 20
c
c     Alpha
c
            if (rho(n,2).lt.tol_rho) goto 25
            rho43 = (2d0*rho(n,2))**F43
            rrho = 0.5d0/rho(n,2)
            rho13 = rho43*rrho
c
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
            gam12 = 2d0*dsqrt(gamma)
            if (gam12.le.tol_rho**2) goto 25
c
            s = gam12/(2d0*ckf*rho43)
            p = s*s


            if(p.lt.1d-8) then
              Fx  = 1.0d0 + mu*p + alpha*b1*mu*p*s - 
     $              (beta + 2d0/3d0)*p*p + 
     $              (alpha*b2 - alpha*b1*(beta + 2d0/3d0))*mu*p*p*s +
     $              (2d0/3d0 + 4d0/3d0*beta + beta**2)*mu*p*p*p
              dFx  = mu + 1.5d0*alpha*b1*mu*s -
     $               2d0*(beta + 2d0/3d0)*p +
     $               2.5d0*(b2 - b1*(beta + 2d0/3d0))*alpha*mu*p*s +
     $               (2d0 + 4d0*beta + 3d0*beta**2)*mu*p*p
              d2Fx = 0.75d0*alpha*b1*mu/s - 2d0*(beta + 2d0/3d0) +
     $               3.75d0*(b2 - b1*(beta + 2d0/3d0))*alpha*mu*s +
     $               2d0*(2d0 + 4d0*beta + 3d0*beta**2)*mu*p
              d3Fx = -0.375d0*alpha*b1*mu/(p*s) + 
     $                1.875d0*(b2 - b1*(beta + 2d0/3d0))*alpha*mu/s +
     $                2d0*(2d0 + 4d0*beta + 3d0*beta**2)*mu
            else
              if(s.gt.19d0) then
                atanhs = 0.25d0*pi
                temp   = 0.5d0*atanhs/s
                dtemp  = -0.25d0*atanhs/(s*p)
                d2temp = 0.375d0*atanhs/(p*p*s)
              else
                atanhs = atan(tanh(s))
                sechs2 = 1.0d0 - tanh(s)**2
                sech2s = sechs2/(1.0d0 + tanhs**2)
                tanh2s = tanh(2*s)
                temp   = 0.5d0*(atanhs/s + sech2s)
                dtemp  = 0.25d0*sech2s/p - 0.25d0*atanhs/(s*p) -
     $                   0.5d0*sech2s*tanh2s/s
                d2temp = 0.125d0*(3d0*atanhs + 
     $              s*sech2s*(4d0*p - 8d0*p*sech2s**2 - 3d0))/(p*p*s)
              endif

              small  = mu*s*atanhs/(1.0d0 + beta*s*atanhs)
              dsmall = mu*temp/(1.0d0 + beta*s*atanhs)**2
              d2small = mu*dtemp/(1.0d0 + beta*s*atanhs)**2 -
     $                  2*beta*temp*dsmall/(1d0 + beta*s*atanhs)
              d3small = mu*d2temp/(1d0 + beta*s*atanhs)**2 -
     $                  2*mu*dtemp/(1d0 + beta*s*atanhs)**3*beta*temp -
     $                  2*beta*dtemp*dsmall/(1d0 + beta*s*atanhs) -
     $                  2*beta*temp*d2small/(1d0 + beta*s*atanhs) +
     $           2*beta*temp*dsmall/(1d0 + beta*s*atanhs)**2*beta*temp

              factor   = b1 + b2*log(1.0d0 + p)
              dfactor  = b2/(1.0d0 + p)
              d2factor = -dfactor/(1d0 + p)
              d3factor = -2d0*d2factor/(1d0 + p)

              large   = 1.0d0 + alpha*s*factor
              dlarge  = 0.5d0*alpha*factor/s + alpha*s*dfactor
              d2large = -0.25d0*alpha*factor/(s*p) + alpha*dfactor/s +
     $                  alpha*s*d2factor
              d3large = -0.75d0*alpha*dfactor/(s*p) +
     $                   0.375d0*alpha*factor/(s*p*p) +
     $                   1.5d0*alpha*d2factor/s +
     $                   alpha*s*d3factor
              
              Fx   = 1d0 + small*large
              dFx  = dsmall*large + small*dlarge
              d2Fx = d2small*large + 2d0*dsmall*dlarge + small*d2large
              d3Fx = d3small*large + 3*d2small*dlarge +
     $               3*dsmall*d2large + small*d3large
            endif

            s2overgamma = 4d0/(ckf*ckf*rho43*rho43)
            s2overgamma_rho43 = 4d0/(ckf*ckf*rho43)

            Ex = Ex + 0.5d0*Ax*rho43*Fx*qwght(n)*fac
            if (ldew) func(n) = func(n) + 0.5d0*Ax*rho43*Fx*fac

            Amat(n,1) = Amat(n,1) + F43*Ax*rho13*(Fx - 2*p*dFx)*fac
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) +
     &                       0.5d0*Ax*rho43*dFx*s2overgamma*fac

#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            rhom23 = rho13*rrho
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + 2d0*F49*Ax*rhom23*
     &                          (Fx + 6*p*dFx + 16*p*p*d2Fx)*fac

            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) - F43*Ax*rho13*
     &                           s2overgamma*(dFx + 2*p*d2Fx)*fac

            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + fac*Ax*rho43*
     &                            (d2Fx*s2overgamma**2)*0.5d0
#endif
#ifdef THIRD_DERIV
            rhom53 = rhom23*rrho

            Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA) - 4d0*fac*
     &      F23*F49*Ax*rhom53*(Fx - s*dFx + 18d0*s2*d2Fx +8d0*s2*s*d3Fx)

            Cmat3(n,D3_RA_RA_GAA) = Cmat3(n,D3_RA_RA_GAA) +
     &      F49*Ax*rhom23*d1s*s*(7d0*d2Fx + 4d0*d3Fx*s)*2d0*fac

            Cmat3(n,D3_RA_GAA_GAA) = Cmat3(n,D3_RA_GAA_GAA) -
     &      F43*Ax*rho13*(d2Fx*d1s**2 + d3Fx*d1s**2*s +
     &                     d2Fx*d2s*s)*fac

            Cmat3(n,D3_GAA_GAA_GAA) = Cmat3(n,D3_GAA_GAA_GAA) +
     &      Ax*rho43*(d3Fx*d1s**3 + 3d0*d2Fx*d1s*d2s +
     &                dFx*d3s)*0.5d0*fac

#endif
c
c     Beta
c
 25         continue
            if (rho(n,3).lt.tol_rho) goto 20
            rho43 = (2d0*rho(n,3))**F43
            rrho = 0.5d0/rho(n,3)
            rho13 = rho43*rrho
c
            gamma = delrho(n,1,2)*delrho(n,1,2) +
     &              delrho(n,2,2)*delrho(n,2,2) +
     &              delrho(n,3,2)*delrho(n,3,2)
            gam12 = 2d0*dsqrt(gamma)
            if (gam12.le.tol_rho**2) goto 20
c
            s = gam12/(2d0*ckf*rho43)
            p = s*s

            if(p.lt.1d-8) then
              Fx  = 1.0d0 + mu*p + alpha*b1*mu*p*s - 
     $              (beta + 2d0/3d0)*p*p + 
     $              (alpha*b2 - alpha*b1*(beta + 2d0/3d0))*mu*p*p*s +
     $              (2d0/3d0 + 4d0/3d0*beta + beta**2)*mu*p*p*p
              dFx  = mu + 1.5d0*alpha*b1*mu*s -
     $               2d0*(beta + 2d0/3d0)*p +
     $               2.5d0*(b2 - b1*(beta + 2d0/3d0))*alpha*mu*p*s +
     $               (2d0 + 4d0*beta + 3d0*beta**2)*mu*p*p
              d2Fx = 0.75d0*alpha*b1*mu/s - 2d0*(beta + 2d0/3d0) +
     $               3.75d0*(b2 - b1*(beta + 2d0/3d0))*alpha*mu*s +
     $               2d0*(2d0 + 4d0*beta + 3d0*beta**2)*mu*p
              d3Fx = -0.375d0*alpha*b1*mu/(p*s) + 
     $                1.875d0*(b2 - b1*(beta + 2d0/3d0))*alpha*mu/s +
     $                2d0*(2d0 + 4d0*beta + 3d0*beta**2)*mu
            else
              if(s.gt.19d0) then
                atanhs = 0.25d0*pi
                temp   = 0.5d0*atanhs/s
                dtemp  = -0.25d0*atanhs/(s*p)
                d2temp = 0.375d0*atanhs/(p*p*s)
              else
                atanhs = atan(tanh(s))
                sechs2 = 1.0d0 - tanh(s)**2
                sech2s = sechs2/(1.0d0 + tanhs**2)
                tanh2s = tanh(2*s)
                temp   = 0.5d0*(atanhs/s + sech2s)
                dtemp  = 0.25d0*sech2s/p - 0.25d0*atanhs/(s*p) -
     $                   0.5d0*sech2s*tanh2s/s
                d2temp = 0.125d0*(3d0*atanhs + 
     $              s*sech2s*(4d0*p - 8d0*p*sech2s**2 - 3d0))/(p*p*s)
              endif

              small  = mu*s*atanhs/(1.0d0 + beta*s*atanhs)
              dsmall = mu*temp/(1.0d0 + beta*s*atanhs)**2
              d2small = mu*dtemp/(1.0d0 + beta*s*atanhs)**2 -
     $                  2*beta*temp*dsmall/(1d0 + beta*s*atanhs)
              d3small = mu*d2temp/(1d0 + beta*s*atanhs)**2 -
     $                  2*mu*dtemp/(1d0 + beta*s*atanhs)**3*beta*temp -
     $                  2*beta*dtemp*dsmall/(1d0 + beta*s*atanhs) -
     $                  2*beta*temp*d2small/(1d0 + beta*s*atanhs) +
     $           2*beta*temp*dsmall/(1d0 + beta*s*atanhs)**2*beta*temp

              factor   = b1 + b2*log(1.0d0 + p)
              dfactor  = b2/(1.0d0 + p)
              d2factor = -dfactor/(1d0 + p)
              d3factor = -2d0*d2factor/(1d0 + p)

              large   = 1.0d0 + alpha*s*factor
              dlarge  = 0.5d0*alpha*factor/s + alpha*s*dfactor
              d2large = -0.25d0*alpha*factor/(s*p) + alpha*dfactor/s +
     $                  alpha*s*d2factor
              d3large = -0.75d0*alpha*dfactor/(s*p) +
     $                   0.375d0*alpha*factor/(s*p*p) +
     $                   1.5d0*alpha*d2factor/s +
     $                   alpha*s*d3factor
              
              Fx   = 1d0 + small*large
              dFx  = dsmall*large + small*dlarge
              d2Fx = d2small*large + 2d0*dsmall*dlarge + small*d2large
              d3Fx = d3small*large + 3*d2small*dlarge +
     $               3*dsmall*d2large + small*d3large
            endif

            s2overgamma = 4d0/(ckf*ckf*rho43*rho43)
            s2overgamma_rho43 = 4d0/(ckf*ckf*rho43)

            Ex = Ex + 0.5d0*Ax*rho43*Fx*qwght(n)*fac
            if (ldew) func(n) = func(n) + 0.5d0*Ax*rho43*Fx*fac

            Amat(n,2) = Amat(n,2) + F43*Ax*rho13*(Fx - 2*p*dFx)*fac
            Cmat(n,D1_GBB) = Cmat(n,D1_GBB) +
     $                       0.5d0*Ax*rho43*dFx*s2overgamma*fac

#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            d2s = -0.5d0*d1s/gamma
            rhom23 = rho13*rrho
            Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + 2d0*F49*Ax*rhom23*
     &                          (Fx + 6*p*dFx + 16*p*p*d2Fx)*fac

            Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB) - F43*Ax*rho13*
     &                           s2overgamma*(dFx + 2*p*d2Fx)*fac

            Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB) + fac*Ax*rho43*
     &                s2overgamma**2*(d2Fx)*0.5d0
#endif
#ifdef THIRD_DERIV
            rhom53 = rhom23*rrho
            d3s = -1.5d0*d2s/gamma

            Amat3(n,D3_RB_RB_RB) = Amat3(n,D3_RB_RB_RB) - 4d0*fac*
     &      F23*F49*Ax*rhom53*(Fx - s*dFx + 18d0*s2*d2Fx +8d0*s2*s*d3Fx)

            Cmat3(n,D3_RB_RB_GBB) = Cmat3(n,D3_RB_RB_GBB) +
     &      F49*Ax*rhom23*d1s*s*(7d0*d2Fx + 4d0*d3Fx*s)*2d0*fac

            Cmat3(n,D3_RB_GBB_GBB) = Cmat3(n,D3_RB_GBB_GBB) -
     &      F43*Ax*rho13*(d2Fx*d1s**2 + d3Fx*d1s**2*s +
     &                     d2Fx*d2s*s)*fac

            Cmat3(n,D3_GBB_GBB_GBB) = Cmat3(n,D3_GBB_GBB_GBB) +
     &      Ax*rho43*(d3Fx*d1s**3 + 3d0*d2Fx*d1s*d2s +
     &                dFx*d3s)*0.5d0*fac
#endif
c
 20      continue
      endif
c
      return
      end
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "xc_xf3.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "xc_xf3.F"
#endif
c $Id$
