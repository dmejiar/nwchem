
      subroutine xc_xrevwpbe(tol_rho, fac,  rho, delrho, Amat, Cmat, nq,
     &                    ipol, Ex, qwght, ldew, func, funcid)
c
c wPBE functional based on the exchange-hole model by Ernzerhof and Perdew
c 
c M Ernzerhof and JP Perdew, JCP 109, 3313 (1998)
c J Heyd, G Scuseria and M Ernzerhof, JCP 118, 8207 (2013)
c J Heyd and G Scuseria, JCP 120, 7274 (2014)
c              


#include "dft2drv.fh"
      implicit none


      integer :: nq, ipol, funcid
      logical :: ldew
      double precision :: fac,tol_rho
      double precision :: rho(nq,*),delrho(nq,3,ipol)
      double precision :: Ex,func(nq),qwght(nq),Amat(nq,ipol)
      double precision :: Cmat(nq,*)


      if (ipol.eq.1) then
        call xc_xrevwpbe_cs(tol_rho, fac, rho, delrho, Amat,
     $                   Cmat, nq, ipol, Ex, qwght, ldew, func, 1d0,
     $                   funcid)
      else
        call xc_xrevwpbe_cs(tol_rho, fac, rho(1,2), delrho, Amat,
     $                   Cmat(1,D1_GAA), nq, ipol, Ex, qwght, ldew,
     $                   func, 2d0, funcid)
        call xc_xrevwpbe_cs(tol_rho, fac, rho(1,3), delrho(1,1,2),
     $                   Amat(1,2),
     $                   Cmat(1,D1_GBB), nq, ipol, Ex, qwght, ldew,
     $                   func, 2d0, funcid)
      endif

      end subroutine

      subroutine xc_xrevwpbe_cs(tol_rho, fac, rho, delrho, Amat, Cmat,
     $                       nq, ipol, Ex, qwght, ldew, func, afac,
     $                       funcid)

      implicit none
#include "case.fh"

      integer :: nq, ipol, funcid
      logical :: ldew
      double precision :: tol_rho, fac, afac
      double precision :: rho(nq), delrho(nq,3), Amat(nq)
      double precision :: Cmat(nq), qwght(nq), func(nq)
      double precision :: Ex

      integer :: n

      double precision :: a,b,c,d,e
      parameter (a=1.0161144d0, b=-0.37170836d0, c=-0.077215461d0)
      parameter (d=0.57786348d0, e=-0.051955731d0)

      double precision :: a1h, a2h, a3h, a4h, a5h
      parameter (a1h=0.00979681d0, a2h=0.0410834d0, a3h=0.187440d0)
      parameter (a4h=0.00120824d0, a5h=0.0347188d0)

      double precision :: ea1, ea2, ea3, ea4, ea5, ea6, ea7, ea8
      parameter (ea1=-1.128223946706117d0,ea2=1.452736265762971d0)
      parameter (ea3=-1.243162299390327d0,ea4=0.971824836115601d0)
      parameter (ea5=-0.568861079687373d0,ea6=0.246880514820192d0)
      parameter (ea7=-0.065032363850763d0,ea8=0.008401793031216d0)

      double precision :: fc1, fc2
      parameter(fc1 = 6.4753871d0)
      parameter(fc2 = 0.47965830d0)

      double precision :: f13,f43,f83
      parameter (f13=1d0/3d0, f43=4d0/3d0,f83=8d0/3d0)

      double precision :: pi, pi2, pi12
      double precision :: eb1


      double precision :: kf, omega
      double precision :: w, w2, w3, w4, w5, w6, w7, w8

      double precision :: s, p, p2, p3
      double precision :: hnum, hden, h, dh_p
      double precision :: f, df_p
      double precision :: aux1, aux2, aux3, aux4, aux5, aux6
      double precision :: daux1_p, daux2_p, daux3_p, daux4_p, daux5_p
      double precision :: daux6_p, daux3_w, daux4_w, daux5_w, daux6_w
      double precision :: aux1_12,aux1_52,aux1_72,aux1_2,aux1_3
      double precision :: aux3_12,aux3_32,aux3_52,aux3_2,aux1_92
      double precision :: aux4_12,aux4_32,aux4_52,aux4_72,aux4_2
      double precision :: aux4_3,aux4_4
      double precision :: aux6_12,aux6_32,aux6_52,aux6_72,aux6_92
      double precision :: aux6_2, aux6_3, aux6_4, aux6_5

      double precision :: term1, term2, term3, term4, term5
      double precision :: dterm1_p,dterm2_p,dterm3_p,dterm4_p,dterm5_p
      double precision :: dterm1_w,dterm3_w,dterm4_w,dterm5_w
      double precision :: term2_denom,term3_denom,term4_denom
      double precision :: term5_denom

      double precision :: epegp,egp,gb,g,dg_p,degp_p
      double precision :: t10, dt10_p, dt10_w
      double precision :: np1, dnp1_w, np2, dnp2_w
      double precision :: t1, dt1_p, dt1_w
      double precision :: f2,f3,f4,f5,f6,f7,f8,f9
      double precision :: df2_p,df3_p,df4_p,df5_p,df6_p,df7_p,df8_p
      double precision :: df9_p,df2_w,df3_w,df4_w,df5_w,df6_w,df7_w
      double precision :: df8_w,df9_w
      double precision :: t2t9,dt2t9_p,dt2t9_w

      double precision :: fx, dfx_p, dfx_w
      double precision :: dp_n, dp_g, dw_n, dp_p

      double precision :: rhoval,rho13,rho43,gamma

      double precision :: Ax,dga_p,dgb_p,ga
      double precision :: expint, erfcx
      double precision,external :: expint1

      pi = acos(-1d0)
      pi2 = pi**2
      pi12 = sqrt(pi)
      Ax = -3d0/(4d0*pi)*(3d0*pi2)**f13

      do n=1,nq

        if (rho(n).lt.tol_rho) cycle

        rhoval = afac*rho(n)
        rho13 = rhoval**f13
        rho43 = rhoval*rho13

        gamma = delrho(n,1)**2 + delrho(n,2)**2 + delrho(n,3)**2
        gamma = afac**2 * gamma
        if (gamma.lt.tol_rho**2) cycle

        kf = rho13*(3d0*pi2)**f13

        p = gamma/(4d0*kf**2*rhoval**2)

        ! Rescaling introduced in 2004
        if (p.gt.8.3d0**2) then
          p = 8.572844d0**2 + 18.796223d0**2/p**2 - 
     $        2d0*8.572844d0*18.796223d0/p
          dp_p = -2d0*18.796223d0**2/p**3 + 
     %            2d0*8.572844d0*18.796223d0/p**2
        else
          dp_p = 1d0
        endif

        p2 = p**2
        p3 = p2*p

        hnum = p*(a1h + a2h*p)
        hden = (1d0 + a3h*p2 + a4h*p**2.5 + a5h*p3)
        h = hnum/hden
        dh_p = (a1h + 2d0*a2h*p)/hden - 
     $        (2d0*a3h*p + 2.5d0*a4h*p**1.5 + 3d0*a5h*p2)*h/hden

        ! We have renamed f to simplify further computations
        f = 1d0 + p*(fc1*h + fc2)
        df_p = (fc1*h + fc2) + p*fc1*dh_p

        aux1 = h*p + d
        aux2 = 9d0*h*p/(4d0*a)
        daux1_p = h + p*dh_p
        daux2_p = 9d0/(4d0*a) * daux1_p

        aux1_2 = aux1**2
        aux1_3 = aux1*aux1_2
        aux1_52 = aux1**2.5
        aux1_72 = aux1**3.5
        aux1_92 = aux1**4.5

        if (p.gt.0.0001d0) then
          erfcx = exp(aux2)*erfc(sqrt(aux2))
          egp = -1d0*(e + 0.4d0*c*f*aux1 + 4d0/15d0*b*aux1_2 +
     $      8d0/15d0*a*aux1_3 + 0.8d0*pi12*aux1_72*(1d0-sqrt(a)*erfcx))
          degp_p = -0.4d0*(c*df_p*aux1 + 
     $             daux1_p*(c*f + f43*b*aux1 + 4d0*a*aux1_2 + 
     $             pi12*aux1_52*( 7d0 - 7d0*sqrt(a)*erfcx - 
     $             2d0*aux1*sqrt(a)*(erfcx-1d0/pi12/sqrt(aux2)
     $             )*9d0/(4d0*a))))
        else
          egp = -0.02628417880d0*p - 0.07117647788d0*p2 +
     $           0.08534541323d0*p3
          degp_p = -0.02628417880d0 - 2d0*0.07117647788d0*p +
     $              3d0*0.08534541323d0*p2     
        endif
        epegp = e + egp

        term2 = 0.5d0*(b/aux1 + c*f/aux1_2 + 2d0*epegp/aux1_3)
        dterm2_p = -daux1_p*(0.5d0*b + c*f/aux1 + 
     $              3d0*epegp/aux1_2)/aux1_2 + 0.5d0*c*df_p/aux1_2 +
     $               degp_p/aux1_3

        if (funcid.eq.1) then
          ! This the "full" PBE exchange-hole model
          w = 0d0

          expint = expint1(aux2)
          t1 = 0.5d0*a*expint
          dt1_p = 0.5d0*a*(expint - 1d0/aux2)*daux2_p

          t10 = 0.5d0*a*(log(h) + log(p) - log(aux1))
          dt10_p = 0.5d0*a*(dh_p/h + 1d0/p - daux1_p/aux1)

          fx = -8d0/9d0*(t1 + t10 + term2)
          dfx_p = -8d0/9d0*(dt1_p + dt10_p + dterm2_p)*dp_p
          dfx_w = 0d0
        else
          ! This is the short-range PBE exchange-hole model
          w = cam_omega/kf

          if (w.lt.14d0) then
            eb1 = 1.455915450052607d0
          else
            eb1 = 2d0
          endif

          w2 = w*w
          w3 = w2*w
          w4 = w3*w
          w5 = w4*w
          w6 = w5*w
          w7 = w6*w
          w8 = w7*w

          aux3 = aux1 + w2
          aux4 = h*p + eb1*w2
          aux5 = 9d0*aux4/(4d0*a)
          aux6 = d + aux4

          aux1_12 = sqrt(aux1)
          aux1_52 = aux1_2*aux1_12

          aux3_2 = aux3*aux3
          aux3_12 = sqrt(aux3)
          aux3_32 = aux3*aux3_12
          aux3_52 = aux3_2*aux3_12

          aux4_12 = sqrt(aux4)
          aux4_32 = aux4*aux4_12
          aux4_52 = aux4*aux4_32
          aux4_72 = aux4*aux4_52
          aux4_2 = aux4*aux4
          aux4_3 = aux4_2*aux4
          aux4_4 = aux4_3*aux4

          aux6_12 = sqrt(aux6)
          aux6_32 = aux6*aux6_12
          aux6_52 = aux6*aux6_32
          aux6_72 = aux6*aux6_52
          aux6_92 = aux6*aux6_72
          aux6_2 = aux6*aux6
          aux6_3 = aux6*aux6_2
          aux6_4 = aux6*aux6_3
          aux6_5 = aux6*aux6_4

          daux3_p = daux1_p
          daux4_p = daux1_p
          daux5_p = 9d0*daux1_p/(4d0*a)
          daux6_p = daux1_p

          daux3_w = 2d0*w
          daux4_w = 2d0*eb1*w
          daux5_w = 9d0*daux4_w/(4d0*a)
          daux6_w = daux4_w

          term3_denom = 8d0*aux1*aux3_52
          term3 = -w*(4d0*aux3_2*b + 6d0*aux3*c*f + 
     $            15d0*epegp)/term3_denom
          dterm3_p = -w*(8d0*aux3*daux3_p*b + 6d0*daux3_p*c*f +
     $               6d0*aux3*c*df_p + 15d0*degp_p)/
     $               term3_denom -
     $           term3/term3_denom*(8d0*daux1_p*aux3_52 +
     $                              20d0*aux1*aux3_32*daux3_p) 
          dterm3_w = -(4d0*aux3_2*b + 6d0*aux3*c*f +
     $            15d0*epegp)/term3_denom - 
     $            w*(8d0*aux3*daux3_w*b + 6d0*daux3_w*c*f)/
     $            term3_denom -
     $            term3/term3_denom*(20d0*aux1*aux3_32*daux3_w)

          term4_denom = 2d0*aux1_2*aux3_52
          term4 = -w3*(aux3*c*f + 5d0*epegp)/term4_denom
          dterm4_p = -w3*(daux3_p*c*f+aux3*c*df_p+
     $                5d0*degp_p)/term4_denom-
     $            term4/term4_denom*(4d0*aux1*daux1_p*aux3_52 +
     $                                5d0*aux1_2*aux3_32*daux3_p) 
          dterm4_w = -3d0*w2*(aux3*c*f + 5d0*epegp)/term4_denom -
     $            w3*(daux3_w*c*f)/term4_denom -
     $            term4/term4_denom*(5d0*aux1_2*aux3_32*daux3_w) 

          term5_denom = aux1_3*aux3_52
          term5 = -w5*epegp/term5_denom
          dterm5_p = -w5*degp_p/term5_denom - 
     $            term5/term5_denom*(3d0*aux1_2*daux1_p*aux3_52 +
     $                               2.5d0*aux1_3*aux3_32*daux3_p)
          dterm5_w = -5d0*w4*epegp/term5_denom -
     $           term5/term5_denom*(2.5d0*aux1_3*aux3_32*daux3_w) 

          t10 = 0.5d0*a*log(aux4/aux6)
          dt10_p = 0.5d0*a/aux4*(daux4_p - aux4*daux6_p/aux6)
          dt10_w = 0.5d0*a/aux4*(daux4_w - aux4*daux6_w/aux6)

          np1 = -1.5d0*ea1*sqrt(a)*w + 27d0*ea3*w3/(8d0*sqrt(a)) -
     $       243d0*ea5*w5/(32d0*a*sqrt(a)) +
     $       2187d0*ea7*w7/(128d0*a*a*sqrt(a))
          dnp1_w = -1.5d0*ea1*sqrt(a) + 81d0*ea3*w2/(8d0*sqrt(a)) -
     $          5d0*243d0*ea5*w4/(32d0*a*sqrt(a)) +
     $           2187d0*7d0*ea7*w6/(128d0*sqrt(a)*a**2) 

          np2 = -a + 9d0*ea2*w2/4d0 - 81d0*ea4*w4/(16d0*a) +
     $       729d0*ea6*w6/(64d0*a**2) - 6561d0*ea8*w8/(256d0*a**3)
          dnp2_w = 9d0/2d0*ea2*w - 81d0*ea4*w3/(4d0*a) +
     $         6d0*729d0*ea6*w5/(64d0*a**2) - 
     $         6561d0*ea8*w7/(32d0*a**3) 

          erfcx = exp(aux5)*erfc(sqrt(aux5))
          expint = expint1(aux5)
  
          t1 = 0.5d0*(np1*pi*erfcx - np2*expint)
          dt1_p = 0.5d0*daux5_p*(np1*(pi*erfcx - pi12/sqrt(aux5)) - 
     $                       np2*(expint - 1d0/aux5))
          dt1_w = 0.5d0*(dnp1_w*pi*erfcx -dnp2_w*expint + 
     $               daux5_w*(np1*pi*erfcx - np1*pi12/sqrt(aux5) -
     $                        np2*(expint - 1d0/aux5)))

          f2 = 0.5d0*ea1*pi12*a/sqrt(aux6)
          df2_p = -0.5d0*f2/aux6*daux6_p
          df2_w = -0.5d0*f2/aux6*daux6_w

          f3 = 0.5d0*ea2*a/aux6
          df3_p = -f3/aux6*daux6_p
          df3_w = -f3/aux6*daux6_w

          f4 = ea3*pi12*(-9d0/(8d0*aux4_12) + 0.25d0*a/aux6_32)
          df4_p = ea3*pi12*(9d0/(16d0*aux4_32)*daux4_p - 
     $        3d0/8d0*a/aux6_52*daux6_p)
          df4_w = ea3*pi12*(9d0/(16d0*aux4_32)*daux4_w - 
     $        3d0/8d0*a/aux6_52*daux6_w)

          f5 = ea4/128d0*(-144d0/aux4 + 64d0*a/aux6_2)
          df5_p = ea4/128d0*(144d0/aux4_2*daux4_p-
     $                       128d0*a/aux6_3*daux6_p)
          df5_w = ea4/128d0*(144d0/aux4_2*daux4_w-
     $                       128d0*a/aux6_3*daux6_w)

          f6 = 3d0*ea5*pi12/(32d0*a) * (27d0/aux4_12 - 6d0*a/aux4_32 +
     $                              4d0*a**2/aux6_52) 
          df6_p = 3d0*ea5*pi12/(32d0*a) * (-27d0/2d0/aux4_32*daux4_p + 
     $                                   9d0*a/aux4_52*daux4_p -
     $                                  10d0*a**2/aux6_72*daux6_p)
          df6_w = 3d0*ea5*pi12/(32d0*a) * (-27d0/2d0/aux4_32*daux4_w + 
     $                                   9d0*a/aux4_52*daux4_w -
     $                                  10d0*a**2/aux6_72*daux6_w)

          f7 = ea6*(32d0*a/aux6_3 + (-36d0 + 81d0*p*h/a)/aux4_2)/32d0
          df7_p = ea6*(-96d0*a/aux6_4*daux6_p + 81d0*h/a/aux4_2 +
     $              81d0*p*dh_p/a/aux4_2 - 
     $              2d0*(-36d0 + 81d0*p*h/a)/aux4_3*daux4_p)/32d0
          df7_w = ea6*(-96d0*a/aux6_4*daux6_w -
     $              2d0*(-36d0 + 81d0*p*h/a)/aux4_3*daux4_w)/32d0

          f8 = 3d0*pi12*ea7/(128d0*a**2) * (40d0*a**3/aux6_72 - 
     $         9d0*(4d0*a**2/aux4_52 - 6d0*a/aux4_32 + 27d0/aux4_12)) 
          df8_p =  3d0*pi12*ea7/(128d0*a**2)*
     $            (-140d0*a**3/aux6_92*daux6_p-
     $         9d0*daux4_p*(-10d0*a**2/aux4_72 + 9d0*a/aux4_52 -
     $          27d0/2d0/aux4_32))
          df8_w = 3d0*pi12*ea7/(128d0*a**2)*
     $            (-140d0*a**3/aux6_92*daux6_w-
     $         9d0*daux4_w*(-10d0*a**2/aux4_72 + 9d0*a/aux4_52 -
     $          27d0/2d0/aux4_32)) 

          f9 = 3d0*a*ea8/aux6_4 - 9d0/4d0*ea8/aux4_3 + 
     $       81d0*(ea8+ea6*eb1)/(32d0*a*aux4_2) -
     $       729d0*ea8/(128d0*a**2*aux4)   
          df9_p = 27d0*ea8/(4d0*aux4_4)*daux4_p - 
     $          81d0*(ea8 + ea6*eb1)/(16d0*a*aux4_3)*daux4_p +
     $          729d0*ea8/(128d0*a**2*aux4_2)*daux4_p - 
     $          12d0*a*ea8/aux6_5*daux6_p   
          df9_w = 27d0*ea8/(4d0*aux4_4)*daux4_w - 
     $          81d0*(ea8 + ea6*eb1)/(16d0*a*aux4_3)*daux4_w +
     $          729d0*ea8/(128d0*a**2*aux4_2)*daux4_w - 
     $          12d0*a*ea8/aux6_5*daux6_w   

          t2t9 = f2*w+f3*w2+f4*w3+f5*w4+f6*w5+f7*w6+f8*w7+f9*w8 
          dt2t9_p = df2_p*w + df3_p*w2 + df4_p*w3 + df5_p*w4 + 
     $              df6_p*w5 + df7_p*w6 + df8_p*w7 + df9_p*w8
          dt2t9_w= df2_w*w + df3_w*w2 + df4_w*w3 + df5_w*w4 + df6_w*w5 +
     $          df7_w*w6 + df8_w*w7 + df9_w*w8 +
     $          f2 + 2d0*f3*w + 3d0*f4*w2 + 4d0*f5*w3 + 5d0*f6*w4 +
     $          6d0*f7*w5 + 7d0*f8*w6 + 8d0*f9*w7

          if (w.gt.14d0) then
            term1 = -0.5d0*a*(log(aux6) - log(aux4) - expint)
            dterm1_p = -0.5d0*a*(1d0/aux6*daux6_p - 1d0/aux4*daux4_p -
     $                  (expint - 1d0/aux5)*daux5_p)
            dterm1_w = -0.5d0*a*(1d0/aux6*daux6_w - 1d0/aux4*daux4_w -
     $                  (expint - 1d0/aux5)*daux5_w)
          else
            term1 = t1 + t2t9 + t10
            dterm1_p = dt1_p + dt2t9_p + dt10_p
            dterm1_w = dt1_w + dt2t9_w + dt10_w
          endif

          fx = -8d0/9d0*(term1 + term2 + term3 + term4 + term5)
          dfx_p = -8d0/9d0*(dterm1_p + dterm2_p + dterm3_p + dterm4_p +
     $                  dterm5_p)*dp_p
          dfx_w = -8d0/9d0*(dterm1_w + dterm3_w + dterm4_w + dterm5_w) 
        endif

        dp_g = 1d0/(4d0*(3d0*pi2)**(2d0/3d0)*rhoval**(8d0/3d0))

        Ex = Ex + fac*Ax*rho43*fx*qwght(n)/afac
        if (ldew) func(n) = func(n) + fac*Ax*rho43*fx/afac
    
        Amat(n) = Amat(n) + fac*(f43*Ax*fx*rho13-f83*Ax*rho13*dfx_p*p-
     $                      f13*Ax*rho13*dfx_w*w)
        Cmat(n) = Cmat(n) + 2d0*fac*(Ax*rho43*dfx_p*dp_g)
      enddo

      end subroutine

      double precision function expint1(x)
      implicit none
      double precision   :: x
      integer, parameter :: maxit = 200
      double precision, parameter :: eps = epsilon(1d0)
      double precision, parameter :: fpmin = tiny(1d0)*1d2
      double precision, parameter :: euler = 0.57721566490153286d0

      integer :: i, ii
      double precision :: a,b,c,d,del,fact,h,psi,y

      if (x.eq.0d0) then
        y = epsilon(x)
      else
        y = x
      endif

      if (y.gt.1d0) then
              b = y + 1d0
              c = 1d0/fpmin
              d = 1d0/b
              h = d
              do i=1,maxit
                a = -i**2
                b = b + 2d0
                d = 1d0/(a*d + b)
                c = b + a/c
                del = c*d
                h = h*del
                if (abs(del-1d0).lt.eps) then
                        expint1 = h
                        return
                endif
              enddo
      else
              expint1 = -(log(y) + euler)
              fact = 1d0
              do i=1,maxit
                fact = -fact*y/dble(i)
                del = -fact/dble(i)
                expint1 = expint1 + del
                if (abs(del).lt.abs(expint1)*eps) then
                        expint1 = exp(x)*expint1
                        return
                endif
              enddo
      endif

      return
      end function

