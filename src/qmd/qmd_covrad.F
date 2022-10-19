      subroutine qmd_covrad(nat,geom,covrad)
      implicit none
#include "geom.fh"
#include "errquit.fh"
#include "util_params.fh"
      integer nat
      double precision geom(3,nat), covrad(nat)

      integer iat, atn
      double precision chg, c(3)
      character(len=16) tag, ele
      character(len=2) sym
      logical status_tag


      double precision, parameter :: rad(118) = cang2au *  [
     & 0.40d0,0.37d0,  ! H,He
     & 1.30d0,0.99d0,0.84d0,0.75d0,0.71d0,0.64d0,0.60d0,0.62d0, ! Li-Ne
     & 1.60d0,1.40d0,1.24d0,1.14d0,1.09d0,1.04d0,1.00d0,1.01d0, ! Na-Ar
     & 2.00d0,1.74d0,  ! K,Ca
     & 1.59d0,1.48d0,1.44d0,1.30d0,1.29d0,  ! Sc-
     & 1.24d0,1.18d0,1.17d0,1.22d0,1.20d0,  ! -Zn
     & 1.23d0,1.20d0,1.20d0,1.18d0,1.17d0,1.16d0, ! Ga-Kr
     & 2.15d0,1.90d0,  ! Rbd0,Sr
     & 1.76d0,1.64d0,1.56d0,1.46d0,1.38d0,  ! Y-
     & 1.36d0,1.34d0,1.30d0,1.36d0,1.40d0,  ! -Cd
     & 1.42d0,1.40d0,1.40d0,1.37d0,1.36d0,1.36d0, ! In-Xe
     & 2.38d0,2.06d0,  ! Csd0,Ba
     & 1.94d0,1.84d0,1.90d0,1.88d0,1.86d0,1.85d0,1.83d0,! La-Eu
     & 1.82d0,1.81d0,1.80d0,1.79d0,1.77d0,1.77d0,1.78d0, ! Gd-Yb
     & 1.74d0,1.64d0,1.58d0,1.50d0,1.41d0, ! Lu-
     & 1.36d0,1.32d0,1.30d0,1.30d0,1.32d0, ! -Hg
     & 1.44d0,1.45d0,1.50d0,1.42d0,1.48d0,1.46d0,! Tl-Rn
     & 2.42d0,2.11d0,  ! Frd0,Ra
     & 2.01d0,1.90d0,1.84d0,1.83d0,1.80d0,1.80d0, ! Ac-Pu
     & 1.49d0,! Am
     & 1.49d0,1.51d0,1.51d0,1.48d0,1.50d0,1.56d0,1.58d0, ! Cm-No
     & 1.45d0,1.41d0,1.34d0,1.29d0,1.27d0, ! Lr-
     & 1.21d0,1.16d0,1.15d0,1.09d0,1.22d0, ! -Cn
     & 1.36d0,1.43d0,1.46d0,1.58d0,1.48d0,1.57d0 ] !Nh-Og

      do iat = 1,nat
        if (.not.geom_cent_get(geom,iat,tag,c,chg)) call errquit
     &      ('qmd_covrad: geom_cent_get failed:i',911,
     &       GEOM_ERR)
        status_tag = geom_tag_to_element(tag,sym,ele,atn)
        if ((sym.eq.'bq')) cycle
        if (.not.status_tag)call errquit
     &      ('qmd_covrad:geom_tag_to_element failed:i',911,
     &       GEOM_ERR)
        covrad(iat) = rad(atn)*1.75d0
      enddo

      end subroutine

