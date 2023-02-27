      subroutine dftplusu_hamiltonian(g_dens,g_fock)
      implicit none

      integer g_dens(2),g_fock(2)
      double precision energy
#include "bas.fh"
#include "cdft.fh"
#include "global.fh"
#include "dftplusu.fh"
#include "errquit.fh"
#include "mafdecls.fh"

      double precision occ1,occ2,factor
      integer g_tmp1,g_tmp2,isp,ival
      logical status
      integer nprim,ityp,ngen,iat,iacf,iacl,iabf,iabl,ibf,ish
      integer ifirst,ilast,sphcart

      if(.not.ga_create(mt_dbl,nbf_ao,nbf_ao,"noccup",0,0,g_tmp1))
     $    call errquit("dftplusu_hamiltonian",0,GA_ERR)   
      if(.not.ga_create(mt_dbl,nbf_ao,nbf_ao,"noccup",0,0,g_tmp2))
     $    call errquit("dftplusu_hamiltonian",0,GA_ERR)   
      
      energy = 0d0
      factor = 0.5d0*ipol
      plusuiter = plusuiter + 1

      if (response) then
        if (plusuiter.le.1) occs(1,1:nvals) = 0
        if (plusuiter.le.2) occs(2,1:nvals) = 0
        occs(3,1:nvals) = 0
      endif

      do isp=1,ipol
        do ival=1,nvals
          call ga_copy(g_spps(ival),g_tmp1)
          if (mullpop) call ga_symmetrize(g_tmp1)
          occ1 = ga_ddot(g_dens(isp),g_spps(ival))
          call ga_add(0.5d0*uvals(ival),g_tmp1,1d0,g_fock(isp),
     $                g_fock(isp))

          if (response) then
            call ga_add(fdstep(ival),g_tmp1,1d0,g_fock(isp),
     $                  g_fock(isp))

            if (plusuiter.eq.1) then
              occs(1,ival) = occs(1,ival) + occ1
            elseif (plusuiter.eq.2) then
              occs(2,ival) = occs(2,ival) + occ1
            endif
            occs(3,ival) = occs(3,ival) + occ1
            energy = energy + fdstep(ival)*occ1
            
            if(plusuiter.gt.1) then
              energy = energy + 0.5d0*uvals(ival)*occ1
              cycle
            endif
          endif
          !write(*,*) occ1

          call ga_dgemm("n","n",nbf_ao,nbf_ao,nbf_ao,1d0,g_dens(isp),
     $                   g_spps(ival),0d0,g_tmp1)
          call ga_dgemm("n","n",nbf_ao,nbf_ao,nbf_ao,1d0,g_spps(ival),
     $                   g_tmp1,0d0,g_tmp2)
          if (mullpop) call ga_symmetrize(g_tmp2)
          call ga_add(-uvals(ival),g_tmp2,1d0,g_fock(isp),g_fock(isp))
          occ2 = factor*ga_ddot(g_dens(isp),g_tmp2)
          !write(*,*) occ2
          if (response) then
            call ga_add(-uvals(ival),g_tmp2,1d0,g_upot0(isp),
     $                   g_upot0(isp))
          endif
          energy = energy + 0.5d0*uvals(ival)*(occ1 - occ2)
        enddo
        if (response.and.plusuiter.gt.1) then
          call ga_add(1d0,g_upot0(isp),1d0,g_fock(isp),g_fock(isp))
          energy = energy + 0.5d0*ga_ddot(g_upot0(isp),g_dens(isp))
        endif
      enddo
      edftplusu = energy

      if(.not.ga_destroy(g_tmp1))
     $  call errquit("dftplusu_hamiltonian",0,GA_ERR) 
      if(.not.ga_destroy(g_tmp2))
     $  call errquit("dftplusu_hamiltonian",0,GA_ERR)

      end

      double precision function get_edftplusu()
      implicit none
#include "dftplusu.fh"
      get_edftplusu = edftplusu
      end
