      subroutine dftplusu_hamiltonian(g_dens,g_fock)
      implicit none

      integer g_dens(2),g_fock(2)
      double precision energy
#include "cdft.fh"
#include "global.fh"
#include "dftplusu.fh"
#include "errquit.fh"
#include "mafdecls.fh"

      double precision occ1,occ2,factor
      integer g_tmp1,g_tmp2,isp,ival

      if(.not.ga_create(mt_dbl,nbf_ao,nbf_ao,"noccup",0,0,g_tmp1))
     $    call errquit("dftplusu_hamiltonian",0,GA_ERR)   
      if(.not.ga_create(mt_dbl,nbf_ao,nbf_ao,"noccup",0,0,g_tmp2))
     $    call errquit("dftplusu_hamiltonian",0,GA_ERR)   
      
      energy = 0d0
      factor = 0.5d0*ipol

      do isp=1,ipol
        do ival=1,nvals
          occ1 = ga_ddot(g_dens(isp),g_spps(ival))
          !write(*,*) occ1
          call ga_dgemm("n","n",nbf_ao,nbf_ao,nbf_ao,1d0,g_dens(isp),
     $                   g_spps(ival),0d0,g_tmp1)
          call ga_dgemm("n","n",nbf_ao,nbf_ao,nbf_ao,1d0,g_spps(ival),
     $                   g_tmp1,0d0,g_tmp2)
          occ2 = factor*ga_ddot(g_dens,g_tmp2)
          !write(*,*) occ2
          energy = energy + 0.5*uvals(ival)*(occ1 - occ2)
          call ga_add(1d0,g_spps(ival),-2d0*factor,g_tmp2,g_tmp1)
          call ga_add(0.5d0*uvals(ival),g_tmp1,1d0,g_fock(isp),
     $                g_fock(isp))
        enddo
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
