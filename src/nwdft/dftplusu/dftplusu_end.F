      subroutine dftplusu_end(rtdb)
      implicit none
#include "global.fh"
#include "dftplusu.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "cdft.fh"
#include "util_params.fh"

      integer ival,rtdb,jval,j,isp
      double precision chi,chiks

      do ival=1,nvals
        if(.not.ga_destroy(g_spps(ival)))
     $    call errquit("dftplusu_end",0,GA_ERR)
      enddo

      if(response) then
        do isp=1,ipol
          if(.not.ga_destroy(g_upot0(isp)))
     $      call errquit("dftplusu_end",0,GA_ERR)
        enddo

        if(ga_nodeid().eq.0) then
          write(*,*) "             DFT+U Linear Response   "
          write(*,*)
          write(*,*) "  Site      Ref. Occ.    N.I. Occ.    I. Occ "


          do jval=1,nvals
            write(*,"(I5,5X,3F12.6)") jval,occs(1,jval),occs(2,jval),
     $                             occs(3,jval)
          enddo

          write(*,*)
          write(*,*) "     Response Functions (eV^-1) "
          do jval=1,nvals
            if (fdstep(jval).gt.0d0) then
              write(*,"(A,I4)") 
     $              "   I      ChiKS_IJ       Chi_IJ       J=",jval
              do ival=1,nvals
                chi=(occs(3,ival)-occs(1,ival))/fdstep(jval)/cau2ev
                chiks=(occs(2,ival)-occs(1,ival))/fdstep(jval)/cau2ev
                write(*,"(I4,3X,2F12.6)") ival,chiks,chi
              enddo
            endif
          enddo

        endif
      endif

      end
