      subroutine dftplusu_end(rtdb)
      implicit none
#include "global.fh"
#include "dftplusu.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"

      integer ival,rtdb

      do ival=1,nvals
        if(.not.ga_destroy(g_spps(ival)))
     $    call errquit("dftplusu_end",0,GA_ERR)
      enddo

      if(response) then
        if(ga_nodeid().eq.0) then
          write(*,*) "       DFT+U Linear Response"
          do ival=1,nvals
            write(*,"(A,I5)")   "                      Atom: ",
     $        atoms(ival)
            write(*,"(A,F10.4)") "      Reference Occupation: ",
     $        occs(1,ival)
            write(*,"(A,F10.4)") "Non interacting Occupation: ",
     $        occs(2,ival)
            write(*,"(A,F10.4)") "    Interacting Occupation: ",
     $        occs(3,ival)
            write(*,"(A,F10.4,A)") " Linear Response Hubbard U: ",
     $        (0.5d0*uvals(ival)/(occs(2,ival) - occs(1,ival)) -
     $         0.5d0*uvals(ival)/(occs(3,ival) -
     $                                  occs(1,ival)))*27.2114d0," eV"
            write(*,*)
          enddo
        endif
      endif

      end
