      subroutine dftplusu_end
      implicit none
#include "global.fh"
#include "dftplusu.fh"
#include "errquit.fh"

      integer ival

      do ival=1,nvals
        if(.not.ga_destroy(g_spps(ival)))
     $    call errquit("dftplusu_end",0,GA_ERR)
      enddo

      end
