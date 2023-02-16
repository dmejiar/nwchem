      subroutine dftplusu_inp(rtdb,module)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "cdft.fh"
#include "stdio.fh"
      integer rtdb
      character*(*) module

c     Parse DFT+U input
c
      integer num_dirs, ind, uvals, avals, lvals
      parameter(num_dirs=3)
      character*10 dirs(num_dirs)

      integer max_atoms
      parameter(max_atoms=10)

      integer atoms(max_atoms), shells(max_atoms)
      double precision ueff(max_atoms)

      double precision ddummy
      integer idummy

      character*255 test

      data dirs /'ueff','atoms','lvals'/

      uvals = 0
      avals = 0
      lvals = 0

   10 if(.not.inp_a(test)) goto 1999

      if(.not.inp_match(num_dirs, .false., test, dirs, ind)) then
         call errquit("dftplusu_inp: unknown directive",0,INPUT_ERR)
      endif

      goto (1000,2000,3000) ind
      call errquit('dftplusu_inp: unknown directive',ind,INPUT_ERR)
      
c     Effective Hubbard U values
 1000 continue
      if(.not.inp_f(ddummy)) goto 10
      uvals = uvals + 1
      ueff(uvals) = ddummy/27.2114d0
      goto 1000

c     Atoms to apply U
 2000 continue
      if(.not.inp_i(idummy)) goto 10
      avals = avals + 1
      atoms(avals) = idummy
      goto 2000

c     L values to apply U
 3000 continue
      if(.not.inp_i(idummy)) goto 10
      lvals = lvals + 1
      shells(lvals) = idummy
      goto 3000

c     END
 1999 continue

      if (uvals.ne.avals)
     $ call errquit("dftplusu_inp: Ueff and Atoms mismatch",0,INPUT_ERR)
      if (uvals.ne.lvals)
     $ call errquit("dftplusu_inp: Ueff and Lvals mismatch",0,INPUT_ERR)
      
      if(uvals.eq.0) return

      if(.not.rtdb_put(rtdb,'dft:plusu:nvals',mt_int,1,uvals))
     $  call errquit("dftplusu_inp: rtdb_put failed",0,RTDB_ERR)
      if(.not.rtdb_put(rtdb,"dft:plusu:uvals",mt_dbl,uvals,ueff))
     $  call errquit("dftplusu_inp: rtdb_put failed",0,RTDB_ERR)
      if(.not.rtdb_put(rtdb,"dft:plusu:atoms",mt_int,uvals,atoms))
     $  call errquit("dftplusu_inp: rtdb_put failed",0,RTDB_ERR)
      if(.not.rtdb_put(rtdb,"dft:plusu:lvals",mt_int,uvals,shells))
     $  call errquit("dftplusu_inp: rtdb_put failed",0,RTDB_ERR)

      return
      end
