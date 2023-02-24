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
      parameter(num_dirs=5)
      character*10 dirs(num_dirs)

      integer max_atoms
      parameter(max_atoms=10)

      integer atoms(max_atoms), shells(max_atoms)
      double precision ueff(max_atoms),fdstep

      double precision ddummy
      integer idummy

      character*255 test
      logical mullpop
      logical response

      data dirs /'ueff','atoms','lvals','pop','response'/

      uvals = 0
      avals = 0
      lvals = 0
      mullpop = .false.
      response = .false.

   10 if(.not.inp_a(test)) goto 1999

      if(.not.inp_match(num_dirs, .false., test, dirs, ind)) then
         call errquit("dftplusu_inp: unknown directive",0,INPUT_ERR)
      endif

      goto (1000,2000,3000,4000,5000) ind
      call errquit('dftplusu_inp: unknown directive',ind,INPUT_ERR)
      
c     Effective Hubbard U values
 1000 continue
      if(.not.inp_f(ddummy)) then
        if(uvals.eq.0)
     $  call errquit("dftplusu_inp: could not find Ueff",0,INPUT_ERR)   
        goto 10
      endif
      uvals = uvals + 1
      ueff(uvals) = ddummy/27.2114d0
      goto 1000

c     Atoms to apply U
 2000 continue
      if(.not.inp_i(idummy)) then
        if(avals.eq.0)
     $  call errquit("dftplusu_inp: could not find Atom",0,INPUT_ERR)   
        goto 10
      endif
      avals = avals + 1
      atoms(avals) = idummy
      goto 2000

c     L values to apply U
 3000 continue
      if(.not.inp_i(idummy)) then
        if(lvals.eq.0)
     $  call errquit("dftplusu_inp: could not find L",0,INPUT_ERR)   
        goto 10
      endif
      lvals = lvals + 1
      shells(lvals) = idummy
      goto 3000
c
c     Population method
 4000 continue
      if(.not.inp_a(test))
     $  call errquit("dftplusu_inp: could not find Pop",0,INPUT_ERR)   
      if (inp_compare(.false., test, 'mulliken')) then
        mullpop = .true.
      elseif (inp_compare(.false., test, 'loewdin')) then
        mullpop = .false.
      else
       call errquit("dftplusu_inp: Unknown population method",0,
     $                INPUT_ERR)   
      endif
      goto 10
c
c     Linear reponse determination of U
 5000 continue
      response = .true.
      goto 10
      
c     END
 1999 continue

      if (uvals.ne.avals)
     $ call errquit("dftplusu_inp: Ueff and Atoms mismatch",0,INPUT_ERR)
      if (uvals.ne.lvals)
     $ call errquit("dftplusu_inp: Ueff and Lvals mismatch",0,INPUT_ERR)
      
      if(uvals.eq.0) return

      if(.not.rtdb_put(rtdb,'dft:plusu:response',mt_log,1,response))
     $  call errquit("dftplusu_inp: rtdb_put failed",0,RTDB_ERR)
      if(.not.rtdb_put(rtdb,'dft:plusu:nvals',mt_int,1,uvals))
     $  call errquit("dftplusu_inp: rtdb_put failed",0,RTDB_ERR)
      if(.not.rtdb_put(rtdb,"dft:plusu:uvals",mt_dbl,uvals,ueff))
     $  call errquit("dftplusu_inp: rtdb_put failed",0,RTDB_ERR)
      if(.not.rtdb_put(rtdb,"dft:plusu:atoms",mt_int,uvals,atoms))
     $  call errquit("dftplusu_inp: rtdb_put failed",0,RTDB_ERR)
      if(.not.rtdb_put(rtdb,"dft:plusu:lvals",mt_int,uvals,shells))
     $  call errquit("dftplusu_inp: rtdb_put failed",0,RTDB_ERR)
      if(.not.rtdb_put(rtdb,"dft:plusu:mullpop",mt_log,1,mullpop))
     $  call errquit("dftplusu_inp: rtdb_put failed",0,RTDB_ERR)

      return
      end
