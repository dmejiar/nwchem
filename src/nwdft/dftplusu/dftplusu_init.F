      subroutine dftplusu_init(rtdb,g_s,g_sm12,g_sp12)
      implicit none
      integer rtdb,g_s,g_sm12,g_sp12


#include "global.fh"
#include "tcgmsg.fh"
#include "bas.fh"
#include "geom.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "cdft.fh"
#include "dftplusu.fh"

      integer g_ppt,g_tmp1,isp
      integer ival,iat,ish,iacf,iacl
      integer ifirst,ilast,ityp,ngen,sphcart,nprim
      integer basis,nshl,ibf
      integer ilo,ihi,jlo,jhi
      logical status
      double precision ddummy

      plusuiter = 0

      if(.not.rtdb_get(rtdb,"dft:plusu:response",mt_log,1,response))
     $  call errquit("dftplusu_init: rtdb_get failed",0,RTDB_ERR)
      if(.not.rtdb_get(rtdb,"dft:plusu:nvals",mt_int,1,nvals))
     $  call errquit("dftplusu_init: rtdb_get failed",0,RTDB_ERR)
      if(.not.rtdb_get(rtdb,"dft:plusu:mullpop",mt_log,1,mullpop))
     $  call errquit("dftplusu_init: rtdb_get failed",0,RTDB_ERR)
      if(.not.rtdb_get(rtdb,"dft:plusu:lvals",mt_int,nvals,lvals))
     $  call errquit("dftplusu_init: rtdb_get failed",0,RTDB_ERR)
      if(.not.rtdb_get(rtdb,"dft:plusu:atoms",mt_int,nvals,atoms))
     $  call errquit("dftplusu_init: rtdb_get failed",0,RTDB_ERR)
      if(.not.rtdb_get(rtdb,"dft:plusu:uvals",mt_dbl,nvals,uvals))
     $  call errquit("dftplusu_init: rtdb_get failed",0,RTDB_ERR)
      if(.not.rtdb_get(rtdb,"dft:plusu:fdstep",mt_dbl,nvals,fdstep))
     $  call errquit("dftplusu_init: rtdb_get failed",0,RTDB_ERR)

      if(.not.ga_create(mt_dbl,nbf_ao,nbf_ao,"ppt",0,0,g_ppt))
     $  call errquit("dftplus_init: ga_create failed",0,GA_ERR)
      call ga_distribution(g_ppt,ga_nodeid(),ilo,ihi,jlo,jhi)

      if (.not.mullpop) then
        if(.not.ga_create(mt_dbl,nbf_ao,nbf_ao,"tmp1",0,0,g_tmp1))
     $    call errquit("dftplus_init: ga_create failed",0,GA_ERR)
          call ga_dgemm('n','n',nbf_ao,nbf_ao,nbf_ao,1d0,g_sm12,g_s,
     $                   0d0,g_sp12)   
      endif

      if (response) then
        do isp=1,ipol
          if(.not.ga_create(mt_dbl,nbf_ao,nbf_ao,"Upot0",0,0,
     $                      g_upot0(isp)))
     $    call errquit("dftplus_init: ga_create failed",0,GA_ERR)
          call ga_zero(g_upot0(isp))
        enddo
      endif


      basis = ao_bas_han

      do ival=1,nvals
        iat = atoms(ival)
        if(.not.ga_create(mt_dbl,nbf_ao,nbf_ao,"spps",0,0,g_spps(ival)))
     $    call errquit("dftplusu_init: ga_create failed",0,GA_ERR)   
        call ga_zero(g_ppt)

        status = bas_ce2cnr(basis,iat,iacf,iacl)

        nshl = 0
        do ish=iacf,iacl
          status = bas_continfo(basis,ish,ityp,nprim,ngen,sphcart)
          if (ityp.ne.lvals(ival)) cycle
          nshl = nshl + 1
          status = bas_cn2bfr(basis,ish,ifirst,ilast)
          do ibf=ifirst,ilast
            if ((ibf.ge.ilo).and.(ibf.le.ihi).and.(ibf.ge.jlo).and.
     $          (ibf.le.jhi)) then
              call ga_put(g_ppt,ibf,ibf,ibf,ibf,1d0,1)
            endif
          enddo
        enddo
        call ga_sync()

        if (mullpop) then
          call ga_dgemm('n','n',nbf_ao,nbf_ao,nbf_ao,1d0,g_s,g_ppt,
     $                   0d0,g_spps(ival))
        else
          call ga_dgemm('n','n',nbf_ao,nbf_ao,nbf_ao,1d0,g_sp12,g_ppt,
     $                   0d0,g_tmp1)
          call ga_dgemm('n','n',nbf_ao,nbf_ao,nbf_ao,1d0,g_tmp1,g_sp12,
     $                   0d0,g_spps(ival))
        endif
        
      enddo

      status = ga_destroy(g_ppt)
      if (.not.mullpop) status = status.and.ga_destroy(g_tmp1)

      if (.not.status)
     $  call errquit("dftplusu_init: ga_destroy failed",0,GA_ERR) 

      return
      end
      
