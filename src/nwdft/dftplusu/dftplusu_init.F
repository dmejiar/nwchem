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

      integer size1
      integer lovl,kovl,lovl2,kovl2
      integer ival,iat,ish,iacf,iacl,iabf,iabl,llpov,ulpov,llov,ulov
      integer lpovl,kpovl
      integer ifirst,ilast,im,ityp,ncart,ngen,sphcart,nprim,ll,ul
      integer basis,ktmp,ltmp,nbf_max_at,nshl
      logical status


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
      if(.not.ma_push_get(mt_dbl,nbf_ao*nbf_ao,'overlap',lovl,kovl))
     $  call errquit("dftplusu_init: rtdb_get failed",0,RTDB_ERR)
      if(.not.ma_push_get(mt_dbl,nbf_ao*nbf_ao,'tmp',ltmp,ktmp))
     $  call errquit("dftplusu_init: ma_push_get failed",0,MA_ERR)

      if(.not.ma_push_get(mt_dbl,nbf_ao*nbf_ao,'povl',lpovl,kpovl))
     $  call errquit("dftplusu_init: rtdb_get failed",0,RTDB_ERR)


      if (mullpop) then
      else
        call ga_dgemm('n','n',nbf_ao,nbf_ao,nbf_ao,1d0,g_sm12,g_s,0d0,
     $                 g_sp12)
      endif

      basis = ao_bas_han
      nbf_max_at = nbf_ao_mxnbf_ce

      do ival=1,nvals
        size1 = 2*lvals(ival) + 1
        if(.not.ga_create(mt_dbl,nbf_ao,nbf_ao,"spps",0,0,g_spps(ival)))
     $  call errquit("dftplusu_init: ga_create failed",0,GA_ERR)   
        call ga_zero(g_spps(ival))
     
        if(ga_nodeid().eq.0) then

        iat = atoms(ival)
        status = bas_ce2cnr(basis,iat,iacf,iacl)
        status = bas_ce2bfr(basis,iat,iabf,iabl)
        if (mullpop) then
          call ga_get(g_s,1,nbf_ao,iabf,iabl,dbl_mb(kovl),nbf_ao)
        else
          call ga_get(g_sp12,1,nbf_ao,iabf,iabl,dbl_mb(kovl),nbf_ao)
        endif
        
        nshl = 0
        do ish=iacf,iacl
          status = bas_continfo(basis,ish,ityp,nprim,ngen,sphcart)
          if (ityp.ne.lvals(ival)) cycle

          nshl = nshl + 1
          status = bas_cn2bfr(basis,ish,ifirst,ilast)
          do im=0,size1-1
            llpov = kpovl + im*nbf_ao + (nshl-1)*size1*nbf_ao
            ulpov = llpov + nbf_ao
            llov = kovl + (ifirst-iabf+im)*nbf_ao
            ulov = llov + nbf_ao
            if (mullpop) then
              call ga_put(g_spps(ival),1,nbf_ao,ifirst+im,ifirst+im,
     $                    dbl_mb(llov),nbf_ao)         
            else
              call dcopy(nbf_ao,dbl_mb(llov),1,dbl_mb(llpov),1)
c              dbl_mb(llpov:ulpov) = dbl_mb(llov:ulov)
            endif
          enddo
        enddo

        if (.not.mullpop) then
         call dgemm('n','t',nbf_ao,nbf_ao,size1*nshl,1d0,dbl_mb(kpovl),
     $              nbf_ao,dbl_mb(kpovl),nbf_ao,0d0,dbl_mb(ktmp),
     $              nbf_ao)
         call ga_put(g_spps(ival),1,nbf_ao,1,nbf_ao,dbl_mb(ktmp),nbf_ao)
        endif

        endif
        
      enddo

      if(.not.ma_chop_stack(lovl))
     $ call errquit("dftplusu: ma_chop_stack failed",0,MA_ERR)

      return
      end
      
