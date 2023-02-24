      subroutine grad_dftplusu(rtdb, g_dens, g_wdens)
      implicit none

      integer natom, rtdb, g_dens(2)

#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "cdft.fh"
#include "dftplusu.fh"

      integer g_wdens
      integer g_tmp1, g_tmp2, g_tmp3, g_s12, g_s, g_ppt, g_svecs
      integer lsvals,ksvals,me,ish,ityp,nprim,ngen,ival,ibf
      integer ilo,ihi,jlo,jhi,i,j,iabf,iabl,iacf,iacl, sphcart
      integer ifirst,ilast,iat,isp
      logical status,ga_create_atom_blocked
      external ga_create_atom_blocked
      double precision toll_s,fij
      double precision factor
      parameter(toll_s=1d-10)

      me = ga_nodeid()

c     Create AO overlap matrix GA and get distribution
      g_s = ga_create_atom_blocked(geom, ao_bas_han, 'overlap')
      g_wdens = ga_create_atom_blocked(geom, ao_bas_han, 'wdens')
      call ga_distribution(g_s,me,ilo,ihi,jlo,jhi)
      call ga_zero(g_wdens)

c     Compute AO overlap matrix
      call ga_zero(g_s)
      call int_1e_ga(ao_bas_han, ao_bas_han, g_s, 'overlap', oskel)
      if(oskel) call sym_symmetrize(geom,ao_bas_han,.false.,g_s)

c     Create GAs
      status = ga_duplicate(g_s,g_tmp1,"tmp1") .and.
     $         ga_duplicate(g_s,g_tmp2,"tmp2") .and.
     $         ga_duplicate(g_s,g_tmp3,"tmp3")
      if(.not.status)
     $  call errquit("grad_dftplusu: GA duplicate failed",0,GA_ERR) 

c     Compute S^(1/2)
      if (.not.mullpop) then
        ! Allocate Eigenvalue space
        if(.not.ma_push_get(mt_dbl, nbf_ao,'s eigvals',lsvals,ksvals))
     $    call errquit('grad_dftplusu: cannot allocate s eigs',0,MA_ERR)
        call dfill(nbf_ao,0d0,dbl_mb(ksvals),1)

        ! Create GAs
        if(.not.(ga_duplicate(g_s,g_svecs,'svecs').and.
     $           ga_duplicate(g_s,g_s12,'s12')))
     $    call errquit("grad_dftplusu: could not duplicate GA",0,GA_ERR)
        call ga_zero(g_svecs)
        call ga_zero(g_s12)

        ! Diagonalize overlap matrix
#if defined(PARALLEL_DIAG)
        call ga_diag_std(g_s,g_svecs,dbl_mb(ksvals))
#else
        call ga_diag_std_seq(g_s,g_svecs,dbl_mb(ksvals))
#endif

        ! Build S^(1/2)
        call diis_bld12(toll_s,dbl_mb(ksvals),g_svecs,g_s12,g_tmp1,3)
      endif

c     Loop over all U values
      factor = 1d0
      if (ipol.eq.1) factor = 0.5d0
      g_ppt = ga_create_atom_blocked(geom, ao_bas_han, 'ppt')
      do ival=1,nvals
        call ga_zero(g_ppt)

c       Build Projector P*P^T
        iat = atoms(ival)
        status = bas_ce2cnr(ao_bas_han,iat,iacf,iacl)
        status = bas_ce2bfr(ao_bas_han,iat,iabf,iabl)
        do ish=iacf,iacl
          status = bas_continfo(ao_bas_han,ish,ityp,nprim,ngen,
     $                          sphcart)
          if (ityp.ne.lvals(ival)) cycle
          status = bas_cn2bfr(ao_bas_han,ish,ifirst,ilast)
          do ibf=ifirst,ilast
            if ((ibf.ge.ilo).and.(ibf.le.ihi).and.(ibf.ge.jlo).and.
     $          (ibf.le.jhi)) then
              call ga_put(g_ppt,ibf,ibf,ibf,ibf,1d0,1)
            endif
          enddo
        enddo
        call ga_sync()
c
        do isp=1,ipol

c       Add contribution to energy-weighted density matrix
        if (mullpop) then
        ! Linear term
          call ga_dgemm('n','n',nbf_ao,nbf_ao,nbf_ao,factor,g_dens(isp),
     $                   g_ppt,0d0,g_tmp1)
          call ga_symmetrize(g_tmp1)
          call ga_add(-0.5d0*uvals(ival)/factor,g_tmp1,1d0,g_wdens,
     $                g_wdens)

        ! Idempotency term
          call ga_dgemm('n','n',nbf_ao,nbf_ao,nbf_ao,factor,g_dens(isp),
     $                   g_ppt,0d0,g_tmp1)
          call ga_dgemm('n','n',nbf_ao,nbf_ao,nbf_ao,1d0,g_tmp1,g_s,0d0,
     $                   g_tmp2)
          call ga_dgemm('n','n',nbf_ao,nbf_ao,nbf_ao,2d0,g_tmp2,g_tmp1,
     $                   0d0,g_tmp3)
          call ga_symmetrize(g_tmp3)
          call ga_add(0.5d0*uvals(ival)/factor,g_tmp3,1d0,g_wdens,
     $                g_wdens)

        else
          ! TMP1 -> P * P^T * S^(1/2)
          call ga_dgemm('n','n',nbf_ao,nbf_ao,nbf_ao,1d0,g_ppt,g_s12,
     $                   0d0,g_tmp1)

          ! TMP2 -> P * P^T * S^(1/2) * D
          call ga_dgemm('n','n',nbf_ao,nbf_ao,nbf_ao,factor,g_tmp1,
     $                   g_dens(isp),0d0,g_tmp2)

          ! TMP3 -> I - S^(1/2) * P * P^T * S^(1/2) * D
          call ga_zero(g_tmp3)
          do i=ilo,ihi
            do j=jlo,jhi
              if (i.eq.j) then
                call ga_put(g_tmp3,i,i,i,i,1d0,1)
              endif
            enddo
          enddo
          call ga_sync()
          call ga_dgemm('n','n',nbf_ao,nbf_ao,nbf_ao,-2d0,g_s12,g_tmp2,
     $                   1d0,g_tmp3)

          ! S -> 2*P*P^T*S^(1/2)*D * (I - S^(1/2) * P * P^T * S^(1/2) * D)
          call ga_dgemm('n','n',nbf_ao,nbf_ao,nbf_ao,2d0,g_tmp2,g_tmp3,
     $                   0d0,g_s)

          ! Symmetrize
          call ga_symmetrize(g_s)

          ! First Eigenvector Pass
          call ga_dgemm('n','n',nbf_ao,nbf_ao,nbf_ao,1d0,g_s,g_svecs,
     $                   0d0,g_tmp2)
          call ga_dgemm('t','n',nbf_ao,nbf_ao,nbf_ao,1d0,g_svecs,g_tmp2,
     $                   0d0,g_s)

          ! Scale with Eigenvalues
          do i=ilo,ihi
            do j=jlo,jhi
              call ga_get(g_s,i,i,j,j,fij,1)
              fij = fij/
     $           (sqrt(dbl_mb(ksvals+i-1))+sqrt(dbl_mb(ksvals+j-1)))
              call ga_put(g_s,i,i,j,j,fij,1)
            enddo
          enddo

          ! Last eigenvectors pass
          call ga_sync()
          call ga_dgemm('n','n',nbf_ao,nbf_ao,nbf_ao,1d0,g_svecs,g_s,
     $                   0d0,g_tmp1)
          call ga_dgemm('n','t',nbf_ao,nbf_ao,nbf_ao,1d0,g_tmp1,g_svecs,
     $                   0d0,g_s)   
          call ga_add(-0.5d0*uvals(ival)/factor,g_s,1d0,g_wdens,
     $                 g_wdens)
        endif
        enddo
      enddo

      status = ga_destroy(g_tmp1) .and. ga_destroy(g_tmp2) .and.
     $         ga_destroy(g_tmp3) .and. ga_destroy(g_s) .and.
     $         ga_destroy(g_ppt) 
      if (.not.mullpop) then
        status = status .and.
     $           ga_destroy(g_svecs) .and. ga_destroy(g_s12)   
        if(.not.ma_pop_stack(lsvals))
     $    call errquit("grad_dftplusu: could not pop stack",0,MA_ERR) 
      endif

      if (.not.status)
     $  call errquit("grad_dftplusu: GA destroy failed",0,GA_ERR) 

      return
      end



      subroutine grad1_dftplusu(rtdb, geom, basis, nat, g_wdens, force, 
     $                          max_at_bf, nproc, oskel,
     $                          wdens, H, lbuf, scr, lscr)
      implicit none
      integer rtdb, g_wdens, nat, lbuf, lscr
      integer basis, geom, nproc, max_at_bf
      double precision force(3,nat), wdens(max_at_bf,max_at_bf)
      double precision H(*), scr(*)
      logical oskel

#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "msgids.fh"

      integer ijatom, next, iat1, iat2, iat3, ish1, ish2,
     $     iab1f, iab1l, iab2f, iab2l, iac1f, iac1l, iac2f, iac2l,
     $     if1, il1, if2, il2,
     $     icart, ic, nint, ip1, ip2

      double precision crd1, crd2 ! atomic coordinates
      dimension crd1(3), crd2(3)

      integer idatom
      dimension idatom(2)

      double precision dE, dx, dy, dz, qfac, fact, q1, q2

      logical status, pointforce,dobq

      character*16 name

      integer nxtask, task_size
      external nxtask

      task_size = 1
      call hf_print_set(1)

      ijatom = -1
      next = nxtask(nproc,task_size)
      do 90, iat1 = 1, nat
        do 80, iat2 = 1, iat1

          ijatom = ijatom + 1
          if ( ijatom .eq. next ) then

            status = bas_ce2bfr(basis,iat1,iab1f,iab1l)
            status = bas_ce2bfr(basis,iat2,iab2f,iab2l)

            if (iab1f.le.0 .or. iab2f.le.0) then
c     
c     At least one center has no functions on it ... next atom
c
              goto 1010
            endif

            if (oskel) then
               if (.not. sym_atom_pair(geom, iat1, iat2, qfac))
     $              goto 1010
            else
               qfac = 1.0d0
            endif

            status = bas_ce2cnr(basis,iat1,iac1f,iac1l)
            status = bas_ce2cnr(basis,iat2,iac2f,iac2l)

            call ga_get(g_wdens,iab1f,iab1l,iab2f,iab2l,wdens,max_at_bf)
c
            do 70, ish1 = iac1f, iac1l
              if ( iat1.eq.iat2 ) iac2l = ish1
              do 60, ish2 = iac2f, iac2l

C               shell block in atomic (D/Dw)-matrix block
                status = bas_cn2bfr(basis,ish1,if1,il1)
                if1 = if1 - iab1f + 1
                il1 = il1 - iab1f + 1
                status = bas_cn2bfr(basis,ish2,if2,il2)
                if2 = if2 - iab2f + 1
                il2 = il2 - iab2f + 1

                nint = ( il1 - if1 + 1 ) * ( il2 - if2 + 1 )

C               overlap derivatives
                call intd_1eov(basis,ish1,basis,ish2,lscr,scr,
     &               lbuf,H,idatom)

C     Dw x S

                if ( idatom(1) .ge. 1 ) then
C               idatom(1).ge.0 <=> idatom(2).ge.0 (no check necessary)
                  ic = 1
                  do 28, icart = 1, 3
                    de = 0.D0
                    do 22, ip1 = if1, il1
                      do 20, ip2 = if2, il2
                        dE = dE + wdens(ip1,ip2) * H(ic)
                        ic = ic + 1
 20                   continue
 22                 continue
                    dE = dE * qfac
                    force(icart,idatom(1)) = force(icart,idatom(1)) 
     $                                      - dE - dE
                    force(icart,idatom(2)) = force(icart,idatom(2)) 
     $                                      + dE + dE
 28               continue
                endif

 60           continue
 70         continue

 1010       continue

            next = nxtask(nproc,task_size)
          endif

 80     continue
 90   continue
      next = nxtask(-nproc,task_size)

      call ga_dgop(msg_grad_xc, force,  3*nat, '+')

      end subroutine
