      subroutine grad_dftplusu(rtdb, g_dens, g_wdens)
      implicit none

      integer rtdb, g_dens(2), g_wdens

#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "cdft.fh"
#include "dftplusu.fh"

      integer g_tmp1, g_tmp2, g_tmp3, g_s12, g_s, g_ppt, g_svecs
      integer lsvals,ksvals,me,ish,ityp,nprim,ngen,ival,ibf
      integer ilo,ihi,jlo,jhi,i,j,iabf,iabl,iacf,iacl, sphcart
      integer ifirst,ilast,iat
      logical status,ga_create_atom_blocked
      external ga_create_atom_blocked
      double precision toll_s,fij
      double precision factor
      parameter(toll_s=1d-10)

      me = ga_nodeid()

c     Create AO overlap matrix GA and get distribution
      g_s = ga_create_atom_blocked(geom, ao_bas_han, 'overlap')
      call ga_distribution(g_s,me,ilo,ihi,jlo,jhi)

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

c       Add contribution to energy-weighted density matrix
        if (mullpop) then
        ! Linear term
          call ga_dgemm('n','n',nbf_ao,nbf_ao,nbf_ao,factor,g_dens(1),
     $                   g_ppt,0d0,g_tmp1)
          call ga_symmetrize(g_tmp1)
          call ga_add(0.5d0*uvals(ival)/factor,g_tmp1,1d0,g_wdens,
     $                g_wdens)

        ! Idempotency term
          call ga_dgemm('n','n',nbf_ao,nbf_ao,nbf_ao,factor,g_dens(1),
     $                   g_ppt,0d0,g_tmp1)
          call ga_dgemm('n','n',nbf_ao,nbf_ao,nbf_ao,1d0,g_tmp1,g_s,0d0,
     $                   g_tmp2)
          call ga_dgemm('n','n',nbf_ao,nbf_ao,nbf_ao,2d0,g_tmp2,g_tmp1,
     $                   0d0,g_tmp3)
          call ga_symmetrize(g_tmp3)
          call ga_add(-0.5d0*uvals(ival)/factor,g_tmp3,1d0,g_wdens,
     $                g_wdens)

          if (ipol.gt.1) then
           ! Linear term
           call ga_dgemm('n','n',nbf_ao,nbf_ao,nbf_ao,1d0,g_dens(2),
     $                   g_ppt,0d0,g_tmp1)
           call ga_symmetrize(g_tmp1)
           call ga_add(0.5d0*uvals(ival),g_tmp1,1d0,g_wdens,g_wdens)

           ! Idempotency term
           call ga_dgemm('n','n',nbf_ao,nbf_ao,nbf_ao,factor,g_dens(2),
     $                   g_ppt,0d0,g_tmp1)
           call ga_dgemm('n','n',nbf_ao,nbf_ao,nbf_ao,1d0,g_tmp1,g_s,
     $                   0d0,g_tmp2)
           call ga_dgemm('n','n',nbf_ao,nbf_ao,nbf_ao,2d0,g_tmp2,g_tmp1,
     $                   0d0,g_tmp3)
           call ga_symmetrize(g_tmp3)
           call ga_add(-0.5d0*uvals(ival),g_tmp3,1d0,g_wdens,g_wdens)
          endif
        else
          ! TMP1 -> P * P^T * S^(1/2)
          call ga_dgemm('n','n',nbf_ao,nbf_ao,nbf_ao,1d0,g_ppt,g_s12,
     $                   0d0,g_tmp1)

          ! TMP2 -> P * P^T * S^(1/2) * D
          call ga_dgemm('n','n',nbf_ao,nbf_ao,nbf_ao,factor,g_tmp1,
     $                   g_dens(1),0d0,g_tmp2)

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

          if (ipol.gt.1) then
            ! TMP2 -> P * P^T * S^(1/2) * D
            call ga_dgemm('n','n',nbf_ao,nbf_ao,nbf_ao,1d0,g_tmp1,
     $                     g_dens(2),0d0,g_tmp2)

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
            call ga_dgemm('n','n',nbf_ao,nbf_ao,nbf_ao,-2d0,g_s12,
     $                     g_tmp2,1d0,g_tmp3)

            ! S -> 2*P*P^T*S^(1/2)*D * (I - S^(1/2) * P * P^T * S^(1/2) * D)
            call ga_dgemm('n','n',nbf_ao,nbf_ao,nbf_ao,2d0,g_tmp2,
     $                     g_tmp3,1d0,g_s)
          endif

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
          call ga_add(0.5d0*uvals(ival)/factor,g_s,1d0,g_wdens,g_wdens)
        endif
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


      end subroutine
