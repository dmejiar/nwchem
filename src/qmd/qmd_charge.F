      subroutine qmd_charge(nat,istep,charges)
      implicit none
      integer nat,istep,iat,ipos,npos,unitno
      double precision charges(nat,50), charge(nat)
      parameter(unitno=66)

      open(unit=unitno,status='unknown',form='formatted',
     $     file='charge',position='rewind')
      do iat=1,nat
        read(unitno,*) charge(iat)
      enddo
      close(unit=unitno,status='keep')

      npos = min(istep,50)
      ipos = mod(istep,50)
      if (ipos.eq.0) ipos = 50
      charges(:,ipos) = charge

      open(unit=unitno,status='unknown',form='formatted',
     $     file='charges',position='rewind')
      do iat=1,nat
        write(unitno,*) sum(charges(iat,1:npos))/dble(npos)
      enddo
      call util_flush(unitno)
      close(unit=unitno,status='keep')
      end subroutine

