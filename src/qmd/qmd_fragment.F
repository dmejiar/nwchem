      subroutine qmd_fragment(nat,coord,rcov,fragment,fragstep,stopnow,
     $                        dt)

      implicit none

      integer nat, istep, fragstep
      logical fragment, stopnow
      double precision coord(3,nat), rcov(nat), dt

      integer idfrag(nat)
      logical bonded(nat,nat), finished
      double precision distance
      integer iat, jat, nfrag, totat

      do iat=1, nat
        bonded(iat,iat) = .true.
        do jat=iat+1, nat
          distance = dsqrt(sum((coord(1:3,iat) - coord(1:3,jat))**2))
          if (distance.lt.rcov(iat)+rcov(jat)) then
            bonded(iat,jat) = .true.
            bonded(jat,iat) = .true.
          else
            bonded(iat,jat) = .false.
            bonded(jat,iat) = .false.
          endif
        enddo
      enddo

      nfrag = 0
      idfrag(:) = 0
      idfrag(1) = 1
      totat = 1
      nfrag = nfrag + 1
      do while (totat.lt.nat)
        finished = .false.
        do while (.not.finished)
          finished = .true.
          do iat=1,nat
            if (idfrag(iat).ne.nfrag) cycle
            do jat=1,nat
              if (idfrag(jat).ne.0) cycle
              if (.not.bonded(iat,jat)) cycle
              finished = .false.
              idfrag(jat) = nfrag
              totat = totat + 1
            enddo
          enddo
        enddo

        do iat=1,nat
          if (idfrag(iat).eq.0) then
            idfrag(iat) = nfrag + 1
            totat = totat + 1
            nfrag = nfrag + 1
            exit
          endif
        enddo
      enddo

      nfrag = maxval(idfrag)
      if (nfrag.eq.1) then
        if (fragment) then
          fragment = .false.
          fragstep = 0
        endif
      else
        fragment = .true.
        fragstep  = fragstep + 1
        if (fragstep*dt.ge.50.0d0) stopnow = .true.
      endif

      if (nfrag.gt.1) then
        write(*,*)
        write(*,*) '@! Fragmentation detected, ',fragstep
        if (stopnow) then
          do iat=1,nat
            write(*,*) iat,idfrag(iat)
          enddo
        endif
        write(*,*)
      endif

      end subroutine
