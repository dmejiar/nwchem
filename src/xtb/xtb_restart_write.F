      subroutine xtb_restart_write(mol, calc, wfn)
      use mctc_io, only : structure_type
      use tblite_xtb_calculator, only: xtb_calculator
      use tblite_wavefunction_type, only : wavefunction_type
      implicit none
#include "xtb.fh"

      type(structure_type) :: mol
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn

      integer i, unitno, ispin, iat
      parameter(unitno=66)

      open(unit=unitno,status='unknown',form='unformatted',
     $     file=restart_out,position='rewind')
      do ispin=1,nspin
c        do i=1,nao
c          call swrite(unitno, dmat(1,i,ispin), int(calc%bas%nao))
c        enddo
        call swrite(unitno, wfn%qsh(1,ispin), int(calc%bas%nsh))
        call swrite(unitno, wfn%qat(1,ispin), int(mol%nat))
        call swrite(unitno, wfn%dpat(1,1,ispin), int(3*mol%nat))
        call swrite(unitno, wfn%qpat(1,1,ispin), int(6*mol%nat))
      enddo
      call util_flush(unitno)
      close(unit=unitno,status='keep')

      open(unit=unitno,status='unknown',form='formatted',
     $     file='charge',position='rewind')
      do iat=1,mol%nat
        write(unitno,*) wfn%qat(iat,1)
      enddo
      call util_flush(unitno)
      close(unit=unitno,status='keep')
      end
