      SUBROUTINE OFFSET_ccsdt_lambda1_18_5_1(l_a_offset,k_a_offset,size)
C     $Id$
C     This is a Fortran77 program generated by Tensor Contraction Engine v.1.0
C     Copyright (c) Battelle & Pacific Northwest National Laboratory (2002)
C     i2 ( h2 h10 p13 h12 h14 h16 )_vt
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "errquit.fh"
#include "tce.fh"
      INTEGER l_a_offset
      INTEGER k_a_offset
      INTEGER size
      INTEGER length
      INTEGER addr
      INTEGER h2b
      INTEGER h10b
      INTEGER p13b
      INTEGER h12b
      INTEGER h14b
      INTEGER h16b
      length = 0
      DO h2b = 1,noab
      DO h10b = 1,noab
      DO p13b = noab+1,noab+nvab
      DO h12b = 1,noab
      DO h14b = h12b,noab
      DO h16b = h14b,noab
      IF (int_mb(k_spin+h2b-1)+int_mb(k_spin+h10b-1)+int_mb(k_spin+p13b-
     &1) .eq. int_mb(k_spin+h12b-1)+int_mb(k_spin+h14b-1)+int_mb(k_spin+
     &h16b-1)) THEN
      IF (ieor(int_mb(k_sym+h2b-1),ieor(int_mb(k_sym+h10b-1),ieor(int_mb
     &(k_sym+p13b-1),ieor(int_mb(k_sym+h12b-1),ieor(int_mb(k_sym+h14b-1)
     &,int_mb(k_sym+h16b-1)))))) .eq. ieor(irrep_v,irrep_t)) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+h2b-1)+int_mb(k_spin+h10b-
     &1)+int_mb(k_spin+p13b-1)+int_mb(k_spin+h12b-1)+int_mb(k_spin+h14b-
     &1)+int_mb(k_spin+h16b-1).ne.12)) THEN
      length = length + 1
      END IF
      END IF
      END IF
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      IF (.not.MA_PUSH_GET(mt_int,2*length+1,'noname',l_a_offset,k_a_off
     &set)) CALL ERRQUIT('ccsdt_lambda1_18_5_1',0,MA_ERR)
      int_mb(k_a_offset) = length
      addr = 0
      size = 0
      DO h2b = 1,noab
      DO h10b = 1,noab
      DO p13b = noab+1,noab+nvab
      DO h12b = 1,noab
      DO h14b = h12b,noab
      DO h16b = h14b,noab
      IF (int_mb(k_spin+h2b-1)+int_mb(k_spin+h10b-1)+int_mb(k_spin+p13b-
     &1) .eq. int_mb(k_spin+h12b-1)+int_mb(k_spin+h14b-1)+int_mb(k_spin+
     &h16b-1)) THEN
      IF (ieor(int_mb(k_sym+h2b-1),ieor(int_mb(k_sym+h10b-1),ieor(int_mb
     &(k_sym+p13b-1),ieor(int_mb(k_sym+h12b-1),ieor(int_mb(k_sym+h14b-1)
     &,int_mb(k_sym+h16b-1)))))) .eq. ieor(irrep_v,irrep_t)) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+h2b-1)+int_mb(k_spin+h10b-
     &1)+int_mb(k_spin+p13b-1)+int_mb(k_spin+h12b-1)+int_mb(k_spin+h14b-
     &1)+int_mb(k_spin+h16b-1).ne.12)) THEN
      addr = addr + 1
      int_mb(k_a_offset+addr) = h16b - 1 + noab * (h14b - 1 + noab * (h1
     &2b - 1 + noab * (p13b - noab - 1 + nvab * (h10b - 1 + noab * (h2b 
     &- 1)))))
      int_mb(k_a_offset+length+addr) = size
      size = size + int_mb(k_range+h2b-1) * int_mb(k_range+h10b-1) * int
     &_mb(k_range+p13b-1) * int_mb(k_range+h12b-1) * int_mb(k_range+h14b
     &-1) * int_mb(k_range+h16b-1)
      END IF
      END IF
      END IF
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      RETURN
      END
