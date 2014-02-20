!
! Copyright (C) 2006-2014 Houtan Bastani and Luca Guerrieri
!
!
! This free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! It is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
!

SUBROUTINE exact_shift(h,q,iq,qrows,qcols,neq,nexact)

IMPLICIT NONE

INTEGER, INTENT(IN) :: qrows, qcols, neq
INTEGER, INTENT(INOUT) :: iq
INTEGER, INTENT(OUT) :: nexact
DOUBLE PRECISION, DIMENSION(neq,neq*3), INTENT(INOUT) :: h
DOUBLE PRECISION, DIMENSION(qrows,qcols), INTENT(INOUT) :: q

LOGICAL, ALLOCATABLE, DIMENSION(:) :: zerorows
INTEGER, ALLOCATABLE, DIMENSION(:) :: selector
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: hpart
INTEGER status
INTEGER lleft,rleft
INTEGER lright,rright
INTEGER nz

DOUBLE PRECISION :: zerotol

EXTERNAL get_zerorows, get_selector, shiftright


nexact = 0
lleft = 1
rleft = qcols
lright = qcols+1
rright = qcols+neq

zerotol = 0.0000000001d0

ALLOCATE(zerorows(neq),STAT=status)
CALL get_zerorows(h(:,lright:rright),neq,neq,zerotol,zerorows,nz)
!WRITE(*,*) "zerorows"
!WRITE(*,*) zerorows

DO
   IF ( ANY(zerorows) .AND. (iq .LE. qrows) ) THEN
      ALLOCATE(selector(nz),STAT=status)
      CALL get_selector(zerorows,neq,selector,nz)


      q(iq+1:iq+nz,:) = h(selector,lleft:rleft)
      iq = iq +nz
      nexact = nexact+nz
!      write(*,*) 'selector'
!      write(*,*) selector
      ALLOCATE(hpart(nz,3*neq),STAT=status)
      hpart(:,:) = 0.0d0
      CALL shiftright(h(selector,:),nz,3*neq,neq,hpart)
      h(selector,:) = hpart
      DEALLOCATE(hpart,STAT=status)

      CALL get_zerorows(h(:,lright:rright),neq,neq,zerotol,zerorows,nz)
!      WRITE(*,*) "ZEROROWS"
!      WRITE(*,*) zerorows
      DEALLOCATE(selector,STAT=status)
   ELSE
      EXIT
   END IF
END DO

DEALLOCATE(zerorows, STAT=status)

END SUBROUTINE

