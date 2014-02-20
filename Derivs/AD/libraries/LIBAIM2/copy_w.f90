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

SUBROUTINE copy_w(q,qrows,qcols,w,ia,iq,js)

IMPLICIT NONE

INTEGER, INTENT(IN) :: qrows, qcols, ia, iq
DOUBLE PRECISION, DIMENSION(qrows,qcols), INTENT(INOUT) :: q
DOUBLE PRECISION, DIMENSION(ia,ia), INTENT(IN) :: w
INTEGER, DIMENSION(ia), INTENT(IN) :: js

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: wt

INTEGER :: status

EXTERNAL transp
ALLOCATE(wt(qrows-iq,ia),STAT=status)

IF (iq<qrows) THEN
   CALL transp(w(:,ia-(qrows-iq)+1:ia),ia,qrows-iq,wt)
   q(iq+1:qrows,js) = wt
END IF

DEALLOCATE(wt,STAT=status)
END SUBROUTINE

