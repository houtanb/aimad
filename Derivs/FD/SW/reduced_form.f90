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

SUBROUTINE reduced_form(q,qrows,qcols,neq,b,info)

IMPLICIT NONE

INTEGER, INTENT(IN) :: qrows, qcols, neq
DOUBLE PRECISION, DIMENSION(qrows,qcols), INTENT(IN) :: q
DOUBLE PRECISION, DIMENSION(neq,neq), INTENT(OUT) :: b
INTEGER, INTENT(OUT) :: info

! variables needed for dgesv
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: aq
INTEGER, ALLOCATABLE, DIMENSION(:) :: ipiv

! local variables
INTEGER :: status

EXTERNAL :: dgesv

ALLOCATE(aq(neq,neq),ipiv(neq),STAT = status)
b = q(:,1:neq)
aq = -q(:,neq+1:qcols)

CALL dgesv(neq,neq,aq,neq,ipiv,b,neq,info)


DEALLOCATE(aq,ipiv)

END SUBROUTINE
