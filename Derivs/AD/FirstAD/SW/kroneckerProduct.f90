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

SUBROUTINE kroneckerproduct( a, arows, acols, b, brows, bcols, c )

IMPLICIT NONE
INTEGER, INTENT(IN) :: arows, acols, brows, bcols
DOUBLE PRECISION, INTENT(IN), DIMENSION(arows,acols) :: a
DOUBLE PRECISION, INTENT(IN), DIMENSION(brows,bcols) :: b
DOUBLE PRECISION, INTENT(OUT), DIMENSION(arows*brows,acols*bcols) :: c

INTEGER :: i, j

DO j=1,acols
   DO i=1,arows
      c( (i-1)*brows+1:i*brows, (j-1)*bcols+1:j*bcols ) = a( i, j ) * b
   END DO
END DO

END SUBROUTINE
