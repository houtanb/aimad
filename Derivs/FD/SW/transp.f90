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

SUBROUTINE transp( a, m, n, at )

IMPLICIT NONE

INTEGER, INTENT(IN) :: m,n
DOUBLE PRECISION, INTENT(IN), DIMENSION(m,n) :: a
DOUBLE PRECISION, INTENT(OUT), DIMENSION(n,m) :: at

INTEGER :: i, j

DO i=1,m
   DO j=1,n
      at(j,i) = a(i,j)
   END DO
END DO

END SUBROUTINE
