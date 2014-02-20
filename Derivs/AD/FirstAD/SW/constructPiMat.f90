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

SUBROUTINE constructpimat( pi, pirows, picols, mat )

IMPLICIT NONE
INTEGER, INTENT(IN) :: pirows, picols
DOUBLE PRECISION, INTENT(IN), DIMENSION(pirows,picols) :: pi
DOUBLE PRECISION, INTENT(OUT), DIMENSION(pirows+picols,pirows+picols) :: mat

INTEGER :: n

n = pirows + picols

mat(:,:) = 0.0d0
mat( picols+1:n, 1:picols ) = pi
mat( 1:picols, picols+1:n ) = -1.0d0 * TRANSPOSE(pi)

END SUBROUTINE
