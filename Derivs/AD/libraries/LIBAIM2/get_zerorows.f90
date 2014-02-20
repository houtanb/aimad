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

SUBROUTINE get_zerorows(hpart,nrows,ncols,zerotol,zerorows,nz)

IMPLICIT NONE

INTEGER, INTENT(IN) :: nrows, ncols
DOUBLE PRECISION, INTENT(IN), DIMENSION(nrows,ncols) :: hpart
DOUBLE PRECISION, INTENT(IN) :: zerotol
LOGICAL, INTENT(OUT), DIMENSION(nrows) :: zerorows
INTEGER, INTENT(OUT) :: nz
INTEGER :: indxi, indxj
DOUBLE PRECISION :: row_sum

nz = 0
DO indxi=1,nrows

   row_sum = 0.0d0
   DO indxj=1,ncols
      row_sum = row_sum + ABS(hpart(indxi,indxj))
   END DO
   IF (row_sum .LE. zerotol) THEN
      zerorows(indxi) = .TRUE.
      nz = nz+1
   ELSE
      zerorows(indxi) = .FALSE.
   END IF

END DO

END SUBROUTINE


