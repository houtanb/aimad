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

SUBROUTINE get_zerocols(hpart,nrows,ncols,zerotol,zerocols,nz)


! zerocols needs to be initialized to all false before the first call


IMPLICIT NONE

! declaration of arguments
INTEGER, INTENT(IN) :: nrows, ncols
DOUBLE PRECISION, INTENT(IN), DIMENSION(nrows,ncols) :: hpart
DOUBLE PRECISION, INTENT(IN) :: zerotol
LOGICAL,  DIMENSION(nrows), INTENT(INOUT) :: zerocols
INTEGER, INTENT(INOUT) :: nz


! local variables
INTEGER :: indxi, indxj
DOUBLE PRECISION :: col_sum

DO indxj=1,ncols
   IF (.NOT. zerocols(indxj)) THEN
      col_sum = 0.0d0
      DO indxi=1,nrows
         IF (.NOT. zerocols(indxi)) THEN
            !skip the rows correspending to cols previously eliminated
            col_sum = col_sum + ABS(hpart(indxi,indxj))
         END IF
      END DO
      IF (col_sum .LE. zerotol) THEN
         zerocols(indxj) = .TRUE.
         nz = nz+1
      END IF
   END IF
END DO

END SUBROUTINE
