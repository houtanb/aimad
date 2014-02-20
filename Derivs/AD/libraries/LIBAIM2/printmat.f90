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

SUBROUTINE printmat(mat,nrows,ncols)

IMPLICIT NONE

INTEGER, INTENT(IN) :: nrows,ncols
DOUBLE PRECISION, INTENT(IN), DIMENSION(nrows,ncols) :: mat

DOUBLE PRECISION :: zerotol = 0.000000000000001
INTEGER :: indxi, indxj


DO indxj = 1,ncols
   DO indxi = 1,nrows
      IF (ABS(mat(indxi,indxj))>zerotol) THEN
         WRITE(*,100) indxi,indxj,mat(indxi,indxj)
      END IF
   END DO
END DO

100 FORMAT('cofbf(',I4,',',I4,')',' =  ',ES30.16,';')

END SUBROUTINE
