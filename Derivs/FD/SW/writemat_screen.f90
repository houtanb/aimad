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

SUBROUTINE writemat_screen(mat, rows, cols, strow, endrow, stcol, endcol)

IMPLICIT NONE
INTEGER, INTENT(IN) :: strow, endrow, stcol, endcol, rows, cols
DOUBLE PRECISION, DIMENSION(rows,cols), INTENT(IN) :: mat
INTEGER :: i, j
DOUBLE PRECISION :: zerotol = 0.000000000000001
CHARACTER(len=500) :: output


DO j=stcol,endcol
   DO i=strow,endrow
      IF (ABS(mat(i,j))>zerotol) THEN
         WRITE (*,100) i,j,mat(i,j)
		 !CALL mexPrintf(trim(output)//char(13))
      END IF
     END DO
END DO

100 FORMAT ('mat(', I4, ',', I4, ')=', ES30.15, ';')


END SUBROUTINE
