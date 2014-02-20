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

SUBROUTINE writematint(mat, rows, cols, strow, endrow, stcol, endcol, filename)

IMPLICIT NONE
INTEGER, INTENT(IN) :: strow, endrow, stcol, endcol, rows, cols
INTEGER, DIMENSION(rows,cols), INTENT(IN) :: mat
INTEGER :: i, j
CHARACTER(len=100), INTENT(IN) :: filename

OPEN (UNIT=4, FILE=filename, STATUS='OLD', ACTION='WRITE', POSITION='APPEND')

DO j=stcol,endcol
   DO i=strow,endrow
      WRITE (4,100) i,j,mat(i,j)
     END DO
END DO

100 FORMAT ('mat(', I4, ',', I4, ')=', I5, ';')

CLOSE(UNIT=4)
END SUBROUTINE
