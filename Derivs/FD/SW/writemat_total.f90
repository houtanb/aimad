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

SUBROUTINE writemat_total(mat, rows, cols, filename, filenamelen)

IMPLICIT NONE
INTEGER, INTENT(IN) :: rows, cols
DOUBLE PRECISION, DIMENSION(rows,cols), INTENT(IN) :: mat
INTEGER :: i, j
INTEGER, INTENT(IN) :: filenamelen
CHARACTER(len=filenamelen), INTENT(IN) :: filename
DOUBLE PRECISION :: zerotol = 0.000000000000001

OPEN (UNIT=4, FILE=filename, STATUS='REPLACE', ACTION='WRITE')

WRITE (4,*) 'mat=zeros(', rows, ',', cols, ');'
DO j=1,cols
   DO i=1,rows
      IF (ABS(mat(i,j))>zerotol) THEN
         WRITE (4,100) i,j,mat(i,j)
      END IF
     END DO
END DO
100 FORMAT ('mat(', I4, ',', I4, ')=', ES30.15, ';')

CLOSE(UNIT=4)
END SUBROUTINE
