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

SUBROUTINE loadcof(cof,neq,filename)

IMPLICIT NONE

! define arguments
INTEGER, INTENT(IN) :: neq
DOUBLE PRECISION, DIMENSION(neq,neq*3), INTENT(OUT) :: cof

! define working variables
CHARACTER(len=100), INTENT(IN) :: filename

DOUBLE PRECISION :: matsum
INTEGER :: indxi, indxj

OPEN (UNIT=4, FILE=filename, STATUS='OLD', ACTION='READ')

WRITE (*,*) "Opened file:"
WRITE (*,*) filename

matsum = 0.

! now read the various pieces of the sparse matrix
DO indxj=1,neq*3
   DO indxi=1,neq
   READ (4,100) cof(indxi,indxj)
   matsum = matsum + cof(indxi,indxj)
   END DO
END DO

CLOSE (UNIT=4)

WRITE (*,*) 'The sum is: '
WRITE (*,100) matsum

100 FORMAT (ES30.16)

END SUBROUTINE
