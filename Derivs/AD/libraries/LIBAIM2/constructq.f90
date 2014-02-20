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

SUBROUTINE constructq( m, n, a, tau, q )

IMPLICIT NONE

INTEGER, INTENT(IN) :: m, n
DOUBLE PRECISION, INTENT(IN), DIMENSION(m,n) :: a
DOUBLE PRECISION, INTENT(IN), DIMENSION(MIN(m,n)) :: tau
DOUBLE PRECISION, INTENT(OUT), DIMENSION(m,n) :: q

INTEGER :: row, col
DOUBLE PRECISION, DIMENSION(MIN(m,n),MIN(m,n)) :: eye
DOUBLE PRECISION, DIMENSION(MIN(m,n),1) :: v
DOUBLE PRECISION, DIMENSION(1,MIN(m,n)) :: vt
DOUBLE PRECISION, DIMENSION(MIN(m,n),MIN(m,n)) :: vv
DOUBLE PRECISION, DIMENSION(m,MIN(m,n)) :: h
DOUBLE PRECISION, DIMENSION(m,n) :: qt
EXTERNAL matrixmult, transp


DO row=1,m
   DO col=1,n
      IF (row==col) THEN
         eye(row,col) = 1.d0
         q(row,col) = 1.d0
      ELSE
         eye(row,col) = 0.d0
         q(row,col) = 0.d0
      END IF
   END DO
END DO


DO col=1,MIN(m,n)
   v(1:col-1,1) = 0.d0
   v(col,1) = 1.d0
   v(col+1:m,1) = a(col+1:m,col)

   CALL transp(v,n,1,vt)

   CALL matrixmult(v,n,1,vt,n,vv)

   h = eye - tau(col) * vv

   CALL matrixmult(q,m,n,h,m,qt)

   q = qt

END DO

END SUBROUTINE
