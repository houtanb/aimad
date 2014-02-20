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

FUNCTION checkaim( neq, nlag, nlead, h, b )

IMPLICIT NONE

INTEGER :: neq, nlag, nlead
DOUBLE PRECISION, DIMENSION( neq, neq ) :: b
DOUBLE PRECISION, DIMENSION( neq, neq*3 ) :: h

! Local vars
INTEGER :: status, integer, i, j, info
INTEGER, ALLOCATABLE, DIMENSION( : ) :: rows, minus, plus, ipiv
DOUBLE PRECISION, ALLOCATABLE, DIMENSION( :, : ) :: b_new, eye, q, qshift, a, x, err
DOUBLE PRECISION :: checkaim

ALLOCATE( b_new(neq,neq*2), eye(neq,neq), &
          q(neq*(nlead+1), neq*(nlag+1+nlead)), &
          qshift(neq, neq*(nlag+1+nlead)), &
          rows(neq), minus(neq*nlag), plus((neq*(nlag+1+nlead)-(neq*nlag+1)+1)), &
          a(neq*(nlead+1), (neq*(nlag+1+nlead)-(neq*nlag+1)+1)), x(neq*(nlead+1),neq*nlag), &
          ipiv(neq*(nlead+1)), err(neq, neq*nlag), STAT=status )
IF( status .NE. 0 ) THEN
   WRITE(*,*) 'In checkaim() :: could not allocate memory'
   RETURN
END IF

checkaim = 0.0d0
CALL eyemat( eye, neq )

q(:,:) = 0.0d0
q(1:neq,1:neq) = b
q(1:neq,neq+1:neq*2) = -1.0d0*eye

DO i=1,nlead
   qshift(:,:) = 0.0d0

   DO j=1,neq
      rows(j) = i*neq + j
   END DO

   CALL shiftright( q(rows - neq, :), neq, neq*(nlag+1+nlead), neq, qshift )
   q( rows, : ) = qshift
END DO

DO i=1,neq*nlag
   minus(i) = i
END DO

j=1
DO i=(neq*nlag+1),(neq*(nlag+1+nlead))
   plus(j) = i
   j = j + 1
END DO


x = q( :, minus )
a = q( :, plus  )
CALL dgesv( neq*(nlead+1), neq*nlag, a, neq*(nlead+1), ipiv, x, neq*(nlead+1), info )
IF( info .NE. 0 ) THEN
   WRITE(*,*) 'In checkaim() :: could not calculate inverse'
END IF
q( :, minus ) = -1.0d0 * x

err = h(:,minus) + MATMUL( h(:,plus), q(:,minus) )
DO i=1,neq
   DO j=1,neq*nlag
      checkaim = MAX( checkaim, err(i,j) )
   END DO
END DO

DEALLOCATE( b_new, eye, q, qshift, rows, minus, plus, a, x, ipiv, err, STAT=status )
IF( status .NE. 0 ) THEN
   WRITE(*,*) 'In checkaim() :: could not deallocate memory'
   RETURN
END IF

RETURN
END FUNCTION
