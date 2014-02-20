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

SUBROUTINE pseudoinv( a, neq, pseudoinvmat )

IMPLICIT NONE
INTEGER :: neq
DOUBLE PRECISION, INTENT(IN), DIMENSION(neq,neq) :: a
DOUBLE PRECISION, INTENT(OUT), DIMENSION(neq,neq) :: pseudoinvmat

INTEGER :: info, lwork, status, i
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: d, work
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: u, vt, recipd
EXTERNAL :: dgesvd


!
! Allocate memory
!
pseudoinvmat(:,:) = 0.0d0
lwork = MAX( 1, 3*neq+neq, 5*neq )
ALLOCATE( u(neq,neq), vt(neq,neq), recipd(neq,neq), d(neq), work(lwork), STAT=status )
IF (status .NE. 0) THEN
   WRITE(*,*)
   WRITE(*,*) '**** IN PSEUDOINV(): COULD NOT ALLOCATE MEMORY ****'
   RETURN
END IF


CALL dgesvd( 'A', 'A', neq, neq, a, neq, d, u, neq, vt, neq, work, lwork, info )
IF( info .NE. 0 ) THEN
   WRITE(*,*)
   WRITE(*,*) '**** IN PSEUDOINV():: Call to dgesvd() returned info != 0'
END IF
recipd(:,:) = 0.0d0
DO i=1,neq
   IF( ABS(d(i)) > 10e-12 ) THEN
      recipd(i,i) = 1/d(i)
   END IF
END DO

pseudoinvmat = MATMUL(MATMUL(TRANSPOSE(vt),recipd),TRANSPOSE(u))


!
! Deallocate memory
!
DEALLOCATE( u, vt, d, work, recipd, STAT=status )
IF (status .NE. 0) THEN
   WRITE(*,*)
   WRITE(*,*) '**** IN PSEUDOINV(): COULD NOT DEALLOCATE MEMORY ****'
   RETURN
END IF
END SUBROUTINE
