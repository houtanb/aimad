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

SUBROUTINE inverse_d( a, ad, lda )

IMPLICIT NONE
INTEGER, INTENT(IN) :: lda
DOUBLE PRECISION, INTENT(INOUT), DIMENSION( lda, lda ) :: a, ad

!
! Local vars
INTEGER :: status, lwork, info
INTEGER, ALLOCATABLE, DIMENSION(:) :: ipiv, iwork
CHARACTER*1 :: norm
DOUBLE PRECISION :: anorm, rcond, eps, zerotol, dlange
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: work

EXTERNAL :: dlange, dgetrf, dgecon, dgetri

!
! Assign constants and allocate memory
!
eps = 10e-14
norm = 'I'
lwork = 5*lda
ALLOCATE( work(lwork), ipiv(lda), iwork(lda), STAT=status )
IF (status .NE. 0) THEN
   WRITE(*, *) '**** IN INVERSE_D(): COULD NOT ALLOCATE MEMORY ****'
   RETURN
END IF



!
! Calculate derivative using same calls as in Matlab
!

! norm(A)
anorm = dlange( norm, lda, lda, a, lda, work )
zerotol = 2*anorm*eps

! LU Decomposition
CALL dgetrf( lda, lda, a, lda, ipiv, info )
IF( info .NE. 0 ) THEN
!   WRITE(*,*)
!   WRITE(*,*) '**** IN INVERSE_D():: Call to dgetrf() returned info != 0'
END IF

! Reciprocal of condition number
CALL dgecon( norm, lda, a, lda, anorm, rcond, work, iwork, info )
IF( info .NE. 0 ) THEN
!   WRITE(*,*)
!   WRITE(*,*) '**** IN INVERSE_D():: Call to dgecon() returned info != 0'
END IF

!IF( rcond > zerotol ) THEN
IF( .TRUE. ) THEN
   ! A^-1
   CALL dgetri( lda, a, lda, ipiv, work, lwork, info )
   IF( info .NE. 0 ) THEN
!      WRITE(*,*)
!      WRITE(*,*) '**** IN INVERSE_D():: Call to dgetri() returned info != 0'
   END IF

   !
   ! Calculate derivative of f(x) = inv(a) ...
   ! from Edelman ( 18.325, Handout #2, Example #4 )
   !

   ! f'(x) = -inv(a) * da * inv(a)
   ad = MATMUL( -a, MATMUL( ad, a ) )
ELSE
   WRITE(*,*)
   WRITE(*,*) '**** IN INVERSE_D():: RCOND too large !!'
   WRITE(*,*) 'RCOND = ',rcond
   WRITE(*,*) 'ZEROTOL = ',zerotol
   WRITE(*,*) 'ANORM = ',anorm
   WRITE(*,*) '**** IN INVERSE_D():: INV(A) NOT COMPUTED !!'
   a(:,:) = 0.0d0
   ad(:,:) = 0.0d0
   RETURN
END IF

DEALLOCATE( work, ipiv, iwork, STAT=status )
IF (status .NE. 0) THEN
   WRITE(*, *) '**** IN INVERSE_D(): COULD NOT DEALLOCATE MEMORY ****'
   RETURN
END IF

END SUBROUTINE
