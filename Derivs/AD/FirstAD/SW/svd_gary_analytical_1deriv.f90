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

SUBROUTINE svd_gary_analytical_1deriv( a, da, neq, y, dy, lambda, x, nbdirs )

!!
!!
!! FROM: Anderson, Gary. Journal of Economic Dynamics and Control 11 (1987) 465-481
!!       specifically, pp. 473-474
!!
!!

IMPLICIT NONE

INTEGER, INTENT(IN) :: neq, nbdirs
DOUBLE PRECISION, INTENT(IN), DIMENSION(neq,neq) :: a, y, x
DOUBLE PRECISION, INTENT(IN), DIMENSION(neq) :: lambda

DOUBLE PRECISION, INTENT(IN), DIMENSION(nbdirs,neq,neq) :: da
DOUBLE PRECISION, INTENT(OUT), DIMENSION(nbdirs,neq,neq) :: dy

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: invlambda
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: da_dir, theta1_j, bigtheta

DOUBLE PRECISION :: zerotol

INTEGER :: nnz, nz, status, i, j

dy(:,:,:) = 0.0d0
zerotol = 1.0e-12

nnz = 0

DO i=1,neq
   IF( lambda(i) .GT. zerotol ) THEN
      nnz = nnz + 1
   END IF
END DO

nz = neq - nnz

ALLOCATE( da_dir(neq,neq), bigtheta(neq,neq), invlambda(nnz,nnz), &
          theta1_j(nnz,nz), STAT=status )
IF( status .NE. 0 ) THEN
   WRITE(*,*) 'ERROR in svd_gary_analytical_2deriv(): Could not allocate memory'
   RETURN
END IF

invlambda = 0.d0

DO i=1,nnz
   invlambda(i,i) = 1.d0/lambda(i)
END DO


DO j=1,nbdirs
   theta1_j = 0.d0
   da_dir(:,:) = da(j,:,:)
   theta1_j = MATMUL(MATMUL(MATMUL(invlambda,TRANSPOSE(y(:,nz+1:neq))),da_dir),x(:,1:nz))

   bigtheta = 0.d0
   bigtheta(1:nz,nz+1:neq) = -1.0d0*TRANSPOSE(theta1_j)
   bigtheta(nz+1:neq,1:nz) = theta1_j
   dy(j,:,:) = MATMUL( y, bigtheta )

END DO

DEALLOCATE( invlambda, da_dir, theta1_j, bigtheta, STAT=status )
IF( status .NE. 0 ) THEN
   WRITE(*,*) 'ERROR in svd_gary_analytical_2deriv(): Could not deallocate memory'
   RETURN
END IF
END SUBROUTINE
