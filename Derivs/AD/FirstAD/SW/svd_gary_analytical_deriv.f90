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

SUBROUTINE svd_gary_analytical_deriv( a, da, neq, y, dy, lambda, x, nbdirs )
USE DIFFSIZES
IMPLICIT NONE
INTEGER, INTENT(IN) :: neq, nbdirs

DOUBLE PRECISION, INTENT(IN), DIMENSION(neq,neq) :: a, y, x
DOUBLE PRECISION, INTENT(IN), DIMENSION(nbdirsmax,neq,neq) :: da
DOUBLE PRECISION, INTENT(IN), DIMENSION(neq) :: lambda
DOUBLE PRECISION, INTENT(OUT), DIMENSION(nbdirsmax,neq,neq) :: dy

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: y1, y2, x1, x2, invlambda
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: da_dir, theta, themat
INTEGER :: zerotol, nnz, nz, status, i

EXTERNAL :: dgesv, eyemat


dy(:,:,:) = 0.0d0
zerotol = 1.0e-12

nnz = 0
DO i=1,neq
   IF( lambda(i) .GT. zerotol ) THEN
      nnz = nnz + 1
   ELSE
      EXIT
   END IF
END DO

nz = neq - nnz

ALLOCATE( y1(neq,nz), y2(neq,nnz), &
     x1(neq,nz), x2(neq,nnz), &
     da_dir(neq,neq), &
     invlambda(nnz,nnz), &
     theta(nnz,nz), &
     themat(neq,neq), STAT=status )
IF( status .NE. 0 ) THEN
   WRITE(*,*) 'ERROR in svd_gary_analytical_deriv(): Could not allocate memory'
   RETURN
END IF
!y1 = y( :, 1:nz )
y2 = y( :, nz+1:neq )

x1 = x( :, 1:nz )
!x2 = x( :, nz+1:neq )

invlambda=0.d0

DO i=1,nnz
   invlambda(i,i) = 1.0d0/lambda(i)
END DO


DO i=1,nbdirs
   da_dir(:,:) = da(i,:,:)

   theta(:,:) = 0.0d0
   !theta = MATMUL( invlambda, MATMUL( TRANSPOSE(y2), MATMUL( da_dir, x1 ) ) )
   !theta = MATMUL(MATMUL(MATMUL(invlambda,TRANSPOSE(y2)),da_dir),x1)
   theta = MATMUL( invlambda, MATMUL( MATMUL( TRANSPOSE(y2), da_dir), x1 ) )

   ! dy = y * [ zeros(nz,nz) -theta'; theta zeros(nnz,nnz) ]
   themat(:,:) = 0.0d0
   themat( 1:nz, nz+1:neq ) = -1.0d0 * TRANSPOSE(theta)
   themat( nz+1:neq, 1:nz ) = theta

   dy(i,:,:) = MATMUL( y, themat )
END DO

DEALLOCATE( y1, y2, x1, x2, invlambda, da_dir, theta, themat, STAT=status )
IF( status .NE. 0 ) THEN
   WRITE(*,*) 'ERROR in svd_gary_analytical_deriv(): Could not deallocate memory'
   RETURN
END IF
END SUBROUTINE
