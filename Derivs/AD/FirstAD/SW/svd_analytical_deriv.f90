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

SUBROUTINE svd_analytical_deriv( a, da_3d, neq, u, du_3d, d, v, nbdirs )
USE DIFFSIZES
IMPLICIT NONE
INTEGER, INTENT(IN) :: neq, nbdirs
DOUBLE PRECISION, INTENT(IN), DIMENSION(neq,neq) :: a, u, v
DOUBLE PRECISION, INTENT(IN), DIMENSION(nbdirsmax,neq,neq) :: da_3d
DOUBLE PRECISION, INTENT(IN), DIMENSION(neq) :: d
DOUBLE PRECISION, INTENT(OUT), DIMENSION(nbdirsmax,neq,neq) :: du_3d

! Local variables
INTEGER :: i, j, k, l, flag, status, dindx
DOUBLE PRECISION :: zerotol, det
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: omegau, rhs, singvalmat, inverse, da, du
EXTERNAL :: pseudoinv

!
! Zero out deriv matrices
!
du_3d(:,:,:) = 0.0d0


!
! Allocate memory
!
ALLOCATE( omegau(neq,neq), rhs(2,1), singvalmat(2,2), inverse(2,2), da(neq,neq), du(neq,neq), STAT=status )
IF (status .NE. 0) THEN
   WRITE(*,*)
   WRITE(*,*) '**** IN SVD_ANALYTICAL_DERIV(): COULD NOT DEALLOCATE MEMORY ****'
   RETURN
END IF


DO dindx=1,nbdirs
flag = 0
du(:,:) = 0.0d0
da(:,:) = 0.0d0

da(:,:) = da_3d(dindx,:,:)

!
! Begin calculation of jacobian from:
! Papadopoulo & Lourakis (2000), Equation (8)
!
zerotol = 1e-8
DO j=1,neq
   DO i=1,neq

      IF( ABS(da(i,j)) > zerotol ) THEN
         omegau(:,:) = 0.0d0

         DO l=1,neq
            DO k=l+1,neq

               rhs(1,1) =  u(i,k)*v(j,l)
               rhs(2,1) = -u(i,l)*v(j,k)

               IF( d(l) .EQ. d(k) ) THEN
                  IF( flag .NE. 0 ) THEN
                     !calculate pseudoinverse
                     singvalmat(1,1) = d(l) !a
                     singvalmat(1,2) = d(k) !b
                     singvalmat(2,1) = d(k) !c
                     singvalmat(2,2) = d(l) !d

                     CALL pseudoinv( singvalmat, 2, inverse )
                     flag = 1
                  END IF
               ELSE
                  det = d(l)*d(l) - d(k)*d(k)

                  IF( abs(det) > zerotol ) THEN
                     !calculate inverse
                     inverse(1,1) = d(l)/det
                     inverse(1,2) = -d(k)/det
                     inverse(2,1) = -d(k)/det
                     inverse(2,2) = d(l)/det
                  ELSE
                     !calculate pseudoinverse
                     singvalmat(1,1) = d(l) !a
                     singvalmat(1,2) = d(k) !b
                     singvalmat(2,1) = d(k) !c
                     singvalmat(2,2) = d(l) !d

                     CALL pseudoinv( singvalmat, 2, inverse )
                  END IF
                  flag = 0
               END IF
               omegau(k,l) = inverse(1,1) * rhs(1,1) + inverse(1,2) * rhs(2,1)
               omegau(l,k) = -omegau(k,l)
            END DO
         END DO
         du = du + da(i,j)*MATMUL(u,omegau)
      END IF
   END DO
END DO

du_3d(dindx,:,:) = du(:,:)

END DO


!
! Deallocate memory
!
DEALLOCATE( omegau, rhs, singvalmat, inverse, da, du, STAT=status )
IF (status .NE. 0) THEN
   WRITE(*,*)
   WRITE(*,*) '**** IN SVD_ANALYTICAL_DERIV(): COULD NOT DEALLOCATE MEMORY ****'
   RETURN
END IF

END SUBROUTINE
