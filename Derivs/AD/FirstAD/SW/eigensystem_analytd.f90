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

SUBROUTINE eigensystem_analytd(ad,ia,x,t,nsmall,xd,info)

! implements differentiation formula in section 4.2 of
! Anderson, "A procedure for differentiating perfect-foresight-model
!            reduced-form coefficients", in JEDC (1987).

IMPLICIT NONE

INTEGER, INTENT(IN) :: ia, nsmall
INTEGER, INTENT(OUT) :: info
DOUBLE PRECISION, DIMENSION(ia, ia), INTENT(IN) :: ad, t, x
DOUBLE PRECISION, DIMENSION(ia, ia), INTENT(OUT) :: xd


DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: x1, x2, e11, e22, eye_11, eye_22
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: part1_1, part1_2, bigpi
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: part2

INTEGER, ALLOCATABLE, DIMENSION(:) :: ipiv

INTEGER :: status, npart1

EXTERNAL kroneckerproduct, constructpimat, dgesv, eyematnod

ALLOCATE(x1(ia,nsmall),x2(ia,ia-nsmall),e11(nsmall,nsmall),e22(ia-nsmall,ia-nsmall), &
         eye_11(ia-nsmall,ia-nsmall),eye_22(nsmall,nsmall), &
         part1_1((ia-nsmall)*nsmall,(ia-nsmall)*nsmall), part1_2((ia-nsmall)*nsmall,(ia-nsmall)*nsmall), &
         part2((ia-nsmall)*nsmall), ipiv((ia-nsmall)*nsmall), &
         bigpi(ia,ia), STAT=status)
IF( status .NE. 0 ) THEN
   WRITE(*,*)
   WRITE(*,*) "***In EIGENSYSTEM_ANALYTD: Could not allocate memory***"
END IF

! define partitions
x1 = x(:,1:nsmall)
x2 = x(:,nsmall+1:ia)
e11 = t(1:nsmall,1:nsmall)
e22 = t(nsmall+1:ia,nsmall+1:ia)


! define the identity matrices
CALL eyematnod(eye_11,ia-nsmall)
CALL eyematnod(eye_22,nsmall)


! define part1
! part1 = kron(e11',eye_11) - kron(eye_22,e22)
! part1_1 = kron(e11',eye_11)
! part1_2 = kron(eye_22,e22)


CALL kroneckerproduct(TRANSPOSE(e11),nsmall,nsmall,eye_11,ia-nsmall,ia-nsmall,part1_1)
CALL kroneckerproduct(eye_22,nsmall,nsmall,e22,ia-nsmall,ia-nsmall,part1_2)

part1_1 = part1_1 -part1_2

! define part2
! part2 = vec(x2'*ad*x1)
part2 = RESHAPE( MATMUL(TRANSPOSE(x2),MATMUL(ad,x1) )  ,(/(ia-nsmall)*nsmall/) )

npart1 =(ia-nsmall)*nsmall


! construct part1^(-1) * part2, result stored in part2
CALL dgesv(npart1,1,part1_1,npart1,ipiv,part2,npart1,info)



! construct pi
CALL constructpimat(RESHAPE(part2,(/(ia-nsmall), nsmall/)),ia-nsmall,nsmall,bigpi)


! get final product xd=x*pi


xd = MATMUL(x,bigpi)


DEALLOCATE(x1, x2, e11, e22, eye_11, eye_22, part1_1, &
     part1_2, part2, ipiv, bigpi, STAT=status)
IF( status .NE. 0 ) THEN
   WRITE(*,*)
   WRITE(*,*) "***In EIGENSYSTEM_ANALYTD: Could not deallocate memory***"
END IF

END SUBROUTINE


