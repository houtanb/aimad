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

SUBROUTINE kalman( f, r, h, n, y, t, xi10, q, ntrain, likel,explswitch )

USE moderror

IMPLICIT NONE

INTEGER, INTENT(IN) :: r, n, t, ntrain
INTEGER, INTENT(IN) :: explswitch ! if set to 1, imposes a penalty if likelihood explodes
                                  ! otherwise, skips the update for any observation that would cause explosion
DOUBLE PRECISION, INTENT(IN), DIMENSION(r,r) :: f
DOUBLE PRECISION, INTENT(IN), DIMENSION(r,n) :: h
DOUBLE PRECISION, INTENT(IN), DIMENSION(n,t) :: y
DOUBLE PRECISION, INTENT(IN), DIMENSION(r,1) :: xi10
DOUBLE PRECISION, INTENT(IN), DIMENSION(r,r) :: q
DOUBLE PRECISION, INTENT(OUT) :: likel

! Local Vars
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: xi10history, xi11history, p10, hp
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: part1, tmprbyn, hpp10h, b, tmprbyn_
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: tmprby1, tmprbyr, tmprbyr_, p11, tmp1byn
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: part2, origa
DOUBLE PRECISION :: pi, nhalflog2pi, e, detl, detu, det_hpp10h, likelhd
INTEGER, ALLOCATABLE, DIMENSION(:) :: ipiv
INTEGER :: i, j, status, info

INTEGER :: likelcount
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: likeldata

EXTERNAL matrixmult, eyemat, transp, dgesv, dgemm

!CALL writemat(f,r,r,1,r,1,r,"kfmat.txt")

!CALL writemat(h,r,n,1,r,1,n,"khmat.txt")


!CALL writemat(q,r,r,1,r,1,r,"kqmat.txt")

! if subroutine successfully completes task, info is set to 0
! info is set to -1 if memory allocation was a problem
! info is set to -2 if likelihood exploded
! info is set to -3 if there were problem with inverse computation
! info is set to -4 if more than one of the above occurred

!
! Begin Processing
!
e       = 2.718281828459046d0
pi      = 3.141592653589793d0
nhalflog2pi = -DBLE(n)/2.0d0*LOG(2.0d0*pi)
info    = 0
memerr  = 0
inverr  = 0
explerr = 0
likel   = 0.0d0
likelhd = 0.0d0

ALLOCATE(likeldata(t,1))
likelcount = 0

ALLOCATE( xi10history(r,1), xi11history(r,1), p10(r,r), hp(n,r), &
          tmprbyn(r,n), hpp10h(n,n), origa(n,n), part1(n,1), ipiv(n), b(n,n), &
          tmprbyn_(r,n), tmprby1(r,1), tmprbyr(r,r), tmprbyr_(r,r), &
          p11(r,r), tmp1byn(1,n), part2(1,1), STAT=status )
IF( status .NE. 0) THEN
   WRITE(*,*)
   WRITE(*,*) "***In Kalman Filter: Could not allocate memory***"
   memerr = -1
END IF


!
! begin kalman filter, (references to Hamilton)
!
! 13.2.15
xi11history = xi10

! 13.2.17
CALL matrixmult( f, r, r, xi11history, 1, xi10history )

! 13.2.21
p10 = q

DO i=2,t
   CALL transp( h, r, n, hp )

   ! part1 = y(:, i) - h'*xi10history
   CALL matrixmult( hp, n, r, xi10history, 1, part1 )
   part1(:,1) = y(:,i) - part1(:,1)

   ! b = (h'*p10*h)^-1
   CALL matrixmult( p10, r, r, h, n, tmprbyn )
   CALL matrixmult( hp, n, r, tmprbyn, n, hpp10h )
   CALL eyemat(b,n)
   origa = hpp10h
   CALL dgesv( n, n, hpp10h, n, ipiv, b, n, info)
   CALL invert( origa, n )
   b = origa

   IF ( info .NE. 0 ) THEN
!      WRITE(*,*)
!      WRITE(*,*) "***In Kalman Filter: Error on return from dgesv()***"
      inverr = -1
   END IF

   IF ( i > ntrain ) THEN
      !
      ! Note: do not assign to tmprbyn in this if statement
      !
      ! part2 = -0.5 * part1' * inv( hpp10h ) * part1
      tmp1byn(:,:) = 0.0d0
      CALL dgemm( 'T', 'N', 1, n, n, -0.5d0, part1, n, b, n, 1, tmp1byn, 1 );
      CALL matrixmult( tmp1byn, 1, n, part1, 1, part2 )

      ! 13.4.1
      ! det_hpp10h = det( hpp10h )
      detl = 1.0d0
      detu = 1.0d0
      DO j=1,n
         IF ( ipiv(j) .NE. j ) THEN
            detl = detl * (-1.0d0)
         END IF
         detu = detu * hpp10h(j,j)
      END DO
      det_hpp10h = detl * detu
      !write(*,*) 'det_hpp10h'
      !write(*,*) det_hpp10h
      ! 13.4.1 (continued)
      ! likelhd = 2PI^-n/2 * det_hpp10h^-0.5 * e^part2
      IF (det_hpp10h > 10.0d0**(-300)) THEN

         likelhd = nhalflog2pi -0.5d0*LOG(det_hpp10h) + part2(1,1)
         likelcount = likelcount+1
         likeldata(likelcount,1) = likelhd
         !WRITE(*,*) likelhd

         IF(likel+likelhd > 10.0d0**300 ) THEN
            IF (explswitch == 1) THEN
               likel = 10.0d0**300
            END IF
            explerr = -1
         ELSE
            likel = likel+likelhd
			!write(*,*)  'likel'
			!write(*,*) likel
         END IF

      ELSE
         IF (explswitch == 1) THEN
            likel = 10.0d0**300
         END IF
         explerr = -1
      END IF
   END IF


   ! 13.2.15
   ! xi11history = xi10history + p10*h*b*part1
   CALL matrixmult( tmprbyn, r, n, b, n, tmprbyn_ )
   CALL matrixmult( tmprbyn_, r, n, part1, 1, tmprby1 )
   xi11history = xi10history + tmprby1

   ! 13.2.16
   ! p11 = p10 - p10*h*b*h'*p10
   CALL matrixmult( tmprbyn_, r, n, hp, r, tmprbyr_ )
   CALL matrixmult( tmprbyr_, r, r, p10, r, tmprbyr )
   p11 = p10 - tmprbyr

   ! 13.2.17
   ! xi10history = f * xi11history
   CALL matrixmult( f, r, r, xi11history, 1, xi10history )

   ! 13.2.21
   ! p10 = f*p11*f' + q
   p10 = q
   CALL matrixmult( f, r, r, p11, r, tmprbyr )
   CALL dgemm( 'N', 'T', r, r, r, 1.0d0, tmprbyr, r, f, r, 1.0d0, p10, r )
   !write(*,*) 'p10'
   !call writemat(p10,r,r,1,r,1,r,'p10.txt')
END DO

IF ((explerr .NE. -1) .AND. (explswitch .EQ. 1)) THEN
   likel = -1.0d0 * likel
END IF

!CALL writemat(likeldata,t,1,1,likelcount,1,1,'likeldata.txt')
DEALLOCATE(likeldata)
DEALLOCATE( xi10history, xi11history, p10, hp, &
          tmprbyn, hpp10h, part1, ipiv, b, &
          tmprbyn_, tmprby1, tmprbyr, tmprbyr_, &
          p11, tmp1byn, part2, origa, STAT=status )
IF( status .NE. 0) THEN
   WRITE(*,*)
   WRITE(*,*) "***In Kalman Filter: Could not deallocate memory***"
   memerr = -1
END IF


END SUBROUTINE
