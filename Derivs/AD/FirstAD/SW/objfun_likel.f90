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

SUBROUTINE objfun_likel(mode, n, x, objf, objgrd, nstate, iuser, user )

USE modminimizer

IMPLICIT NONE
INTEGER, INTENT(IN) :: n
DOUBLE PRECISION, INTENT(IN), DIMENSION(n) :: x
INTEGER, INTENT(INOUT) :: mode, nstate
INTEGER, INTENT(INOUT) :: iuser(*)
DOUBLE PRECISION, INTENT(INOUT) :: user(*)
DOUBLE PRECISION, INTENT(OUT) :: objf
DOUBLE PRECISION, INTENT(OUT), DIMENSION(n) :: objgrd
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: x_new
INTEGER :: i
EXTERNAL GET_HMAT_DV, likelihood_dv
INTEGER, EXTERNAL :: MCLOCK

objfuncalls = objfuncalls + 1

ALLOCATE( x_new(n), STAT=status )
IF (status .NE. 0) THEN
   WRITE(*, *)
   WRITE(*, *) '***In objfun_likel: Could not allocate memeory***'
END IF
x_new = MIN( MAX( x, lwrbndx ), uprbndx )
DO i=1,n
   psi(iuser(i)) = x_new(i)
END DO

objf = 0.0d0
objgrd(:) = 0.0d0

h(:,:) = 0.0d0
hd(:,:,:) = 0.0d0

psid(:,:) = 0.0d0
DO i=1,nbdirs
   psid(i,iuser(i)) = 1.0d0
END DO

CALL GET_HMAT_DV( h, hd, hrows, hcols, psi, psid, psin, nbdirs )

call writemat(h,hrows,hcols,1,hrows,1,hcols,"hmat_dv.m", LEN("hmat_dv.m"))
call writemat_total_3d(hd,nbdirs,hrows,hcols,'hmatd_dv.m', LEN('hmatd_dv.m'))

WRITE(*,*) 'LOADED h, ABOUT TO CALL LIKELIHOOD'
CALL likelihood_dv( objf, objgrd, &
     endog, nendogvars, &
     errlist, nerrvars, &
     observed, nobservedvars, &
     zdata, tt, &
     h, hd, hrows, hcols, &
     nlead, nlag, neq, ntrain, nbdirs )

nendogvars = nendogvarsbak
endog = endogbak

DO i=1,nbdirs
   IF( .NOT.(x_new(i) .EQ. x(i)) ) THEN
      objgrd(i) = 100000000.0d0 * ( x(i) - x_new(i)  )/ABS( x(i) - x_new(i)  )
   END IF
END DO

objf = objf + 100000000.0d0 * SUM( ABS( x_new - x ) )

DEALLOCATE( x_new, STAT=status )
IF (status .NE. 0) THEN
   WRITE(*, *)
   WRITE(*, *) '***In objfun_likel: Could not deallocate memeory***'
END IF
END SUBROUTINE
