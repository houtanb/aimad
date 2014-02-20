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

SUBROUTINE objfun_likel_2pt_fd(n, x, objf, objgrd, iuser, user )

USE modminimizer

IMPLICIT NONE
INTEGER, INTENT(IN) :: n
DOUBLE PRECISION, INTENT(IN), DIMENSION(n) :: x
INTEGER, INTENT(INOUT) :: iuser(*)
DOUBLE PRECISION, INTENT(INOUT) :: user(*)
DOUBLE PRECISION, INTENT(OUT) :: objf
DOUBLE PRECISION, INTENT(OUT), DIMENSION(n) :: objgrd
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: x_new
DOUBLE PRECISION :: psitemp, objfone, objftwo
INTEGER :: i
EXTERNAL GET_HMAT, likelihood

objfuncalls = objfuncalls + 1

ALLOCATE( x_new(n), STAT=status )
IF (status .NE. 0) THEN
   WRITE(*, *)
   WRITE(*, *) '***In objfun_likel: Could not allocate memory***'
END IF
x_new = MIN( MAX( x, lwrbndx ), uprbndx )

DO i=1,n
   psi(iuser(i)) = x_new(i)
END DO

objf = 0.0d0
objgrd(:) = 0.0d0

DO i=1,n
   psitemp = psi(iuser(i))

   psi(iuser(i)) = psitemp + theh

   h(:,:) = 0.0d0
   CALL GET_HMAT( h, hrows, hcols, psi, psin )

   !CALL writemat(h,hrows,hcols,1,hrows,1,hcols,"d:/fortran/sigma2/hmat.txt")
   !CALL writemat_screen(h,hrows,hcols,1,hrows,1,hcols)
   CALL likelihood( objfone, &
        endog, nendogvars, &
        errlist, nerrvars, &
        observed, nobservedvars, &
        zdata, tt, &
        h, hrows, hcols, &
        nlead, nlag, neq, ntrain )
   nendogvars = nendogvarsbak
   endog = endogbak


   h(:,:) = 0.0d0
   psi(iuser(i)) = psitemp - theh
   CALL GET_HMAT( h, hrows, hcols, psi, psin )
   CALL likelihood( objftwo, &
        endog, nendogvars, &
        errlist, nerrvars, &
        observed, nobservedvars, &
        zdata, tt, &
        h, hrows, hcols, &
        nlead, nlag, neq, ntrain )
   nendogvars = nendogvarsbak
   endog = endogbak
   objgrd(i) = (objfone - objftwo)/(2.0d0 * theh)

   psi(iuser(i)) = psitemp

   IF( .NOT.(x_new(i) .EQ. x(i)) ) THEN
      objgrd(i) = 100000000.0d0 * ( x(i) - x_new(i) )/ABS( x(i) - x_new(i) )
   END IF
END DO



h(:,:) = 0.0d0
CALL GET_HMAT( h, hrows, hcols, psi, psin )
CALL likelihood( objf, &
                 endog, nendogvars, &
                 errlist, nerrvars, &
                 observed, nobservedvars, &
                 zdata, tt, &
                 h, hrows, hcols, &
                 nlead, nlag, neq, ntrain )

nendogvars = nendogvarsbak
endog = endogbak

objf = objf + 100000000.0d0 * SUM( ABS( x_new - x ) )

DEALLOCATE( x_new, STAT=status )
IF (status .NE. 0) THEN
   WRITE(*, *)
   WRITE(*, *) '***In objfun_likel: Could not deallocate memory***'
END IF
END SUBROUTINE
