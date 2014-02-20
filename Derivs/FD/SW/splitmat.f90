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

SUBROUTINE splitmat( cofb, endog, nendogvars, errlist, nerrvars, amat, bmat )
!!
!! given law of motion of the form
!! [x_t]' = cofb [x_{t-1}]'
!! this function produces
!! amat and bmat such that
!! y_t = amat y_{t-1} + bmat e_t
!!
!! the names of the variables in x_t are listed in endog
!! the names of the variables in e_t are listed in errlist
!!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nerrvars
INTEGER, INTENT(INOUT) :: nendogvars
CHARACTER(len=*), DIMENSION(nendogvars), INTENT(INOUT) :: endog
CHARACTER(len=*), DIMENSION(nerrvars), INTENT(IN) :: errlist
DOUBLE PRECISION, INTENT(IN), DIMENSION(nendogvars, nendogvars) :: cofb
DOUBLE PRECISION, INTENT(OUT), DIMENSION(nendogvars-nerrvars,nendogvars-nerrvars) :: amat
DOUBLE PRECISION, INTENT(OUT), DIMENSION(nendogvars-nerrvars,nerrvars) :: bmat

!LOCAL VARS
INTEGER, DIMENSION(nendogvars-nerrvars) :: endogposlist
INTEGER, DIMENSION(nerrvars) :: errposlist
INTEGER :: i, endogindx, errindx
INTEGER, EXTERNAL :: strmatch

errindx         = 0
errposlist(:)   = 0
endogindx       = 0
endogposlist(:) = 0
DO i=1,nendogvars
   IF( strmatch( endog(i), errlist, nerrvars ) < 0 ) THEN
      endogindx   = endogindx + 1
      endogposlist(endogindx) = i
   ELSE
      errindx   = errindx + 1
      errposlist(errindx) = i
   END IF
END DO

amat(:,:) = 0.0d0
bmat(:,:) = 0.0d0
amat  = cofb(endogposlist,endogposlist)
bmat  = cofb(endogposlist,errposlist)

DO i=1,nendogvars
   IF ( i > nendogvars-nerrvars ) THEN
      endog(i) = ''
   ELSE
      endog(i) = endog(endogposlist(i))
   END IF
END DO
nendogvars = nendogvars - nerrvars

END SUBROUTINE
