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

SUBROUTINE shrinkspacesmall(amat,bmat,nendog,nerrvars,observedpos,nobs,selector,nselected)
!CALL shrinkspacesmall(amat,bmat,nendogvarscopy,nerrvars,observedpos,nobservedvars,selector,nselected)

IMPLICIT NONE

INTEGER, INTENT(IN) :: nendog, nerrvars, nobs
DOUBLE PRECISION, DIMENSION(nendog,nendog), INTENT(IN) :: amat
DOUBLE PRECISION, DIMENSION(nendog,nerrvars), INTENT(IN) :: bmat
INTEGER, DIMENSION(nobs), INTENT(IN) :: observedpos
INTEGER, INTENT(OUT) :: nselected
INTEGER, DIMENSION(nendog), INTENT(OUT) :: selector

DOUBLE PRECISION :: tolzero = 0.0000000001
INTEGER :: indxi
INTEGER, EXTERNAL :: array_search



nselected = 0
DO indxi = 1, nendog
     IF ( MAXVAL(ABS(amat(:,indxi)))>tolzero .OR. array_search(indxi,observedpos,nobs)>0 ) THEN
         nselected=nselected+1
         selector(nselected) = indxi
     END IF
END DO




END SUBROUTINE
