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

SUBROUTINE shiftright(x,nrows,ncols,nshift,newx)

IMPLICIT NONE

INTEGER, INTENT(IN) :: nrows,ncols,nshift
DOUBLE PRECISION, DIMENSION(nrows,ncols), INTENT(IN) :: x
DOUBLE PRECISION, DIMENSION(nrows,ncols), INTENT(OUT) :: newx

INTEGER :: lleft,rleft
INTEGER :: lright,rright

INTEGER :: indxi,indxj

lleft = 1
rleft = ncols-nshift
lright = nshift+1
rright = ncols

newx(:,lright:rright) = x(:,lleft:rleft)

DO indxj=1,nshift
   DO indxi=1,nrows
      newx(indxi,indxj) = 0.d0
   END DO
END DO

END SUBROUTINE
