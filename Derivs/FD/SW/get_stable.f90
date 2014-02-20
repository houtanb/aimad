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

SUBROUTINE get_stable(a,nrows,condn,stablepos,ns)

IMPLICIT NONE

INTEGER, INTENT(IN) :: nrows
DOUBLE PRECISION, DIMENSION(nrows,nrows), INTENT(IN) :: a
DOUBLE PRECISION, INTENT(IN) :: condn
INTEGER, INTENT(OUT) :: ns
LOGICAL, DIMENSION(nrows), INTENT(OUT) :: stablepos

INTEGER :: indxi
ns = 0
DO indxi=1,nrows
   IF(ABS(a(indxi,indxi)) .LE. condn) THEN
      stablepos(indxi)=.TRUE.
      ns=ns+1
   ELSE
      stablepos(indxi)=.FALSE.
   END IF
END DO

END SUBROUTINE
