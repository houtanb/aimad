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

SUBROUTINE eyemat(eye_n,n)

IMPLICIT NONE

INTEGER, INTENT(IN) :: n

DOUBLE PRECISION, DIMENSION(n,n), INTENT(OUT) :: eye_n

INTEGER indxi, indxj

eye_n(:,:) = 0.0d0

DO indxi=1,n
   eye_n(indxi,indxi) = 1.0d0
END DO


END SUBROUTINE
