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

SUBROUTINE existstrue( mat, n, ans )

IMPLICIT NONE
INTEGER, INTENT(IN) :: n
LOGICAL, INTENT(IN), DIMENSION(n) :: mat
LOGICAL, INTENT(OUT) :: ans

INTEGER :: i

ans = .FALSE.
DO i=1,n
   IF( mat(i) ) THEN
      ans = .TRUE.
   END IF
END DO

END SUBROUTINE
