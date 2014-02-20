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

SUBROUTINE get_selector(zerorows,neq,selector,nz)

! turn the logical zerorows array into a selector
! that holds the positions where zerorows is true

IMPLICIT NONE
INTEGER, INTENT(IN) :: neq,nz
LOGICAL, DIMENSION(neq), INTENT(IN) :: zerorows
INTEGER, DIMENSION(nz), INTENT(OUT) :: selector

INTEGER :: indxi, indxk


indxk = 1

DO indxi=1,neq
   IF (zerorows(indxi)) THEN
      selector(indxk) = indxi
      indxk=indxk+1
   END IF
END DO

END SUBROUTINE
