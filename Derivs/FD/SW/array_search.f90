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

FUNCTION array_search(match,int_array,size)
! search int_array for match
! if int_array contains match, then return position of first occurence
! otherwise return zero

IMPLICIT NONE

INTEGER, INTENT(IN) :: match, size
INTEGER, DIMENSION(size), INTENT(IN) :: int_array

INTEGER :: array_search

INTEGER :: indxi

array_search = 0

DO indxi = 1,size

   IF (int_array(indxi)==match) THEN
      array_search = indxi
      CONTINUE
   END IF

END DO



END FUNCTION
