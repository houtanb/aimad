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

FUNCTION strmatch( str1, strarr, ldaarr )

IMPLICIT NONE
INTEGER :: strmatch
INTEGER, INTENT(IN) :: ldaarr
CHARACTER(len=*), INTENT(IN) :: str1
CHARACTER(len=*), DIMENSION(ldaarr), INTENT(IN) :: strarr

INTEGER :: i, lenstr1

strmatch = -1
lenstr1  = LEN_TRIM(str1)
DO i=1,ldaarr
   IF( (INDEX(str1, strarr(i))>0) .OR. (INDEX(strarr(i), str1)>0) ) THEN
      IF( LEN_TRIM(strarr(i)) == lenstr1 ) THEN
         strmatch = i
         RETURN
      END IF
   END IF
END DO
END FUNCTION
