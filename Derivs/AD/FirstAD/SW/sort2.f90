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

SUBROUTINE SSORT (X, IX, n)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: n

  DOUBLE PRECISION, DIMENSION(n), INTENT(INOUT) :: X
  INTEGER, DIMENSION(n), INTENT(INOUT) ::  IX

  DOUBLE PRECISION :: tmp
  INTEGER i, itmp, tt

  LOGICAL :: changed

  i = 0
  tt = 0
  changed = .FALSE.
  DO WHILE (.TRUE.)
     i = i+1
     IF (i.EQ.n) THEN
        IF (changed.EQV..FALSE.) THEN
           WRITE(*,*) '!!!!!!!!!!!!!!!!'
           WRITE(*,*) tt
           RETURN
        END IF
        i = 0
        tt = tt + 1
        changed = .FALSE.
     ELSE
        IF (X(i) .GT. 4) THEN
           WRITE(*,*) X(i)
           WRITE(*,*) X(i+1)
        END IF
        IF( X(i).GT.X(i+1)) THEN
           changed = .TRUE.
           tmp    = X(i)
           X(i)   = X(i+1)
           X(i+1) = tmp

           itmp    = IX(i)
           IX(i)   = IX(i+1)
           IX(i+1) = itmp
        END IF
     END IF
  END DO
END SUBROUTINE SSORT
