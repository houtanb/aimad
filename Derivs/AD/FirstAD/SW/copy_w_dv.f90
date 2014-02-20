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

!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade - Version 2.2 (r1239) - Wed 28 Jun 2006 04:59:55 PM CEST
!
!  Differentiation of copy_w in forward (tangent) mode: (multi-directional mode)
!   variations  of output variables: q
!   with respect to input variables: q w
SUBROUTINE COPY_W_DV(q, qd, qrows, qcols, w, wd, ia, iq, js, nbdirs)
  USE DIFFSIZES
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ia
  INTEGER, INTENT(IN) :: iq
  INTEGER, DIMENSION(ia), INTENT(IN) :: js
  INTEGER :: nbdirs
  INTEGER, INTENT(IN) :: qcols
  INTEGER, INTENT(IN) :: qrows
  DOUBLE PRECISION, DIMENSION(qrows, qcols), INTENT(INOUT) :: q
  DOUBLE PRECISION :: qd(nbdirsmax, qrows, qcols)
  DOUBLE PRECISION, DIMENSION(ia, ia), INTENT(IN) :: w
  DOUBLE PRECISION, DIMENSION(nbdirsmax, ia, ia), INTENT(IN) :: wd
  INTEGER :: nd, status
  DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: wt
  DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: wtd

  ALLOCATE( wt(qrows-iq, ia), wtd(nbdirsmax,qrows-iq, ia), STAT=status )
  IF (status .NE. 0) THEN
     WRITE(*,*) '**** IN COPY_W_DV():: COULD NOT ALLOCATE MEMORY ****'
     RETURN
  END IF


  IF (iq .LT. qrows) THEN
    wtd(:,:,:) = 0.0d0
    CALL TRANSP_DV(w(:, ia-(qrows-iq)+1:ia), wd(:, :, ia-(qrows-iq)+1:ia), ia, qrows - iq, wt, wtd, nbdirs)
    DO nd=1,nbdirs
      qd(nd, iq+1:qrows, js) = wtd(nd, :, :)
    END DO
    q(iq+1:qrows, js) = wt
  END IF

  DEALLOCATE( wt, wtd, STAT=status )
  IF (status .NE. 0) THEN
     WRITE(*,*) '**** IN COPY_W_DV():: COULD NOT DEALLOCATE MEMORY ****'
     RETURN
  END IF
END SUBROUTINE COPY_W_DV