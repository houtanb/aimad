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

SUBROUTINE EIGENSYSTEM_D(asmall, asmalld, ia, uprbnd, w, wd, nlarge, &
&  info)


IMPLICIT NONE
INTEGER, INTENT(IN) :: ia
DOUBLE PRECISION, INTENT(IN) :: uprbnd
DOUBLE PRECISION, DIMENSION(ia,ia), INTENT(IN) :: asmall
DOUBLE PRECISION, DIMENSION(ia, ia), INTENT(IN) :: asmalld
DOUBLE PRECISION, DIMENSION(ia,ia), INTENT(OUT) :: w
DOUBLE PRECISION, DIMENSION(ia, ia), INTENT(OUT) :: wd
INTEGER, INTENT(OUT) :: nlarge, info

! variables needed for the call to dgeev
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: acopy
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: wr, wi
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: work, iwork
LOGICAL, ALLOCATABLE, DIMENSION(:) :: bwork
INTEGER :: lwork, liwork, nsmall
DOUBLE PRECISION :: rcondv,rconde


! local variables
INTEGER :: indxi, indxj, status

LOGICAL, EXTERNAL :: selectf
EXTERNAL :: dgeesx, transp, printmat, eigensystem_analytd


lwork = ia*ia
liwork = ia*ia
ALLOCATE(acopy(ia,ia), wr(ia),wi(ia),work(lwork),iwork(liwork),bwork(ia),STAT=status)

acopy = asmall


CALL dgeesx('V','S',selectf,'B',ia,acopy,ia,nsmall,wr,wi,w,ia,rconde,rcondv,work,lwork,iwork,liwork,bwork,info)

nlarge = ia -nsmall

CALL eigensystem_analytd(asmalld,ia,w,acopy,nsmall,wd,info)


DEALLOCATE(acopy,wr,wi,work,iwork,bwork)


END SUBROUTINE
