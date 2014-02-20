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

SUBROUTINE eigensystem(asmall,ia,uprbnd,w,nlarge,info)

! computes eigenvalues and eigenvectors of asmall
IMPLICIT NONE
INTEGER, INTENT(IN) :: ia
DOUBLE PRECISION, INTENT(IN) :: uprbnd
DOUBLE PRECISION, DIMENSION(ia,ia), INTENT(IN) :: asmall
DOUBLE PRECISION, DIMENSION(ia,ia), INTENT(OUT) :: w
INTEGER, INTENT(OUT) :: nlarge, info

! variables needed for the call to dgeev
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: acopy
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: wr, wi
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: work, iwork
LOGICAL, ALLOCATABLE, DIMENSION(:) :: bwork
INTEGER :: lwork, liwork, nsmall
DOUBLE PRECISION :: rcondv,rconde


! local variables
INTEGER :: status

LOGICAL, EXTERNAL :: selectf
EXTERNAL :: dgeesx


lwork = ia*ia
liwork = ia*ia
ALLOCATE(acopy(ia,ia), wr(ia),wi(ia),work(lwork),iwork(liwork),bwork(ia),STAT=status)

acopy = asmall
CALL dgeesx('V','S',selectf,'B',ia,acopy,ia,nsmall,wr,wi,w,ia,rconde,rcondv,work,lwork,iwork,liwork,bwork,info)
nlarge = ia - nsmall

DEALLOCATE( acopy, wr, wi, work, iwork, bwork, STAT=status )

END SUBROUTINE
