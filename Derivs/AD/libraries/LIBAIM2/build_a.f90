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

SUBROUTINE build_a(h,qcols,neq,a,ia,js)

IMPLICIT NONE

! variables in calling structure
INTEGER, INTENT(IN) :: qcols, neq
DOUBLE PRECISION, DIMENSION(neq,3*neq), INTENT(IN) :: h
DOUBLE PRECISION, DIMENSION(qcols,qcols), INTENT(OUT) :: a
INTEGER, INTENT(OUT) :: ia   ! Effective dimension of a AND js
INTEGER, DIMENSION(qcols), INTENT(OUT) :: js

! variables needed for call to get_zerocols
INTEGER :: nz
DOUBLE PRECISION :: zerotol = 0.000000001d0
LOGICAL, ALLOCATABLE, DIMENSION(:) :: zerocols


! local variables
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: ah
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: bh, c
INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv
INTEGER :: lleft, rleft
INTEGER :: lright, rright
INTEGER :: indxi, indxj, newnz
INTEGER :: status, info

EXTERNAL :: dgesv, get_zerocols

lleft = 1
rleft = qcols
lright = qcols+1
rright = qcols+neq

! construct hout
!
! hout = -h(:,right)\h(:,left) = bh\ah
! to construct this

! get the qr factorization of ah so that bh*P = Q*R
! then compute hout as P*(R\(Q'*ah))
! use DORMQR to construct Q'*ah
! use DTRTRS to construct R\(Q'*ah)
! finally get hout as P*((R\(Q'*ah))

! construct appropriate partitions of h
info = 0
ALLOCATE(ah(neq,qcols),bh(neq,neq),ipiv(neq), c(neq, qcols),STAT=status)
IF (status .NE. 0) THEN
   !WRITE(*,*) '**** COULD NOT ALLOCATE MEMORY ****'
   !RETURN
END IF
ah = -h(:, lleft:rleft)
bh = h(:, lright:rright)

!CALL dgesv(neq, qcols, bh, neq, ipiv, ah, neq, info)
CALL invert( bh, neq )
c(:,:) = 0.0d0
CALL dgemm( 'N', 'N', neq, qcols, neq, 1.0d0, bh, neq, ah, neq, 0.0d0, c, neq )
ah = c


! build the big transition matrix
a(:,:) = 0.0d0
DO indxi = 1,neq
   a(indxi,neq+indxi) = 1.d0
END DO

a(neq+1:qcols,:) = ah


!  Delete inessential lags and build index array js.  js indexes the
!  columns in the big transition matrix that correspond to the
!  essential lags in the model.  They are the columns of q that will
!  get the unstable left eigenvectors.

ALLOCATE(zerocols(qcols),STAT=status)
DO indxi = 1,qcols
   zerocols(indxi) = .FALSE.
END DO

nz = 0  ! holds number of identified zero columns
newnz = nz
!WRITE(*,*) 'a'
!CALL printmat(a,qcols,qcols)

CALL get_zerocols(a,qcols,qcols,zerotol,zerocols,newnz)

DO
   IF (newnz>nz .AND. newnz<qcols) THEN
      nz = newnz
      CALL get_zerocols(a,qcols,qcols,zerotol,zerocols,newnz)
   ELSE
      EXIT
   END IF
END DO

! define js so that the first entries hold the position
! of the essential lags (there will be ia = qcols-nz of them)
indxj = 1
DO indxi = 1,qcols
   IF (.NOT. zerocols(indxi)) THEN
      js(indxj) = indxi
      indxj=indxj+1
   END IF
END DO
ia = qcols - newnz

DEALLOCATE(ah,bh,ipiv,zerocols,c,STAT=status)
IF (status .NE. 0) THEN
   !WRITE(*,*) '**** IN BUILD_A(): COULD NOT DEALLOCATE MEMORY ****'
   !RETURN
END IF
END SUBROUTINE
