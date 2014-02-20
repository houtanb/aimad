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

SUBROUTINE numeric_shift(h,q,iq,qrows,qcols,neq,condn,nnumeric)

!IMPLICIT NONE

!INTEGER, INTENT(IN) :: qrows, qcols, neq
!INTEGER, INTENT(INOUT) :: iq
!INTEGER, INTENT(OUT) :: nnumeric
INTEGER :: qrows, qcols, neq, iq, nnumeric

!DOUBLE PRECISION, DIMENSION(neq,neq*3), INTENT(INOUT) :: h
!DOUBLE PRECISION, DIMENSION(qrows,qcols), INTENT(INOUT) :: q
DOUBLE PRECISION, DIMENSION(neq,neq*3) :: h
DOUBLE PRECISION, DIMENSION(qrows,qcols):: q
!DOUBLE PRECISION, INTENT(IN) :: condn
DOUBLE PRECISION :: condn

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: hpart, hpartcopy, v, hpart2, hprod, u, htmp
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: d, work
LOGICAL, ALLOCATABLE, DIMENSION(:) :: stablerows
INTEGER, ALLOCATABLE, DIMENSION(:) :: selector, order
INTEGER status, lleft, rleft, lright, rright, ns, lwork, info, i, indx, j

EXTERNAL shiftright, dgesvd, get_stable_singvals, get_selector, dgemm
INTRINSIC ANY, MAX

info = 0
lleft = 1
nnumeric = 0
rleft = qcols
lright = qcols+1
rright = qcols+neq
lwork = MAX(1,3*neq+neq,5*neq)
status = 0

ALLOCATE( order(neq), hpart(neq,neq), hpartcopy(neq,neq), v(neq,neq), u(neq,neq), &
          d(neq), work(lwork), stablerows(neq), hprod(neq,3*neq), STAT=status )


hpart = h(:,lright:rright)
hpartcopy = hpart


! order is an array with:
! order(1) = neq;
! order(2) = neq-1;
! ...
! order(neq) = 1
order = (/ ( i, i = neq, 1, -1 ) /)

!
! CALL SVD
!
CALL dgesvd( 'A', 'A', neq, neq, hpart, neq, d, u, neq, v, neq, work, lwork, info )
hpart=hpartcopy
u = u( :, order )
d = d( order )

! if the ith entry on the diagonal of hpart is <= condn, the ith entry of
! stablerows is true
CALL get_stable_singvals( d, neq, condn, stablerows, ns )

DO WHILE( ANY(stablerows) .AND. (iq .LE. qrows) )

       ALLOCATE( selector(ns), STAT=status )
       selector(:) = 0.0d0
       CALL get_selector( stablerows, neq, selector, ns )

       ! h = u'*h
       CALL dgemm( 'T', 'N', neq, 3*neq, neq, 1.d0, u, neq, h, neq, 0.d0, hprod, neq )
       h = hprod

       IF (iq+ns .LE. qrows) THEN
          DO j=lleft,rleft
             DO i=1,ns
                indx = selector(i)
                q( iq + i, j-lleft+1 ) = h( indx, j )
             END DO
          END DO
       ELSE
          RETURN
       END IF


       ALLOCATE( hpart2(ns,3*neq), htmp(ns,3*neq), STAT=status )
       hpart2(:,:) = 0.0d0

       DO j=1,3*neq
          DO i=1,ns
             indx = selector(i)
             htmp( i, j ) = h( indx, j )
          END DO
       END DO

       CALL shiftright( htmp, ns, 3*neq, neq, hpart2 )

       DO j=1,3*neq
          DO i=1,ns
             indx = selector(i)
             h( indx, j ) = hpart2( i, j )
          END DO
       END DO

       DEALLOCATE( hpart2, selector, htmp, STAT=status )

       iq = iq + ns
       nnumeric = nnumeric + ns

       hpart = h(:,lright:rright)
       hpartcopy = hpart
       info = 0

       !
       ! CALL SVD
       !
       CALL dgesvd( 'A', 'A', neq, neq, hpart, neq, d, u, neq, v, neq, work, lwork, info )
       hpart = hpartcopy
       u = u( :, order )
       d = d( order )

       CALL get_stable_singvals( d, neq, condn, stablerows, ns )
END DO

DEALLOCATE( order, hpart, hpartcopy, u, d, v, work, stablerows, hprod, STAT=status)


END SUBROUTINE
