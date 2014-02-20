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

SUBROUTINE numeric_shiftsvd(h,q,iq,qrows,qcols,neq,condn,nnumeric)

IMPLICIT NONE

INTEGER, INTENT(IN) :: qrows, qcols, neq
INTEGER, INTENT(INOUT) :: iq
INTEGER, INTENT(OUT) :: nnumeric
DOUBLE PRECISION, DIMENSION(neq,neq*3), INTENT(INOUT) :: h
DOUBLE PRECISION, DIMENSION(qrows,qcols), INTENT(INOUT) :: q
DOUBLE PRECISION, INTENT(IN) :: condn

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: hpart, hpartcopy, v, hpart2, hprod, u
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: d, work
LOGICAL, ALLOCATABLE, DIMENSION(:) :: stablerows
INTEGER, ALLOCATABLE, DIMENSION(:) :: selector, order
INTEGER status, lleft, rleft, lright, rright, ns, lwork, info, i
DOUBLE PRECISION :: zerotol
EXTERNAL shiftright, dgesvd, constructq, get_stable_singvals, get_selector, dgemm

info = 0
lleft = 1
nnumeric = 0
rleft = qcols
lright = qcols+1
rright = qcols+neq
zerotol = 0.00000000010d0
lwork = MAX( 1, 3*neq+neq, 5*neq )

ALLOCATE( order(neq), hpart(neq,neq), hpartcopy(neq,neq), &
     v(neq,neq), u(neq,neq), d(neq), work(lwork), stablerows(neq), &
     hprod(neq,3*neq), STAT=status )
IF (status .NE. 0) THEN
   WRITE(*,*) '**** IN NUMERIC_SHIFT():: COULD NOT ALLOCATE MEMORY ****'
   STOP
END IF

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
IF( info .NE. 0 ) THEN
   WRITE(*,*)
   WRITE(*,*) '**** IN NUMERIC_SHIFT():: Call to dgesvd() returned info != 0'
END IF
hpart=hpartcopy
u = u( :, order )
d = d( order )
! v=TRANSPOSE(v) NOTE: only necessary for derivative::


! if the ith entry on the diagonal of hpart is <= condn, the ith entry of
! stablerows is true
CALL get_stable_singvals( d, neq, condn, stablerows, ns )
DO
   IF ( ANY(stablerows) .AND. (iq .LE. qrows) ) THEN

       ALLOCATE( selector(ns), STAT=status )
       IF( status .NE. 0 ) THEN
          WRITE(*,*)
          WRITE(*,*) '**** IN NUMERIC_SHIFT():: COULD NOT ALLOCATE MEMORY ****'
          STOP
       END IF
       CALL get_selector(stablerows,neq,selector,ns)

       ! h = u'*h
       CALL dgemm( 'T', 'N', neq, 3*neq, neq, 1.d0, u, neq, h, neq, 0.d0, hprod, neq )
       h = hprod


       IF (iq+ns .LE. qrows) THEN
          q( iq+1:iq+ns, : ) = h( selector, lleft:rleft )
       ELSE
          WRITE(*,*)
          WRITE(*,*) '*****  problem in numeric_shift *****'
          WRITE(*,*)
          STOP
       END IF

       ALLOCATE( hpart2(ns,3*neq), STAT=status )
       IF (status .NE. 0) THEN
          WRITE(*,*)
          WRITE(*,*) '**** IN NUMERIC_SHIFT():: COULD NOT ALLOCATE MEMORY ****'
          STOP
       END IF

       CALL shiftright( h(selector,:), ns, 3*neq, neq, hpart2 )
       h(selector,:) = hpart2

       DEALLOCATE( hpart2, selector, STAT=status )
       IF (status .NE. 0) THEN
          WRITE(*,*)
          WRITE(*,*) '**** IN NUMERIC_SHIFT():: COULD NOT DEALLOCATE MEMORY ****'
          STOP
       END IF


       iq = iq + ns
       nnumeric = nnumeric + ns

       hpart = h(:,lright:rright)
       hpartcopy = hpart
       info = 0

       !
       ! CALL SVD
       !
       CALL dgesvd( 'A', 'A', neq, neq, hpart, neq, d, u, neq, v, neq, work, lwork, info )
       IF( info .NE. 0 ) THEN
          WRITE(*,*)
          WRITE(*,*) '**** IN NUMERIC_SHIFT():: Call in loop to dgesvd() returned info != 0'
       END IF
       hpart = hpartcopy
       u = u( :, order )
       d = d( order )
       ! v=TRANSPOSE(v) NOTE: probably not necessary to calculate
       ! dd or dv since they don't feed into the end product

       CALL get_stable_singvals( d, neq, condn, stablerows, ns )
   ELSE
      EXIT
   END IF
END DO


DEALLOCATE( order, hpart, hpartcopy, u, d, v, work, stablerows, hprod, STAT=status)
IF (status .NE. 0) THEN
   WRITE(*,*) '**** IN NUMERIC_SHIFT():: COULD NOT DEALLOCATE MEMORY ****'
   STOP
END IF

END SUBROUTINE
