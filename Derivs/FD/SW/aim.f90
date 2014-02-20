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

SUBROUTINE aim(horig,neq,nlag,nlead,cofb,condn,uprbnd,info)

IMPLICIT NONE

INTEGER, INTENT(IN) :: neq, nlag, nlead
DOUBLE PRECISION, INTENT(IN) :: condn, uprbnd
DOUBLE PRECISION, INTENT(IN), DIMENSION(neq,3*neq) :: horig
DOUBLE PRECISION, INTENT(OUT), DIMENSION(neq, neq) :: cofb

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: h, q, a
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: asmall, w
INTEGER, ALLOCATABLE, DIMENSION(:) :: js

INTEGER :: qrows, qcols
INTEGER :: ia

INTEGER :: iq,nexact,nnumeric

INTEGER :: status
INTEGER :: nlarge

INTEGER :: info


EXTERNAL exact_shift, numeric_shift, build_a, eigensystem, reduced_form, copy_w

qrows = neq
qcols = neq*2

iq = 0

ALLOCATE(h(neq,neq*3),q(qrows,qcols),STAT=status)

h=horig
q(:,:) = 0.0d0

CALL exact_shift(h,q,iq,qrows,qcols,neq,nexact)

CALL numeric_shift(h,q,iq,qrows,qcols,neq,condn,nnumeric)


ALLOCATE(a(qcols,qcols),js(qcols),STAT=status)
CALL build_a(h,qcols,neq,a,ia,js)

ALLOCATE(asmall(ia,ia),w(ia,ia),STAT=status)

IF (ia > 0) THEN
   asmall = a(js(1:ia),js(1:ia))
   CALL eigensystem(asmall,ia,uprbnd,w,nlarge,info)
ELSE
   WRITE(*,*) '**** Problem occurred ****'
END IF

CALL copy_w(q,qrows,qcols,w,ia,iq,js(1:ia))

IF( (nnumeric+nlarge+nexact) .NE. qrows ) THEN
!   WRITE(*,*) 'IN AIM(): (nnumeric+nlarge+nexact) .NE. qrows )'
!   WRITE(*,*) 'nnumeric+nlarge+nexact = ',(nnumeric+nlarge+nexact)
!   WRITE(*,*) 'qrows = ',qrows
   RETURN
END IF

CALL reduced_form(q,qrows,qcols,neq,cofb,info)
call writemat(cofb,neq,neq,1,neq,1,neq,'cofb.txt')

IF (info .NE. 0) THEN
   WRITE(*,*) ' **** problems in reduced_form ****'
END IF

DEALLOCATE(h,q,a,js,asmall,w,STAT=status)

END SUBROUTINE
