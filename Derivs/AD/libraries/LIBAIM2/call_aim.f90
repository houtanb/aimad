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

PROGRAM call_aim

IMPLICIT NONE

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: h, cofb

INTEGER :: neq, nlag, nlead, indxi
INTEGER :: status
INTEGER :: info
DOUBLE PRECISION :: condn, uprbnd
REAL :: time0, time1
!CHARACTER(len=100) :: filename = 'houtan.txt'
CHARACTER(len=100) :: filename = 'sigma2_estim.txt'


EXTERNAL loadcof, aim, printmat

!rbc_website
!neq = 4

! houtan
!neq = 24

!sigma2_estim
neq = 142

nlag = 1
nlead = 1

condn = 0.000000001
uprbnd = 1+0.0000001



ALLOCATE(h(neq,neq*3),cofb(neq,neq),STAT=status )
!,q(qrows,qcols),STAT=status)

CALL loadcof(h,neq,filename)

CALL CPU_TIME(time0)

DO indxi = 1,10
WRITE(*,*) indxi
CALL aim(h,neq,nlag,nlead,cofb,condn,uprbnd,info)
END DO
CALL CPU_TIME(time1)
WRITE(*,*) 'TIME0',time0
WRITE(*,*) 'TIME1',time1
WRITE(*,*) 'TIME ELAPSED',time1-time0

CALL printmat(cofb,neq,neq)


DEALLOCATE(h,cofb,STAT=status)


END PROGRAM
