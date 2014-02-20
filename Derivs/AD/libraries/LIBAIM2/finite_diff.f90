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

PROGRAM finite_diff

IMPLICIT NONE

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: h, cofb, cofbp, cofbm, cofbd

INTEGER :: neq, nlag, nlead, rowderiv, colderiv
INTEGER :: status
INTEGER :: info

DOUBLE PRECISION :: condn, uprbnd, stepsize
REAL :: time0, time1
!CHARACTER(len=100) :: filename = 'rbc_website.txt'
!CHARACTER(len=100) :: filename = 'houtan.txt'
CHARACTER(len=100) :: filename = 'sigma2_estim.txt'

EXTERNAL loadcof, aim, gettime



! for rbc_website
!neq = 4
!rowderiv = 3
!colderiv = 2

! for houtan
!neq = 24
!rowderiv = 4
!colderiv = 5

! for sigma2_estim
neq = 142
rowderiv = 2
colderiv = 2


nlag = 1
nlead = 1

condn = 0.000000001
uprbnd = 1+0.0000001

stepsize = 0.0000000001

ALLOCATE(h(neq,neq*3),cofb(neq,neq),STAT=status )
!,q(qrows,qcols),STAT=status)

CALL loadcof(h,neq,filename)

CALL CPU_TIME(time0)
!CALL gettime(time0)

ALLOCATE(cofbp(neq,neq),cofbm(neq,neq),cofbd(neq,neq),STAT=status)



CALL aim(h,neq,nlag,nlead,cofb,condn,uprbnd,info)
h(rowderiv,colderiv) = h(rowderiv,colderiv) + stepsize
CALL aim(h,neq,nlag,nlead,cofbp,condn,uprbnd,info)
h(rowderiv,colderiv) = h(rowderiv,colderiv) - 2*stepsize
CALL aim(h,neq,nlag,nlead,cofbm,condn,uprbnd,info)


cofbd = (cofbp-cofbm)/(2*stepsize)


CALL CPU_TIME(time1)
!CALL gettime(time1)

WRITE(*,*) 'TIME0',time0
WRITE(*,*) 'TIME1',time1
WRITE(*,*) 'TIME ELAPSED',time1-time0


DEALLOCATE(h,cofb,STAT=status)


END PROGRAM
