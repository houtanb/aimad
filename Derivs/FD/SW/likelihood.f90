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

!/********************************************************************
! *
! * SUBROUTINE likelihood():
! *
! * This routine evaluates the likelihood function for a given
! * economic model (as defined by the h matrix below) at a point
! * defined by the parallel arrays psiV, psiL, and the given data on
! * a set of variables determined by the model (contained in the
! * matrix zdata).
! *
! * The model matrix h is written as a function of a parameter
! * vector, psiV.
! *
! *
! ********************************************************************
! * CALLING STRUCTURE:
! *
! * int likelihood( double *likel,
! *		   char **endog, int nEndogVars,
! *		   char **errlist, int nErrVars,
! *	           char **observed, int nObservedVars,
! *		   double *zdata, int TT,
! *		   double *h, int hrows, int hcols,
! *		   int leads, int lags, int neq, int ntrain )
! *
! ********************************************************************/
SUBROUTINE likelihood( likel,                     &
                       endog, nendogvars,         &
                       errlist, nerrvars,         &
                       observed, nobservedvars,   &
                       zdata, tt,                 &
                       h, hrows, hcols,           &
                       leads, lags, neq, ntrain )

IMPLICIT NONE
INTEGER, INTENT(IN) :: tt, hrows, hcols, leads, lags, neq, ntrain
INTEGER, INTENT(IN) :: nendogvars, nerrvars, nobservedvars
DOUBLE PRECISION, INTENT(OUT) :: likel
DOUBLE PRECISION, INTENT(IN), DIMENSION(nobservedvars,tt) :: zdata
DOUBLE PRECISION, INTENT(IN), DIMENSION(hrows,hcols) :: h
CHARACTER(len=*), INTENT(IN), DIMENSION(nendogvars) :: endog
CHARACTER(len=*), INTENT(IN), DIMENSION(nerrvars) :: errlist
CHARACTER(len=*), INTENT(IN), DIMENSION(nobservedvars) :: observed

INTEGER :: status, info, err, diffinvars, zcount, i, nendogvarscopy
DOUBLE PRECISION :: condn, uprbnd
INTEGER, ALLOCATABLE, DIMENSION(:) :: selector
INTEGER :: nselected
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: cofb, Z, Q, Qcompare
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: amat, bmat, xi10, amatsmall, bmatsmall
CHARACTER(len=20), ALLOCATABLE, DIMENSION(:) :: endogcopy
INTEGER, ALLOCATABLE, DIMENSION(:) :: observedpos
INTEGER, EXTERNAL :: strmatch

!DOUBLE PRECISION :: error, checkaim
EXTERNAL aim, splitmat, kalman, dgemm

condn = 0.000000001d0
uprbnd = 1.0000001d0

nendogvarscopy = nendogvars

diffinvars = nendogvars - nerrvars
! Allocate memory
!
ALLOCATE( cofb(neq,neq), amat(diffinvars,diffinvars), &
          bmat(diffinvars, nerrvars), endogcopy(nendogvarscopy), &
          selector(nendogvarscopy), observedpos(nobservedvars), STAT=status )
IF( status .NE. 0 ) THEN
   WRITE(*,*)
   WRITE(*,*) "***In Likelihood: Could not allocate memeory***"
   info = -1
   err  = -1
END IF
endogcopy = endog

!!
!! CALL AIM
!!

cofb(:,:) = 0.0d0
info = 0


!CALL printmat( h, neq, neq*3 )

CALL aim( h, neq, lags, leads, cofb, condn, uprbnd, info )
!error = checkaim( neq, lags, leads, h, cofb )


!CALL writemat(cofb,neq,neq,1,neq,1,neq,"cofb.txt")
!!
!!CALL SPLITMAT
!!
CALL splitmat( cofb, endogcopy, nendogvarscopy, errlist, nerrvars, amat, bmat )
DEALLOCATE(cofb, STAT=status)

IF( status .NE. 0 ) THEN
   WRITE(*,*)
   WRITE(*,*) "***In Likelihood: Could not deallocate memory for cofb***"
   info = -1
   err  = -1
END IF
!CALL writemat(amat,neq-nerrvars,neq-nerrvars,1,neq-nerrvars,1,neq-nerrvars,"amat.txt")
!CALL writemat(bmat,neq-nerrvars,nerrvars,1,neq-nerrvars,1,nerrvars,"bmat.txt")


DO i=1,nobservedvars
   observedpos(i) =  strmatch(observed(i), endogcopy, nendogvarscopy)
END DO


CALL shrinkspacesmall(amat,bmat,nendogvarscopy,nerrvars,observedpos,nobservedvars,selector,nselected)


ALLOCATE(amatsmall(nselected,nselected),bmatsmall(nselected,nerrvars),STAT=status)
IF( status .NE. 0 ) THEN
   WRITE(*,*)
   WRITE(*,*) "***In Likelihood: Could not allocate memeory for amatsmall,bmatsmall***"
   info = -1
   err  = -1
END IF


amatsmall = amat(selector(1:nselected),selector(1:nselected))
bmatsmall = bmat(selector(1:nselected),:)
nendogvarscopy = nselected

!CALL writemat(amatsmall,nselected,nselected,1,nselected,1,nselected,"amatsmall.txt")
!CALL writemat(bmatsmall,nselected,nerrvars,1,nselected,1,nerrvars,"bmatsmall.txt")



DEALLOCATE(amat,bmat, observedpos, STAT = status)
IF( status .NE. 0 ) THEN
   WRITE(*,*)
   WRITE(*,*) "***In Likelihood: Could not deallocate memory for amat,bmat***"
   info = -1
   err  = -1
END IF


ALLOCATE( Q(nendogvarscopy, nendogvarscopy), Qcompare(nendogvarscopy, nendogvarscopy),STAT=status )
IF( status .NE. 0 ) THEN
   WRITE(*,*)
   WRITE(*,*) "***In Likelihood: Could not allocate memeory for Q***"
   info = -1
   err  = -1
END IF

Qcompare = MATMUL(bmatsmall,TRANSPOSE(bmatsmall))
CALL dgemm( 'N', 'T', nendogvarscopy, nendogvarscopy, nerrvars, 1.0d0, bmatsmall, nendogvarscopy, &
            bmatsmall, nendogvarscopy, 0.0d0, Q, nendogvarscopy )

DEALLOCATE(bmatsmall, Qcompare,STAT=status)
IF( status .NE. 0 ) THEN
   WRITE(*,*)
   WRITE(*,*) "***In Likelihood: Could not deallocate memeory for bmatsmall***"
   info = -1
   err  = -1
END IF



!!
!! zdataPrime
!!
ALLOCATE( Z(nendogvarscopy, nobservedvars), xi10(nendogvarscopy,1), STAT=status )
IF( status .NE. 0 ) THEN
   WRITE(*,*)
   WRITE(*,*) "***In Likelihood: Could not allocate memeory for Z, xi10***"
   info = -1
   err  = -1
END IF
xi10(:,:) = 0.0d0
Z(:,:) = 0.0d0
zcount = 0
DO i=1,nobservedvars
   Z( strmatch(observed(i), endogcopy(selector(1:nselected)), nendogvarscopy), i ) = 1.0d0
END DO

! special line for Smets Wouters model
xi10(strmatch('trend',endogcopy(selector(1:nselected)), nendogvarscopy),1) = 1.0d0


!!
!! Kalman
!!
CALL kalman( amatsmall, nendogvarscopy, Z, nobservedvars, zdata, tt, xi10, Q, ntrain, likel, 1)

DEALLOCATE(amatsmall, Z, Q, xi10, endogcopy, selector, STAT=status)
IF( status .NE. 0 ) THEN
   WRITE(*,*)
   WRITE(*,*) "***In Likelihood: Could not allocate memory***"
   info = -1
   err  = -1
END IF
END SUBROUTINE
