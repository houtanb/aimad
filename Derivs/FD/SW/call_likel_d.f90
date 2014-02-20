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

PROGRAM call_likel_d

USE modminimizer
USE moderror

IMPLICIT NONE
INTEGER :: n, numfiles
INTEGER, ALLOCATABLE, DIMENSION(:) :: iuser, iwork
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: converge, iterarr
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: user, objgrd, work, xvals, xvalsorig
DOUBLE PRECISION, ALLOCATABLE, TARGET, DIMENSION(:,:) :: hdat, zdatadat, hddat
DOUBLE PRECISION, ALLOCATABLE, TARGET, DIMENSION(:) :: psidat, psiddat, uprbndxdat, lwrbndxdat
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: psiminvals1, psiminvals2, psiminvals3, minlikel
DOUBLE PRECISION :: objf
CHARACTER(len=20), ALLOCATABLE, TARGET, DIMENSION(:) :: endogdat, endogbakdat, errlistdat, observeddat
CHARACTER(len=100) :: inputfile
EXTERNAL loaddata


numfiles = 1
! number of parameters being estimated
n = 4

! number of equations in model
neq = 48

! number of innovations
nerrvars = 7

! number of observed series for estimation
nobservedvars = 7

! size of parameter vector
psin = 36

! number of observations
tt = 230

! max lag in model (fixed at 1 for now)
nlag = 1

! max lead in model (fixed at 1 for now)
nlead = 1

! Training sample for likelihood
ntrain = 40

! location of data file
inputfile = './zdatanew_112107.txt'


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
hrows = neq
hcols = neq * 3
nendogvars = neq
nendogvarsbak = nendogvars

ALLOCATE( hdat(hrows,hcols), hddat(hrows,hcols), &
          zdatadat(nobservedvars,tt), &
          iuser(n), user(n), uprbndxdat(n), lwrbndxdat(n), &
          objgrd(n), iwork(n+1), work(13*n), &
          psidat(psin), psiddat(psin), xvals(n), xvalsorig(n), &
          psiminvals1(1,numfiles), psiminvals2(1,numfiles), &
          psiminvals3(1,numfiles), minlikel(1,numfiles), converge(1,numfiles), &
          iterarr(1,numfiles), endogdat(nendogvars), endogbakdat(nendogvarsbak), &
          errlistdat(nerrvars), observeddat(nobservedvars), STAT=status )
IF (status .NE. 0) THEN
   WRITE(*, *)
   WRITE(*, *) '***In Call_nag_....: Could not allocate memeory***'
END IF
endog => endogdat
endogbak => endogbakdat
errlist => errlistdat
observed => observeddat
h => hdat
zdata => zdatadat
psi => psidat
uprbndx => uprbndxdat
lwrbndx => lwrbndxdat

WRITE(*,*) 'IN CALL FILE'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! labels of observed variables
observed(1) = 'dy'
observed(2) = 'dc'
observed(3) = 'dinve'
observed(4) = 'labobs'
observed(5) = 'pinfobs'
observed(6) = 'dw'
observed(7) = 'robs'


! labels of innovations
errlist(1) ='ea'
errlist(2) ='eb'
errlist(3) ='eg'
errlist(4) ='eqs'
errlist(5) ='em'
errlist(6) ='epinf'
errlist(7) ='ew'


! labels of endogenous variables
endog(1) ='rrf'
endog(2) ='zcapf'
endog(3) ='rkf'
endog(4) ='kf'
endog(5) ='invef'
endog(6) ='pkf'
endog(7) ='cf'
endog(8) ='yf'
endog(9) ='labf'
endog(10) ='wf'
endog(11) ='kpf'
endog(12) ='mc'
endog(13) ='zcap'
endog(14) ='rk'
endog(15) ='k'
endog(16) ='inve'
endog(17) ='pk'
endog(18) ='c'
endog(19) ='y'
endog(20) ='lab'
endog(21) ='pinf'
endog(22) ='w'
endog(23) ='r'
endog(24) ='a'
endog(25) ='b'
endog(26) ='g'
endog(27) ='qs'
endog(28) ='ms'
endog(29) ='spinf'
endog(30) ='epinfma'
endog(31) ='sw'
endog(32) ='ewma'
endog(33) ='ea'
endog(34) ='eb'
endog(35) ='eg'
endog(36) ='eqs'
endog(37) ='em'
endog(38) ='epinf'
endog(39) ='ew'
endog(40) ='kp'
endog(41) ='trend'
endog(42) ='dy'
endog(43) ='dc'
endog(44) ='dinve'
endog(45) ='dw'
endog(46) ='pinfobs'
endog(47) ='robs'
endog(48) ='labobs'


endogbakdat = endogdat


! initial values of parameters
psi(:)= 0.0d0
psi(1)=0.478859915582470d0     !std_ea
psi(2)=0.125859564965700d0     !std_eb
psi(3)=0.523582271392230d0     !std_eg
psi(4)=0.574719103350110d0     !std_eqs
psi(5)=0.228878387020840d0     !std_em
psi(6)=0.134472416501310d0     !std_epinf
psi(7)=0.265486650994870d0     !std_ew
psi(8)=0.982655861241380d0     !crhoa
psi(9)=0.763222905918710d0     !crhob
psi(10)=0.989804522122800d0     !crhog
psi(11)=0.885954430842400d0     !crhoqs
psi(12)=0.010520826221328d0     !crhoms
psi(13)=0.991323101658320d0     !crhopinf
psi(14)=0.854954431246090d0     !crhow
psi(15)=0.920072190583200d0     !cmap
psi(16)=0.833520311036920d0     !cmaw
psi(17)=2.000137164117200d0     !csadjcost
psi(18)=1.790515772889400d0     !csigma
psi(19)=0.440751573697940d0     !chabb
psi(20)=0.919898921933850d0     !cprobw
psi(21)=1.868683475830600d0     !csigl
psi(22)=0.736991490615650d0     !cprobp
psi(23)=0.905399006351010d0     !cindw
psi(24)=0.014956549293304d0     !cindp
psi(25)=0.686118976726430d0     !czcap
psi(26)=1.592490054125300d0     !cfc
psi(27)=2.999983870806800d0     !crpi
psi(28)=0.901174891550580d0     !crr
psi(29)=0.192563764019040d0     !cry
psi(30)=0.258690445333910d0     !crdy
psi(31)=1.496470115650900d0     !constepinf
psi(32)=0.001099584597693d0     !constebeta
psi(33)=-3.372042064854000d0    !constelab
psi(34)=0.452186613659850d0     !ctrend
psi(35)=0.582078105008180d0     !cgy
psi(36)=0.190769838500580d0     !calfa


!
! iuser needs indexes of vector of parameters corresponding to the subset of parameters to be estimated.
! Case   :   alphap
! paramlabels = char('rho7','sdev7','gampi','gamy','psip','psiw');
!upperbound =      [ .999   10      10      5      .950   .900 ]';
!lowerbound =      [-.999   .01     1.05    .001   .2     .2   ]';

iuser   = (/ 19,     20,    26,    36/)
uprbndx = (/ 0.99d0, 0.92d0,3.0d0, 1.00d0/)
lwrbndx = (/ 0.001d0,0.30d0,1.0d0, 0.01d0/)

xvals(:) = psi(iuser)
xvalsorig = psi(iuser)


!
! END MODIFY
!






WRITE(*,*) "ABOUT TO LOAD DATA"


!CALL loaddata( zdatadat, nobservedvars, tt, inputfile )



objf = 0.0d0
objgrd(:) = 0.0d0
CALL loaddata( zdata, nobservedvars, tt, inputfile )
WRITE(*,*) "loaded data"

CALL objfun_likel( n, xvals, objf, iuser, user )
WRITE(*,*) 'objf'
WRITE(*,*) objf

CALL objfun_likel_2pt_fd( n, xvals, objf, objgrd, iuser, user )
   WRITE(*,*) 'objf'
   WRITE(*,*) objf
   WRITE(*,*) 'objgrd'
   WRITE(*,*) objgrd

DEALLOCATE( hdat, hddat, &
          zdatadat, &
          iuser, user, uprbndxdat, lwrbndxdat, &
          objgrd, iwork, work, &
          psidat, psiddat, xvals, xvalsorig, &
          psiminvals1, psiminvals2, &
          psiminvals3, minlikel, converge, &
          iterarr, endogdat, endogbakdat, &
          errlistdat, observeddat, STAT=status )
IF (status .NE. 0) THEN
   WRITE(*, *)
   WRITE(*, *) '***In Call_nag_likel: Could not deallocate memeory***'
END IF
END PROGRAM
