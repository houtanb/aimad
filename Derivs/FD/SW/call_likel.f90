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

IMPLICIT NONE
INTEGER :: n, numfiles, indx, mult
INTEGER, ALLOCATABLE, DIMENSION(:) :: iuser, iwork
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: converge, iterarr, objfuncallsarr
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: user, objgrd, work, xvals, xvalsorig, x
DOUBLE PRECISION, ALLOCATABLE, TARGET, DIMENSION(:,:) :: objgrdd, hdat, zdatadat, psiddat, psid0dat
DOUBLE PRECISION, ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: hddat, hd0dat
DOUBLE PRECISION, ALLOCATABLE, TARGET, DIMENSION(:,:,:,:) :: hdddat
DOUBLE PRECISION, ALLOCATABLE, TARGET, DIMENSION(:) :: psidat, uprbndxdat, lwrbndxdat
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: psiminvals1, psiminvals2, minlikel, psiminvals3, psiminvals4, psiminvals5
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: psiminvals6, minlikelgrd
DOUBLE PRECISION :: objf, temp

DOUBLE PRECISION :: condn, uprbnd
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: cofb
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: cofbd, cofbd0
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: cofbdd

! changes to length will have to be manually propagated to other programs down the chain
CHARACTER(len=20), ALLOCATABLE, TARGET, DIMENSION(:) :: endogdat, endogbakdat, errlistdat, observeddat
CHARACTER(len=100) :: inputfile
EXTERNAL objfun_likel, loaddata

!
! BEGIN MODIFY
!
! number of replications for Monte Carlo exercise
! need to have as many data files
numfiles = 1

! number of parameters being estimated
n = 6

! number of derivatives == n
nbdirs    = n
nbdirs0   = n

! number of equations in model
neq = 34

! number of innovations
nerrvars = 5

! number of observed series for estimation
nobservedvars = 4

! size of parameter vector
psin = 36

! number of observations
tt = 200

! max lag in model (fixed at 1 for now)
nlag = 1

! max lead in model (fixed at 1 for now)
nlead = 1

! Training sample for likelihood
ntrain = 10

! location of data file
inputfile = './zdatanew_112107.txt'


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
hrows = neq
hcols = neq * 3
nendogvars = neq
nendogvarsbak = nendogvars

ALLOCATE( hdat(hrows,hcols), hddat(nbdirs,hrows,hcols), hd0dat(nbdirs0,hrows,hcols), &
          hdddat(nbdirs0,nbdirs,hrows,hcols), zdatadat(nobservedvars,tt), x(n), &
          iuser(n), user(n), uprbndxdat(n), lwrbndxdat(n), objgrdd(nbdirs0, nbdirs), &
          objgrd(n), iwork(n+1), work(13*n), psid0dat(nbdirs0,psin),&
          psidat(psin), psiddat(nbdirs,psin), xvals(n), xvalsorig(n), &
          psiminvals1(1,numfiles), psiminvals2(1,numfiles), &
          minlikel(1,numfiles), converge(1,numfiles), objfuncallsarr(1,numfiles), &
          iterarr(1,numfiles), endogdat(nendogvars), endogbakdat(nendogvarsbak), &
          errlistdat(nerrvars), observeddat(nobservedvars), psiminvals3(1,numfiles), &
          psiminvals4(1,numfiles), psiminvals5(1,numfiles), psiminvals6(1,numfiles), &
          minlikelgrd(numfiles,n), STAT=status )
IF (status .NE. 0) THEN
   WRITE(*, *)
   WRITE(*, *) '***In Call_....: Could not allocate memeory***'
END IF

endog => endogdat
endogbak => endogbakdat
errlist => errlistdat
observed => observeddat
h => hdat
hd => hddat
hd0 => hd0dat
hdd => hdddat
zdata => zdatadat
psi => psidat
psid => psiddat
psid0 => psid0dat
uprbndx => uprbndxdat
lwrbndx => lwrbndxdat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! labels of observed variables
! observedlabels = char('ggdp','dpq','dw','rs');
observed(1) = 'ggdp'
observed(2) = 'dpq'
observed(3) = 'dw'
observed(4) = 'rs'


! labels of innovations
errlist(1) ='errgzpt'
errlist(2) ='errtaxl'
errlist(3) ='errtaxk'
errlist(4) ='errgov'
errlist(5) ='errmon'

! labels of endogenous variables
endog(1) ='lamc'
endog(2) ='cc'
endog(3) ='cc2gdp'
endog(4) ='invest'
endog(5) ='invest2gdp'
endog(6) ='tobq'
endog(7) ='kp'
endog(8) ='rentk'
endog(9) ='gdp'
endog(10) ='rs'
endog(11) ='mcq'
endog(12) ='mpl'
endog(13) ='wmarkup'
endog(14) ='mrs'
endog(15) ='dpq'
endog(16) ='dw'
endog(17) ='zetap'
endog(18) ='hw'
endog(19) ='util'
endog(20) ='gbt'
endog(21) ='ggdp'
endog(22) ='gsresid'
endog(23) ='rrs'
endog(24) ='gzt'
endog(25) ='gzpt'
endog(26) ='taxlshk'
endog(27) ='taxkshk'
endog(28) ='govshk'
endog(29) ='monshk'
endog(30) ='errgzpt'
endog(31) ='errtaxl'
endog(32) ='errtaxk'
endog(33) ='errgov'
endog(34) ='errmon'


endogbakdat = endogdat


! initial values of parameters
psi(:)= 0.0d0
psi(1)=0.000000000000000d0      !utilcon
psi(2)=0.010000000000000d0      !mu1
psi(3)=0.997000000000000d0      !beta
psi(4)=1.000000000000000d0      !gzbar
psi(5)=0.025000000000000d0      !delta
psi(6)=0.000000000000000d0      !tauk
psi(7)=0.280000000000000d0      !taul
psi(8)=0.180000000000000d0      !shrgy
psi(9)=0.000000000000000d0      !phik
psi(10)=3.000000000000000d0      !phii
psi(11)=0.500000000000000d0      !phic
psi(12)=10.000000000000000d0      !chi
psi(13)=-1.000000000000000d0      !cessub
psi(14)=1.000000000000000d0      !reltech
psi(15)=0.160000000000000d0      !shriy
psi(16)=0.000000000000000d0      !gami
psi(17)=1.500000000000000d0      !gampi
psi(18)=0.500000000000000d0      !gamy
psi(19)=0.750000000000000d0      !psip
psi(20)=0.750000000000000d0      !psiw
psi(21)=0.000000000000000d0      !lampb
psi(22)=1.000000000000000d0      !lampf
psi(23)=0.000000000000000d0      !lamwb
psi(24)=1.000000000000000d0      !lamwf
psi(25)=0.100000000000000d0      !thetap
psi(26)=0.100000000000000d0      !thetaw
psi(27)=0.000000000000000d0      !rho1
psi(28)=0.980000000000000d0      !rho4
psi(29)=0.970000000000000d0      !rho5
psi(30)=0.980000000000000d0      !rho6
psi(31)=0.950000000000000d0      !rho7
psi(32)=0.949500000000000d0      !sdev1
psi(33)=3.875900000000000d0      !sdev4
psi(34)=0.800000000000000d0      !sdev5
psi(35)=0.300000000000000d0      !sdev6
psi(36)=0.112500000000000d0      !sdev7


!
! iuser needs indexes of vector of parameters corresponding to the subset of parameters to be estimated.
! Case   :   alphap
! paramlabels = char('rho7','sdev7','gampi','gamy','psip','psiw');
!upperbound =      [ .999   10      10      5      .950   .900 ]';
!lowerbound =      [-.999   .01     1.05    .001   .2     .2   ]';

iuser   = (/ 31, 36, 17, 18, 19, 20/)
uprbndx = (/ .999d0,   10.0d0,      10.0d0,      5.0d0,      .950d0,   .900d0 /)
lwrbndx = (/ -.999d0,   .01d0,      1.05d0,     .001d0,      .2d0,     .2d0   /)

! change starting values
psi(31) =   0.957358093397858d0
psi(36) =   0.109832653161642d0
psi(17) =   1.365355836327819d0
psi(18) =   0.378513720551745d0
psi(19) =   0.720151969764918d0
psi(20) =   0.844155838978084d0


xvals(:) = psi(iuser)
xvalsorig = psi(iuser)

!
! END MODIFY
!




!
! CALL MINIMIZER
!
iterarr(:,:) = 0
objfuncallsarr(:,:) = 0
converge(:,:) = 0
minlikel(:,:) = 0.0d0
minlikelgrd(:,:) = 0.0d0
psiminvals1(:,:) = 0.0d0


WRITE(*,*) "ABOUT TO LOAD DATA"


CALL loaddata( zdatadat, nobservedvars, tt, inputfile )


!CALL GET_HMAT_DV_DV( hdat, hd0dat, hddat, hdddat, hrows, hcols, psi, psid0, psid, psin, nbdirs, nbdirs0 )
!CALL writemat_screen(h,neq,3*neq,1,neq,1,3*neq)

!WRITE (*,*) "zdatadat"
!CALL writemat_screen(zdatadat,tt,1,1,tt,1,1)
!WRITE (*,*) "endog"
!WRITE (*,*) endog(1:neq)
!WRITE (*,*) "observeddat"
!WRITE (*,*) observeddat(1:nobservedvars)

condn = 0.000000001d0
uprbnd = 1.0000001d0

ALLOCATE(cofb(neq,neq), cofbd(nbdirs,neq,neq), cofbd0(nbdirs0,neq,neq), cofbdd(nbdirs0,nbdirs,neq,neq), STAT = status)
IF (status .NE. 0) THEN
    WRITE(*,*) " COULD NOT ALLOCATE MEMORY FOR COFB"
END IF

!WRITE(*,*) "hdat"
!CALL writemat(hdat(:,:),neq,3*neq,1,neq,1,3*neq,"D:/fortran/sigma2_d/hmat.txt")

!CALL AIM_DV_DV(hdat, hd0dat, hddat, hdddat, neq, nlag, nlead, cofb, cofbd0, cofbd&
!&           , cofbdd, condn, uprbnd, info, nbdirs, nbdirs0)

!CALL likelihood_dv_dv( objf, objgrd, objgrdd, &
!     endog, neq, &
!     errlist, nerrvars, &
!     observeddat, nobservedvars, &
!     zdatadat, tt, &
!     hdat, hd0dat, hddat, hdddat, hrows, hcols, &
!     nlead, nlag, neq, ntrain, nbdirs, nbdirs0)



!CALL objfun_likel(2, n, xvals, objf, objgrd, objgrdd, 1, iuser, user )

CALL writemat_total(objgrdd,nbdirs0,nbdirs,'objgrdd.txt')

WRITE(*,*) 'likel'
WRITE(*,*) objf
WRITE(*,*) 'likelgrd'
WRITE(*,*) objgrd

!CALL writemat_total(objgrdd,nbdirs0,nbdirs,'objgrdd.txt')

!CALL objfun_likel_2pt_fd(n, xvals, objf, objgrd, iuser, user )




DEALLOCATE( h, hd, hd0, hdd, zdata, psi, psid, psid0, iuser, user, &
            objgrd, iwork, work, uprbndx, lwrbndx, objfuncallsarr, &
            xvals, xvalsorig, psiminvals1, psiminvals2, &
            converge, iterarr, endogdat, minlikelgrd, &
            endogbakdat, errlistdat, observeddat, x, &
            psiminvals3, psiminvals4, psiminvals5, psiminvals6, cofb, cofbd, cofbd0, cofbdd, STAT=status )

IF (status .NE. 0) THEN
   WRITE(*, *)
   WRITE(*, *) '***In Call_nag_likel: Could not deallocate memeory***'
END IF

END PROGRAM
