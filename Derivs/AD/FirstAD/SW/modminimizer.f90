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

MODULE modminimizer
!
! Purpose: To create storage for likelihood_d vars
!          that are needed in objfun
IMPLICIT NONE
SAVE
INTEGER :: neq, nlag, nlead, nobservedvars, nbdirs
INTEGER :: status, ntrain, nerrvars, psin, objfuncalls
INTEGER :: info, tt, nendogvars, nendogvarsbak, hrows, hcols
CHARACTER(len=20), POINTER, DIMENSION(:) :: endog, endogbak, errlist, observed
DOUBLE PRECISION, POINTER, DIMENSION(:,:) :: h, zdata, psid
DOUBLE PRECISION, POINTER, DIMENSION(:,:,:) :: hd
DOUBLE PRECISION, POINTER, DIMENSION(:) :: psi, uprbndx, lwrbndx

END MODULE modminimizer
