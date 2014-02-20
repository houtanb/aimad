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

SUBROUTINE get_hmat( hmat, hrows, hcols, psi_array, psin )

IMPLICIT NONE

!
! Input Variables
!
INTEGER, INTENT(IN) :: hrows, hcols, psin
DOUBLE PRECISION, INTENT(IN), DIMENSION(psin) :: psi_array
DOUBLE PRECISION, INTENT(OUT), DIMENSION(hrows,hcols) :: hmat

!
! Local Variables
!
DOUBLE PRECISION :: std_ea, std_eb
DOUBLE PRECISION :: std_eg, std_eqs
DOUBLE PRECISION :: std_em, std_epinf
DOUBLE PRECISION :: std_ew, crhoa
DOUBLE PRECISION :: crhob, crhog
DOUBLE PRECISION :: crhoqs, crhoms
DOUBLE PRECISION :: crhopinf, crhow
DOUBLE PRECISION :: cmap, cmaw
DOUBLE PRECISION :: csadjcost, csigma
DOUBLE PRECISION :: chabb, cprobw
DOUBLE PRECISION :: csigl, cprobp
DOUBLE PRECISION :: cindw, cindp
DOUBLE PRECISION :: czcap, cfc
DOUBLE PRECISION :: crpi, crr
DOUBLE PRECISION :: cry, crdy
DOUBLE PRECISION :: constepinf, constebeta
DOUBLE PRECISION :: constelab, ctrend
DOUBLE PRECISION :: cgy, calfa

DOUBLE PRECISION :: ctou, clandaw
DOUBLE PRECISION :: cg, curvp
DOUBLE PRECISION :: curvw, cpie
DOUBLE PRECISION :: cgamma, cbeta
DOUBLE PRECISION :: clandap, cbetabar
DOUBLE PRECISION :: cr, crk
DOUBLE PRECISION :: cw, cikbar
DOUBLE PRECISION :: cik, clk
DOUBLE PRECISION :: cky, ciy
DOUBLE PRECISION :: ccy, crkky
DOUBLE PRECISION :: cwhlc, cwly
DOUBLE PRECISION :: conster

!
! Assign psi_array
!
std_ea     = psi_array(1)
std_eb     = psi_array(2)
std_eg     = psi_array(3)
std_eqs    = psi_array(4)
std_em     = psi_array(5)
std_epinf  = psi_array(6)
std_ew     = psi_array(7)
crhoa      = psi_array(8)
crhob      = psi_array(9)
crhog      = psi_array(10)
crhoqs     = psi_array(11)
crhoms     = psi_array(12)
crhopinf   = psi_array(13)
crhow      = psi_array(14)
cmap       = psi_array(15)
cmaw       = psi_array(16)
csadjcost  = psi_array(17)
csigma     = psi_array(18)
chabb      = psi_array(19)
cprobw     = psi_array(20)
csigl      = psi_array(21)
cprobp     = psi_array(22)
cindw      = psi_array(23)
cindp      = psi_array(24)
czcap      = psi_array(25)
cfc        = psi_array(26)
crpi       = psi_array(27)
crr        = psi_array(28)
cry        = psi_array(29)
crdy       = psi_array(30)
constepinf = psi_array(31)
constebeta = psi_array(32)
constelab  = psi_array(33)
ctrend     = psi_array(34)
cgy        = psi_array(35)
calfa      = psi_array(36)


!
! Assign Composite Parameters
!
ctou     = 0.0250000000000000d0
clandaw  = 1.5000000000000000d0
cg       = 0.1800000000000000d0
curvp    = 10.0000000000000000d0
curvw    = 10.0000000000000000d0
cpie     =       1+constepinf/100
cgamma   =       1+ctrend/100
cbeta    =       1/(1+constebeta/100)
clandap  =       cfc
cbetabar =       1/(1+constebeta/100)*(1+ctrend/100)**(-csigma)
cr       =       (1+constepinf/100)*(1+constebeta/100)/(1+ctrend/100)**(- &
csigma)
crk      =       (1+constebeta/100)*(1+ctrend/100)**csigma-39.D0/40.D0
cw       =       (calfa**calfa*(1-calfa)**(1-calfa)/cfc/((1+constebeta/ &
100)*(1+ctrend/100)**csigma-39.D0/40.D0)**calfa)**(1/(1-calfa))
cikbar   =       1-39.D0/40.D0/(1+ctrend/100)
cik      =       (1-39.D0/40.D0/(1+ctrend/100))*(1+ctrend/100)
clk      =       (1-calfa)/calfa*((1+constebeta/100)*(1+ctrend/ &
100)**csigma-39.D0/40.D0)/(calfa**calfa*(1-calfa)**(1-calfa)/cfc/((1+ &
constebeta/100)*(1+ctrend/100)**csigma-39.D0/40.D0)**calfa)**(1/(1-calfa))
cky      =       cfc*((1-calfa)/calfa*((1+constebeta/100)*(1+ctrend/ &
100)**csigma-39.D0/40.D0)/(calfa**calfa*(1-calfa)**(1-calfa)/cfc/((1+ &
constebeta/100)*(1+ctrend/100)**csigma-39.D0/40.D0)**calfa)**(1/(1- &
calfa)))**(-1+calfa)
ciy      =       (1-39.D0/40.D0/(1+ctrend/100))*(1+ctrend/100)*cfc*((1- &
calfa)/calfa*((1+constebeta/100)*(1+ctrend/100)**csigma-39.D0/40.D0)/ &
(calfa**calfa*(1-calfa)**(1-calfa)/cfc/((1+constebeta/100)*(1+ctrend/ &
100)**csigma-39.D0/40.D0)**calfa)**(1/(1-calfa)))**(-1+calfa)
ccy      =       41.D0/50.D0-(1-39.D0/40.D0/(1+ctrend/100))*(1+ctrend/ &
100)*cfc*((1-calfa)/calfa*((1+constebeta/100)*(1+ctrend/100)**csigma-39.D0/ &
40.D0)/(calfa**calfa*(1-calfa)**(1-calfa)/cfc/((1+constebeta/100)*(1+ &
ctrend/100)**csigma-39.D0/40.D0)**calfa)**(1/(1-calfa)))**(-1+calfa)
crkky    =       ((1+constebeta/100)*(1+ctrend/100)**csigma-39.D0/ &
40.D0)*cfc*((1-calfa)/calfa*((1+constebeta/100)*(1+ctrend/100)**csigma- &
39.D0/40.D0)/(calfa**calfa*(1-calfa)**(1-calfa)/cfc/((1+constebeta/100)*(1+ &
ctrend/100)**csigma-39.D0/40.D0)**calfa)**(1/(1-calfa)))**(-1+calfa)
cwhlc    =       (2.D0/3.D0-2.D0/3.D0*calfa)/calfa*((1+constebeta/100)*(1+ &
ctrend/100)**csigma-39.D0/40.D0)*cfc*((1-calfa)/calfa*((1+constebeta/ &
100)*(1+ctrend/100)**csigma-39.D0/40.D0)/(calfa**calfa*(1-calfa)**(1- &
calfa)/cfc/((1+constebeta/100)*(1+ctrend/100)**csigma-39.D0/ &
40.D0)**calfa)**(1/(1-calfa)))**(-1+calfa)/(41.D0/50.D0-(1-39.D0/40.D0/(1+ &
ctrend/100))*(1+ctrend/100)*cfc*((1-calfa)/calfa*((1+constebeta/100)*(1+ &
ctrend/100)**csigma-39.D0/40.D0)/(calfa**calfa*(1-calfa)**(1-calfa)/cfc/ &
((1+constebeta/100)*(1+ctrend/100)**csigma-39.D0/40.D0)**calfa)**(1/(1- &
calfa)))**(-1+calfa))
cwly     =       1-((1+constebeta/100)*(1+ctrend/100)**csigma-39.D0/ &
40.D0)*cfc*((1-calfa)/calfa*((1+constebeta/100)*(1+ctrend/100)**csigma- &
39.D0/40.D0)/(calfa**calfa*(1-calfa)**(1-calfa)/cfc/((1+constebeta/100)*(1+ &
ctrend/100)**csigma-39.D0/40.D0)**calfa)**(1/(1-calfa)))**(-1+calfa)
conster  =       100*(1+constepinf/100)*(1+constebeta/100)/(1+ctrend/ &
100)**(-csigma)-100
!
! Assign Hmat
!
hmat(:,:) = 0.0d0

hmat(5,5)=      -1/(1+cbetabar*cgamma)
hmat(7,7)=      -chabb/cgamma/(1+chabb/cgamma)
hmat(10,7)=      chabb/cgamma/(1-chabb/cgamma)
hmat(23,8)=      -crdy
hmat(4,11)=-1.0000000000000000d0
hmat(11,11)=      -1+cikbar
hmat(16,16)=      -1/(1+cbetabar*cgamma)
hmat(44,16)=1.0000000000000000d0
hmat(18,18)=      -chabb/cgamma/(1+chabb/cgamma)
hmat(22,18)=      1/(1+cbetabar*cgamma)/cprobw/((clandaw-1)*curvw+1)*chabb/ &
cgamma/(1-chabb/cgamma)-1/(1+cbetabar*cgamma)*cbetabar/((clandaw-1)*curvw+ &
1)*chabb/(1-chabb/cgamma)-1/(1+cbetabar*cgamma)/((clandaw-1)*curvw+ &
1)*chabb/cgamma/(1-chabb/cgamma)+cprobw*cbetabar/(1+cbetabar*cgamma)/ &
((clandaw-1)*curvw+1)*chabb/(1-chabb/cgamma)
hmat(43,18)=1.0000000000000000d0
hmat(23,19)=      crdy
hmat(42,19)=1.0000000000000000d0
hmat(21,21)=      -1/(1+cbetabar*cgamma*cindp)*cindp
hmat(22,21)=      -cindw/(1+cbetabar*cgamma)
hmat(22,22)=      -1/(1+cbetabar*cgamma)
hmat(45,22)=1.0000000000000000d0
hmat(23,23)=      -crr
hmat(24,24)=      -crhoa
hmat(25,25)=      -crhob
hmat(26,26)=      -crhog
hmat(27,27)=      -crhoqs
hmat(28,28)=      -crhoms
hmat(29,29)=      -crhopinf
hmat(29,30)=      cmap
hmat(31,31)=      -crhow
hmat(31,32)=      cmaw
hmat(24,33)=      -std_ea
hmat(26,33)=      -cgy*std_ea
hmat(25,34)=      -std_eb
hmat(26,35)=      -std_eg
hmat(27,36)=      -std_eqs
hmat(28,37)=      -std_em
hmat(30,38)=      -std_epinf
hmat(32,39)=      -std_ew
hmat(15,40)=-1.0000000000000000d0
hmat(40,40)=      -1+cikbar
hmat(41,41)=-0.99999900000000d0
hmat(6,49)=1.0000000000000000d0
hmat(7,49)=      1/csigma/(1+chabb/cgamma)-chabb/cgamma/csigma/(1+chabb/ &
cgamma)
hmat(2,50)=1.0000000000000000d0
hmat(4,50)=-1.0000000000000000d0
hmat(8,50)=      -crkky
hmat(1,51)=      -calfa
hmat(2,51)=      -1/czcap*(1-czcap)
hmat(3,51)=1.0000000000000000d0
hmat(3,52)=1.0000000000000000d0
hmat(4,52)=1.0000000000000000d0
hmat(9,52)=      -cfc*calfa
hmat(5,53)=1.0000000000000000d0
hmat(8,53)=      -ciy
hmat(11,53)=      -cikbar
hmat(5,54)=      -1/(1+cbetabar*cgamma)/cgamma**2/csadjcost
hmat(6,54)=1.0000000000000000d0
hmat(7,55)=1.0000000000000000d0
hmat(8,55)=      -ccy
hmat(10,55)=      -1/(1-chabb/cgamma)
hmat(8,56)=1.0000000000000000d0
hmat(9,56)=1.0000000000000000d0
hmat(23,56)=      cry-cry*crr+crdy
hmat(3,57)=-1.0000000000000000d0
hmat(7,57)=      -cwhlc/(1+chabb/cgamma)+cwhlc/csigma/(1+chabb/cgamma)
hmat(9,57)=      -cfc+cfc*calfa
hmat(10,57)=      -csigl
hmat(1,58)=      -1+calfa
hmat(3,58)=-1.0000000000000000d0
hmat(10,58)=1.0000000000000000d0
hmat(11,59)=1.0000000000000000d0
hmat(12,60)=1.0000000000000000d0
hmat(21,60)=      -1/(1+cbetabar*cgamma*cindp)/cprobp/((cfc-1)*curvp+1)+1/ &
(1+cbetabar*cgamma*cindp)*cbetabar*cgamma/((cfc-1)*curvp+1)+1/(1+ &
cbetabar*cgamma*cindp)/((cfc-1)*curvp+1)-1/(1+ &
cbetabar*cgamma*cindp)*cbetabar*cgamma*cprobp/((cfc-1)*curvp+1)
hmat(13,61)=1.0000000000000000d0
hmat(15,61)=-1.0000000000000000d0
hmat(19,61)=      -crkky
hmat(12,62)=      -calfa
hmat(13,62)=      -1/czcap*(1-czcap)
hmat(14,62)=1.0000000000000000d0
hmat(14,63)=1.0000000000000000d0
hmat(15,63)=1.0000000000000000d0
hmat(20,63)=      -cfc*calfa
hmat(16,64)=1.0000000000000000d0
hmat(19,64)=      -ciy
hmat(40,64)=      -cikbar
hmat(44,64)=-1.0000000000000000d0
hmat(16,65)=      -1/(1+cbetabar*cgamma)/cgamma**2/csadjcost
hmat(17,65)=1.0000000000000000d0
hmat(18,66)=1.0000000000000000d0
hmat(19,66)=      -ccy
hmat(22,66)=      -1/(1+cbetabar*cgamma)/cprobw/((clandaw-1)*curvw+1)/(1- &
chabb/cgamma)+1/(1+cbetabar*cgamma)*cbetabar*cgamma/((clandaw-1)*curvw+1)/ &
(1-chabb/cgamma)+1/(1+cbetabar*cgamma)/((clandaw-1)*curvw+1)/(1-chabb/ &
cgamma)-cprobw*cbetabar*cgamma/(1+cbetabar*cgamma)/((clandaw-1)*curvw+1)/ &
(1-chabb/cgamma)
hmat(43,66)=-1.0000000000000000d0
hmat(19,67)=1.0000000000000000d0
hmat(20,67)=1.0000000000000000d0
hmat(23,67)=      -cry+cry*crr-crdy
hmat(42,67)=-1.0000000000000000d0
hmat(14,68)=-1.0000000000000000d0
hmat(18,68)=      -cwhlc/(1+chabb/cgamma)+cwhlc/csigma/(1+chabb/cgamma)
hmat(20,68)=      -cfc+cfc*calfa
hmat(22,68)=      -1/(1+cbetabar*cgamma)/cprobw/((clandaw-1)*curvw+ &
1)*csigl+1/(1+cbetabar*cgamma)*cbetabar*cgamma/((clandaw-1)*curvw+1)*csigl+ &
1/(1+cbetabar*cgamma)/((clandaw-1)*curvw+1)*csigl-cprobw*cbetabar*cgamma/ &
(1+cbetabar*cgamma)/((clandaw-1)*curvw+1)*csigl
hmat(48,68)=-1.0000000000000000d0
hmat(21,69)=1.0000000000000000d0
hmat(22,69)=      1/(1+cbetabar*cgamma)+cbetabar*cgamma*cindw/(1+ &
cbetabar*cgamma)
hmat(23,69)=      -crpi+crpi*crr
hmat(46,69)=-1.0000000000000000d0
hmat(12,70)=      -1+calfa
hmat(14,70)=-1.0000000000000000d0
hmat(22,70)=      1+1/(1+cbetabar*cgamma)/cprobw/((clandaw-1)*curvw+1)-1/ &
(1+cbetabar*cgamma)*cbetabar*cgamma/((clandaw-1)*curvw+1)-1/(1+ &
cbetabar*cgamma)/((clandaw-1)*curvw+1)+cprobw*cbetabar*cgamma/(1+ &
cbetabar*cgamma)/((clandaw-1)*curvw+1)
hmat(45,70)=-1.0000000000000000d0
hmat(17,71)=1.0000000000000000d0
hmat(18,71)=      1/csigma/(1+chabb/cgamma)-chabb/cgamma/csigma/(1+chabb/ &
cgamma)
hmat(23,71)=1.0000000000000000d0
hmat(47,71)=-1.0000000000000000d0
hmat(1,72)=1.0000000000000000d0
hmat(9,72)=      -cfc
hmat(12,72)=1.0000000000000000d0
hmat(20,72)=      -cfc
hmat(24,72)=1.0000000000000000d0
hmat(6,73)=      -1/(1-chabb/cgamma)*csigma*(1+chabb/cgamma)
hmat(7,73)=-1.0000000000000000d0
hmat(17,73)=      -1/(1-chabb/cgamma)*csigma*(1+chabb/cgamma)
hmat(18,73)=-1.0000000000000000d0
hmat(25,73)=1.0000000000000000d0
hmat(8,74)=-1.0000000000000000d0
hmat(19,74)=-1.0000000000000000d0
hmat(26,74)=1.0000000000000000d0
hmat(5,75)=-1.0000000000000000d0
hmat(11,75)=      -cikbar*cgamma**2*csadjcost
hmat(16,75)=-1.0000000000000000d0
hmat(27,75)=1.0000000000000000d0
hmat(40,75)=      -cikbar*cgamma**2*csadjcost
hmat(23,76)=-1.0000000000000000d0
hmat(28,76)=1.0000000000000000d0
hmat(21,77)=-1.0000000000000000d0
hmat(29,77)=1.0000000000000000d0
hmat(29,78)=-1.0000000000000000d0
hmat(30,78)=1.0000000000000000d0
hmat(22,79)=-1.0000000000000000d0
hmat(31,79)=1.0000000000000000d0
hmat(31,80)=-1.0000000000000000d0
hmat(32,80)=1.0000000000000000d0
hmat(33,81)=1.0000000000000000d0
hmat(34,82)=1.0000000000000000d0
hmat(35,83)=1.0000000000000000d0
hmat(36,84)=1.0000000000000000d0
hmat(37,85)=1.0000000000000000d0
hmat(38,86)=1.0000000000000000d0
hmat(39,87)=1.0000000000000000d0
hmat(40,88)=1.0000000000000000d0
hmat(41,89)=1.0000000000000000d0
hmat(42,89)=      -ctrend
hmat(43,89)=      -ctrend
hmat(44,89)=      -ctrend
hmat(45,89)=      -ctrend
hmat(46,89)=      -constepinf
hmat(47,89)=      -conster
hmat(48,89)=      -constelab
hmat(42,90)=1.0000000000000000d0
hmat(43,91)=1.0000000000000000d0
hmat(44,92)=1.0000000000000000d0
hmat(45,93)=1.0000000000000000d0
hmat(46,94)=1.0000000000000000d0
hmat(47,95)=1.0000000000000000d0
hmat(48,96)=1.0000000000000000d0
hmat(6,99)=      -crk/(crk+1-ctou)
hmat(5,101)=      -1/(1+cbetabar*cgamma)*cbetabar*cgamma
hmat(6,102)=      -1/(crk+1-ctou)+ctou/(crk+1-ctou)
hmat(7,103)=      -1/(1+chabb/cgamma)
hmat(7,105)=      cwhlc/(1+chabb/cgamma)-cwhlc/csigma/(1+chabb/cgamma)
hmat(17,110)=      -crk/(crk+1-ctou)
hmat(16,112)=      -1/(1+cbetabar*cgamma)*cbetabar*cgamma
hmat(17,113)=      -1/(crk+1-ctou)+ctou/(crk+1-ctou)
hmat(18,114)=      -1/(1+chabb/cgamma)
hmat(18,116)=      cwhlc/(1+chabb/cgamma)-cwhlc/csigma/(1+chabb/cgamma)
hmat(17,117)=-1.0000000000000000d0
hmat(18,117)=      -1/csigma/(1+chabb/cgamma)+chabb/cgamma/csigma/(1+chabb/ &
cgamma)
hmat(21,117)=      -1/(1+cbetabar*cgamma*cindp)*cbetabar*cgamma
hmat(22,117)=      -1/(1+cbetabar*cgamma)*cbetabar*cgamma
hmat(22,118)=      -1/(1+cbetabar*cgamma)*cbetabar*cgamma


END SUBROUTINE
