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

SUBROUTINE matrixmult(a,nrows,ncols,x,nxcols,y)

INTEGER, INTENT(IN) :: nrows, ncols, nxcols
DOUBLE PRECISION, INTENT(IN), DIMENSION(nrows,ncols) :: a
DOUBLE PRECISION, INTENT(IN), DIMENSION(ncols,nxcols) :: x
DOUBLE PRECISION, INTENT(OUT), DIMENSION(nrows,nxcols) :: y

INTEGER :: indxc
INTEGER :: indxyc, indxyr
DOUBLE PRECISION :: factor


DO indxyr=1,nrows
   DO indxyc=1,nxcols
      factor = 0.d0
      DO indxc =1,ncols
         factor = factor + a(indxyr,indxc)*x(indxc,indxyc)
      END DO
      y(indxyr,indxyc) = factor
   END DO
END DO

!y = MATMUL(a,x)

END SUBROUTINE
