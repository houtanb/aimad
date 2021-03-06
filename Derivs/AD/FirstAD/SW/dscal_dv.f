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

C        Generated by TAPENADE     (INRIA, Tropics team)
C  Tapenade - Version 2.2 (r1239) - Wed 28 Jun 2006 04:59:55 PM CEST
C
C  Differentiation of dscal in forward (tangent) mode: (multi-directional mode)
C   variations  of output variables: dx
C   with respect to input variables: dx da
      SUBROUTINE DSCAL_DV(n, da, dad, dx, dxd, incx, nbdirs)
      USE DIFFSIZES
      IMPLICIT NONE
      INTEGER incx, n
      INTEGER nbdirs
      DOUBLE PRECISION da, dad(nbdirsmax)
      DOUBLE PRECISION dx(*), dxd(nbdirsmax, *)
      INTEGER i, m, mp1, nd, nincx
      INTRINSIC MOD
C     .. Scalar Arguments ..
C     ..
C     .. Array Arguments ..
C     ..
C
C  Purpose
C  =======
C*
C     scales a vector by a constant.
C     uses unrolled loops for increment equal to one.
C     jack dongarra, linpack, 3/11/78.
C     modified 3/93 to return if incx .le. 0.
C     modified 12/3/93, array(1) declarations changed to array(*)
C
C
C     .. Local Scalars ..
C     ..
C     .. Intrinsic Functions ..
C     ..
      IF (n .LE. 0 .OR. incx .LE. 0) THEN
        RETURN
      ELSE IF (incx .EQ. 1) THEN
C
C        code for increment equal to 1
C
C
C        clean-up loop
C
        m = MOD(n, 5)
        IF (.NOT.m .EQ. 0) THEN
          DO i=1,m
            DO nd=1,nbdirs
              dxd(nd, i) = dad(nd)*dx(i) + da*dxd(nd, i)
            ENDDO
            dx(i) = da*dx(i)
          ENDDO
          IF (n .LT. 5) RETURN
        END IF
        mp1 = m + 1
        DO i=mp1,n,5
          DO nd=1,nbdirs
            dxd(nd, i) = dad(nd)*dx(i) + da*dxd(nd, i)
          ENDDO
          dx(i) = da*dx(i)
          DO nd=1,nbdirs
            dxd(nd, i+1) = dad(nd)*dx(i+1) + da*dxd(nd, i+1)
          ENDDO
          dx(i+1) = da*dx(i+1)
          DO nd=1,nbdirs
            dxd(nd, i+2) = dad(nd)*dx(i+2) + da*dxd(nd, i+2)
          ENDDO
          dx(i+2) = da*dx(i+2)
          DO nd=1,nbdirs
            dxd(nd, i+3) = dad(nd)*dx(i+3) + da*dxd(nd, i+3)
          ENDDO
          dx(i+3) = da*dx(i+3)
          DO nd=1,nbdirs
            dxd(nd, i+4) = dad(nd)*dx(i+4) + da*dxd(nd, i+4)
          ENDDO
          dx(i+4) = da*dx(i+4)
        ENDDO
        RETURN
      ELSE
C
C        code for increment not equal to 1
C
        nincx = n*incx
        DO i=1,nincx,incx
          DO nd=1,nbdirs
            dxd(nd, i) = dad(nd)*dx(i) + da*dxd(nd, i)
          ENDDO
          dx(i) = da*dx(i)
        ENDDO
        RETURN
      END IF
      END
