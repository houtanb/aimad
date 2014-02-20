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
C  Differentiation of ddot in forward (tangent) mode: (multi-directional mode)
C   variations  of output variables: ddot
C   with respect to input variables: dx dy
      SUBROUTINE DDOT_DV(n, dx, dxd, incx, dy, dyd, incy, ddot, ddotd,
     +                   nbdirs)
      USE DIFFSIZES
      IMPLICIT NONE
      INTEGER incx, incy, n
      INTEGER nbdirs
      DOUBLE PRECISION ddot, ddotd(nbdirsmax)
      DOUBLE PRECISION dx(*), dxd(nbdirsmax, *), dy(*), dyd(nbdirsmax, *
     +                 )
      DOUBLE PRECISION dtemp, dtempd(nbdirsmax)
      INTEGER i, ix, iy, m, mp1, nd
      INTRINSIC MOD
C     .. Scalar Arguments ..
C     ..
C     .. Array Arguments ..
C     ..
C
C  Purpose
C  =======
C
C     forms the dot product of two vectors.
C     uses unrolled loops for increments equal to one.
C     jack dongarra, linpack, 3/11/78.
C     modified 12/3/93, array(1) declarations changed to array(*)
C
C
C     .. Local Scalars ..
C     ..
C     .. Intrinsic Functions ..
C     ..
      ddot = 0.0d0
      dtemp = 0.0d0
      IF (n .LE. 0) THEN
        DO nd=1,nbdirs
          ddotd(nd) = 0.D0
        ENDDO
        RETURN
      ELSE IF (incx .EQ. 1 .AND. incy .EQ. 1) THEN
C
C        code for both increments equal to 1
C
C
C        clean-up loop
C
        m = MOD(n, 5)
        IF (m .EQ. 0) THEN
          DO nd=1,nbdirs
            dtempd(nd) = 0.D0
          ENDDO
        ELSE
          DO nd=1,nbdirs
            dtempd(nd) = 0.D0
          ENDDO
          DO i=1,m
            DO nd=1,nbdirs
              dtempd(nd) = dtempd(nd) + dxd(nd, i)*dy(i) + dx(i)*dyd(nd
     +          , i)
            ENDDO
            dtemp = dtemp + dx(i)*dy(i)
          ENDDO
          IF (n .LT. 5) GOTO 60
        END IF
        mp1 = m + 1
        DO i=mp1,n,5
          DO nd=1,nbdirs
            dtempd(nd) = dtempd(nd) + dxd(nd, i)*dy(i) + dx(i)*dyd(nd, i
     +        ) + dxd(nd, i+1)*dy(i+1) + dx(i+1)*dyd(nd, i+1) + dxd(nd,
     +        i+2)*dy(i+2) + dx(i+2)*dyd(nd, i+2) + dxd(nd, i+3)*dy(i+3)
     +        + dx(i+3)*dyd(nd, i+3) + dxd(nd, i+4)*dy(i+4) + dx(i+4)*
     +        dyd(nd, i+4)
          ENDDO
          dtemp = dtemp + dx(i)*dy(i) + dx(i+1)*dy(i+1) + dx(i+2)*dy(i+2
     +      ) + dx(i+3)*dy(i+3) + dx(i+4)*dy(i+4)
        ENDDO
 60     DO nd=1,nbdirs
          ddotd(nd) = dtempd(nd)
        ENDDO
        ddot = dtemp
        RETURN
      ELSE
C
C        code for unequal increments or equal increments
C          not equal to 1
C
        ix = 1
        iy = 1
        IF (incx .LT. 0) ix = (-n+1)*incx + 1
        IF (incy .LT. 0) THEN
          iy = (-n+1)*incy + 1
          DO nd=1,nbdirs
            dtempd(nd) = 0.D0
          ENDDO
        ELSE
          DO nd=1,nbdirs
            dtempd(nd) = 0.D0
          ENDDO
        END IF
        DO i=1,n
          DO nd=1,nbdirs
            dtempd(nd) = dtempd(nd) + dxd(nd, ix)*dy(iy) + dx(ix)*dyd(nd
     +        , iy)
          ENDDO
          dtemp = dtemp + dx(ix)*dy(iy)
          ix = ix + incx
          iy = iy + incy
        ENDDO
        DO nd=1,nbdirs
          ddotd(nd) = dtempd(nd)
        ENDDO
        ddot = dtemp
        RETURN
      END IF
      END