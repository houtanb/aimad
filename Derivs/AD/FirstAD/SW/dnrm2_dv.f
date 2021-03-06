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
C  Differentiation of dnrm2 in forward (tangent) mode: (multi-directional mode)
C   variations  of output variables: dnrm2
C   with respect to input variables: x
      SUBROUTINE DNRM2_DV(n, x, xd, incx, dnrm2, dnrm2d, nbdirs)
      USE DIFFSIZES
      IMPLICIT NONE
C
C     End of DNRM2.
C
      INTEGER incx, n
      INTEGER nbdirs
      DOUBLE PRECISION dnrm2, dnrm2d(nbdirsmax)
      DOUBLE PRECISION x(*), xd(nbdirsmax, *)
      DOUBLE PRECISION one, zero
      PARAMETER (one=1.0d+0, zero=0.0d+0)
      DOUBLE PRECISION absxi, absxid(nbdirsmax), norm, normd(nbdirsmax)
     +                 , result1, result1d(nbdirsmax), scale, scaled(
     +                 nbdirsmax), ssq, ssqd(nbdirsmax)
      INTEGER ix, nd
      INTRINSIC ABS, SQRT
C     .. Scalar Arguments ..
C     ..
C     .. Array Arguments ..
C     ..
C
C  Purpose
C  =======
C
C  DNRM2 returns the euclidean norm of a vector via the function
C  name, so that
C
C     DNRM2 := sqrt( x'*x )
C
C
C  -- This version written on 25-October-1982.
C     Modified on 14-October-1993 to inline the call to DLASSQ.
C     Sven Hammarling, Nag Ltd.
C
C
C     .. Parameters ..
C     ..
C     .. Local Scalars ..
C     ..
C     .. Intrinsic Functions ..
C     ..
      IF (n .LT. 1 .OR. incx .LT. 1) THEN
        norm = zero
        DO nd=1,nbdirs
          normd(nd) = 0.D0
        ENDDO
      ELSE IF (n .EQ. 1) THEN
        IF (x(1) .GE. 0.) THEN
          DO nd=1,nbdirs
            normd(nd) = xd(nd, 1)
          ENDDO
          norm = x(1)
        ELSE
          DO nd=1,nbdirs
            normd(nd) = -xd(nd, 1)
          ENDDO
          norm = -x(1)
        END IF
      ELSE
        scale = zero
        ssq = one
        DO nd=1,nbdirs
          scaled(nd) = 0.D0
        ENDDO
        DO nd=1,nbdirs
          ssqd(nd) = 0.D0
        ENDDO
C        The following loop is equivalent to this call to the LAPACK
C        auxiliary routine:
C        CALL DLASSQ( N, X, INCX, SCALE, SSQ )
C
        DO ix=1,1+(n-1)*incx,incx
          IF (x(ix) .NE. zero) THEN
            IF (x(ix) .GE. 0.) THEN
              DO nd=1,nbdirs
                absxid(nd) = xd(nd, ix)
              ENDDO
              absxi = x(ix)
            ELSE
              DO nd=1,nbdirs
                absxid(nd) = -xd(nd, ix)
              ENDDO
              absxi = -x(ix)
            END IF
            IF (scale .LT. absxi) THEN
              DO nd=1,nbdirs
                ssqd(nd) = ssqd(nd)*scale**2/absxi**2 + ssq*2*scale*(
     +            scaled(nd)*absxi-scale*absxid(nd))/absxi**3
                scaled(nd) = absxid(nd)
              ENDDO
              ssq = one + ssq*(scale/absxi)**2
              scale = absxi
            ELSE
              DO nd=1,nbdirs
                ssqd(nd) = ssqd(nd) + 2*absxi*(absxid(nd)*scale-absxi*
     +            scaled(nd))/scale**3
              ENDDO
              ssq = ssq + (absxi/scale)**2
            END IF
          END IF
        ENDDO
        result1 = SQRT(ssq)
        DO nd=1,nbdirs
          IF (ssqd(nd) .EQ. 0.0 .OR. ssq .EQ. 0.0) THEN
            result1d(nd) = 0.D0
          ELSE
            result1d(nd) = ssqd(nd)/(2.0*SQRT(ssq))
          END IF
          normd(nd) = scaled(nd)*result1 + scale*result1d(nd)
        ENDDO
        norm = scale*result1
      END IF
      DO nd=1,nbdirs
        dnrm2d(nd) = normd(nd)
      ENDDO
C
      dnrm2 = norm
      RETURN
      END
