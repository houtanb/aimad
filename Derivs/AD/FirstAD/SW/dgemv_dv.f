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
C  Differentiation of dgemv in forward (tangent) mode: (multi-directional mode)
C   variations  of output variables: y
C   with respect to input variables: alpha x y a
      SUBROUTINE DGEMV_DV(trans, m, n, alpha, alphad, a, ad, lda, x, xd
     +                    , incx, beta, y, yd, incy, nbdirs)
      USE DIFFSIZES
      IMPLICIT NONE
C
C     End of DGEMV .
C
      INTEGER incx, incy, m, n
      INTEGER lda, nbdirs
      DOUBLE PRECISION a(lda, *), ad(nbdirsmax, lda, *), x(*), xd(
     +                 nbdirsmax, *), y(*), yd(nbdirsmax, *)
      DOUBLE PRECISION alpha, alphad(nbdirsmax), beta
      CHARACTER trans
      DOUBLE PRECISION one, zero
      PARAMETER (one=1.0d+0, zero=0.0d+0)
      INTEGER i, info, ix, iy, j, jx, jy, kx, ky, lenx, leny, max1, nd
      LOGICAL LSAME, result1, result2, result3
      DOUBLE PRECISION temp, tempd(nbdirsmax)
      INTRINSIC MAX
      EXTERNAL XERBLA, LSAME
C     .. Scalar Arguments ..
C     ..
C     .. Array Arguments ..
C     ..
C
C  Purpose
C  =======
C
C  DGEMV  performs one of the matrix-vector operations
C
C     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
C
C  where alpha and beta are scalars, x and y are vectors and A is an
C  m by n matrix.
C
C  Arguments
C  ==========
C
C  TRANS  - CHARACTER*1.
C           On entry, TRANS specifies the operation to be performed as
C           follows:
C
C              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
C
C              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
C
C              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C           On entry, M specifies the number of rows of the matrix A.
C           M must be at least zero.
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the number of columns of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - DOUBLE PRECISION.
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
C           Before entry, the leading m by n part of the array A must
C           contain the matrix of coefficients.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, m ).
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of DIMENSION at least
C           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
C           and at least
C           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
C           Before entry, the incremented array X must contain the
C           vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  BETA   - DOUBLE PRECISION.
C           On entry, BETA specifies the scalar beta. When BETA is
C           supplied as zero then Y need not be set on input.
C           Unchanged on exit.
C
C  Y      - DOUBLE PRECISION array of DIMENSION at least
C           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
C           and at least
C           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
C           Before entry with BETA non-zero, the incremented array Y
C           must contain the vector y. On exit, Y is overwritten by the
C           updated vector y.
C
C  INCY   - INTEGER.
C           On entry, INCY specifies the increment for the elements of
C           Y. INCY must not be zero.
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
C     ..
C     .. Local Scalars ..
C     ..
C     .. External Functions ..
C     ..
C     .. External Subroutines ..
C     ..
C     .. Intrinsic Functions ..
C     ..
C
C     Test the input parameters.
C
      info = 0
      result1 = LSAME(trans, 'N')
      result2 = LSAME(trans, 'T')
      result3 = LSAME(trans, 'C')
      IF (.NOT.result1 .AND. (.NOT.result2) .AND. (.NOT.result3)) THEN
        info = 1
      ELSE IF (m .LT. 0) THEN
        info = 2
      ELSE IF (n .LT. 0) THEN
        info = 3
      ELSE
        IF (1 .LT. m) THEN
          max1 = m
        ELSE
          max1 = 1
        END IF
        IF (lda .LT. max1) THEN
          info = 6
        ELSE IF (incx .EQ. 0) THEN
          info = 8
        ELSE IF (incy .EQ. 0) THEN
          info = 11
        END IF
      END IF
      IF (info .NE. 0) THEN
        CALL XERBLA('DGEMV ', info)
        RETURN
      ELSE IF (m .EQ. 0 .OR. n .EQ. 0 .OR. alpha .EQ. zero .AND. beta
     +    .EQ. one) THEN
C
C     Quick return if possible.
C
        RETURN
      ELSE
C
C     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
C     up the start points in  X  and  Y.
C
        result1 = LSAME(trans, 'N')
        IF (result1) THEN
          lenx = n
          leny = m
        ELSE
          lenx = m
          leny = n
        END IF
        IF (incx .GT. 0) THEN
          kx = 1
        ELSE
          kx = 1 - (lenx-1)*incx
        END IF
        IF (incy .GT. 0) THEN
          ky = 1
        ELSE
          ky = 1 - (leny-1)*incy
        END IF
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through A.
C
C     First form  y := beta*y.
C
        IF (beta .NE. one) THEN
          IF (incy .EQ. 1) THEN
            IF (beta .EQ. zero) THEN
              DO i=1,leny
                DO nd=1,nbdirs
                  yd(nd, i) = 0.D0
                ENDDO
                y(i) = zero
              ENDDO
            ELSE
              DO i=1,leny
                DO nd=1,nbdirs
                  yd(nd, i) = beta*yd(nd, i)
                ENDDO
                y(i) = beta*y(i)
              ENDDO
            END IF
          ELSE
            iy = ky
            IF (beta .EQ. zero) THEN
              DO i=1,leny
                DO nd=1,nbdirs
                  yd(nd, iy) = 0.D0
                ENDDO
                y(iy) = zero
                iy = iy + incy
              ENDDO
            ELSE
              DO i=1,leny
                DO nd=1,nbdirs
                  yd(nd, iy) = beta*yd(nd, iy)
                ENDDO
                y(iy) = beta*y(iy)
                iy = iy + incy
              ENDDO
            END IF
          END IF
        END IF
        IF (alpha .EQ. zero) THEN
          RETURN
        ELSE
          result1 = LSAME(trans, 'N')
          IF (result1) THEN
C
C        Form  y := alpha*A*x + y.
C
            jx = kx
            IF (incy .EQ. 1) THEN
              DO j=1,n
                IF (x(jx) .NE. zero) THEN
                  DO nd=1,nbdirs
                    tempd(nd) = alphad(nd)*x(jx) + alpha*xd(nd, jx)
                  ENDDO
                  temp = alpha*x(jx)
                  DO i=1,m
                    DO nd=1,nbdirs
                      yd(nd, i) = yd(nd, i) + tempd(nd)*a(i, j) + temp*
     +                  ad(nd, i, j)
                    ENDDO
                    y(i) = y(i) + temp*a(i, j)
                  ENDDO
                END IF
                jx = jx + incx
              ENDDO
            ELSE
              DO j=1,n
                IF (x(jx) .NE. zero) THEN
                  DO nd=1,nbdirs
                    tempd(nd) = alphad(nd)*x(jx) + alpha*xd(nd, jx)
                  ENDDO
                  temp = alpha*x(jx)
                  iy = ky
                  DO i=1,m
                    DO nd=1,nbdirs
                      yd(nd, iy) = yd(nd, iy) + tempd(nd)*a(i, j) + temp
     +                  *ad(nd, i, j)
                    ENDDO
                    y(iy) = y(iy) + temp*a(i, j)
                    iy = iy + incy
                  ENDDO
                END IF
                jx = jx + incx
              ENDDO
            END IF
          ELSE
C
C        Form  y := alpha*A'*x + y.
C
            jy = ky
            IF (incx .EQ. 1) THEN
              DO j=1,n
                temp = zero
                DO nd=1,nbdirs
                  tempd(nd) = 0.D0
                ENDDO
                DO i=1,m
                  DO nd=1,nbdirs
                    tempd(nd) = tempd(nd) + ad(nd, i, j)*x(i) + a(i, j)*
     +                xd(nd, i)
                  ENDDO
                  temp = temp + a(i, j)*x(i)
                ENDDO
                DO nd=1,nbdirs
                  yd(nd, jy) = yd(nd, jy) + alphad(nd)*temp + alpha*
     +              tempd(nd)
                ENDDO
                y(jy) = y(jy) + alpha*temp
                jy = jy + incy
              ENDDO
            ELSE
              DO j=1,n
                temp = zero
                ix = kx
                DO nd=1,nbdirs
                  tempd(nd) = 0.D0
                ENDDO
                DO i=1,m
                  DO nd=1,nbdirs
                    tempd(nd) = tempd(nd) + ad(nd, i, j)*x(ix) + a(i, j)
     +                *xd(nd, ix)
                  ENDDO
                  temp = temp + a(i, j)*x(ix)
                  ix = ix + incx
                ENDDO
                DO nd=1,nbdirs
                  yd(nd, jy) = yd(nd, jy) + alphad(nd)*temp + alpha*
     +              tempd(nd)
                ENDDO
                y(jy) = y(jy) + alpha*temp
                jy = jy + incy
              ENDDO
            END IF
          END IF
C
          RETURN
        END IF
      END IF
      END