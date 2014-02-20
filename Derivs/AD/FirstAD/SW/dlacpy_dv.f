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
C  Differentiation of dlacpy in forward (tangent) mode: (multi-directional mode)
C   variations  of output variables: b
C   with respect to input variables: a b
      SUBROUTINE DLACPY_DV(uplo, m, n, a, ad, lda, b, bd, ldb, nbdirs)
      USE DIFFSIZES
      IMPLICIT NONE
C
C     End of DLACPY
C
      INTEGER lda, ldb, nbdirs
      DOUBLE PRECISION a(lda, *), ad(nbdirsmax, lda, *), b(ldb, *), bd(
     +                 nbdirsmax, ldb, *)
      INTEGER m, n
      CHARACTER uplo
      INTEGER i, j, min1, nd
      LOGICAL LSAME, result1
      EXTERNAL LSAME
      INTRINSIC MIN
C
C  -- LAPACK auxiliary routine (version 3.1) --
C     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
C     November 2006
C
C     .. Scalar Arguments ..
C     ..
C     .. Array Arguments ..
C     ..
C
C  Purpose
C  =======
C
C  DLACPY copies all or part of a two-dimensional matrix A to another
C  matrix B.
C
C  Arguments
C  =========
C
C  UPLO    (input) CHARACTER*1
C          Specifies the part of the matrix A to be copied to B.
C          = 'U':      Upper triangular part
C          = 'L':      Lower triangular part
C          Otherwise:  All of the matrix A
C
C  M       (input) INTEGER
C          The number of rows of the matrix A.  M >= 0.
C
C  N       (input) INTEGER
C          The number of columns of the matrix A.  N >= 0.
C
C  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
C          The m by n matrix A.  If UPLO = 'U', only the upper triangle
C          or trapezoid is accessed; if UPLO = 'L', only the lower
C          triangle or trapezoid is accessed.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= max(1,M).
C
C  B       (output) DOUBLE PRECISION array, dimension (LDB,N)
C          On exit, B = A in the locations specified by UPLO.
C
C  LDB     (input) INTEGER
C          The leading dimension of the array B.  LDB >= max(1,M).
C
C  =====================================================================
C
C     .. Local Scalars ..
C     ..
C     .. External Functions ..
C     ..
C     .. Intrinsic Functions ..
C     ..
C     .. Executable Statements ..
C
      result1 = LSAME(uplo, 'U')
      IF (result1) THEN
        DO j=1,n
          IF (j .GT. m) THEN
            min1 = m
          ELSE
            min1 = j
          END IF
          DO i=1,min1
            DO nd=1,nbdirs
              bd(nd, i, j) = ad(nd, i, j)
            ENDDO
            b(i, j) = a(i, j)
          ENDDO
        ENDDO
      ELSE
        result1 = LSAME(uplo, 'L')
        IF (result1) THEN
          DO j=1,n
            DO i=j,m
              DO nd=1,nbdirs
                bd(nd, i, j) = ad(nd, i, j)
              ENDDO
              b(i, j) = a(i, j)
            ENDDO
          ENDDO
        ELSE
          DO j=1,n
            DO i=1,m
              DO nd=1,nbdirs
                bd(nd, i, j) = ad(nd, i, j)
              ENDDO
              b(i, j) = a(i, j)
            ENDDO
          ENDDO
        END IF
      END IF
      RETURN
      END