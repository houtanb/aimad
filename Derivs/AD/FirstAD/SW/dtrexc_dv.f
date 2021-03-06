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
C  Differentiation of dtrexc in forward (tangent) mode: (multi-directional mode)
C   variations  of output variables: q t work
C   with respect to input variables: q t work
      SUBROUTINE DTREXC_DV(compq, n, t, td, ldt, q, qd, ldq, ifst, ilst
     +                     , work, workd, info, nbdirs)
      USE DIFFSIZES
      IMPLICIT NONE
C
C     End of DTREXC
C
      CHARACTER compq
      INTEGER ifst, ilst, info, n
      INTEGER ldq, ldt, nbdirs
      DOUBLE PRECISION q(ldq, *), qd(nbdirsmax, ldq, *), t(ldt, *), td(
     +                 nbdirsmax, ldt, *), work(*), workd(nbdirsmax, *)
      DOUBLE PRECISION zero
      PARAMETER (zero=0.0d+0)
      INTEGER arg1, here, max1, max2, nbf, nbl, nbnext
      LOGICAL LSAME, result1, wantq
      INTRINSIC MAX
      EXTERNAL XERBLA, LSAME
C
C  -- LAPACK routine (version 3.1) --
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
C  DTREXC reorders the real Schur factorization of a real matrix
C  A = Q*T*Q**T, so that the diagonal block of T with row index IFST is
C  moved to row ILST.
C
C  The real Schur form T is reordered by an orthogonal similarity
C  transformation Z**T*T*Z, and optionally the matrix Q of Schur vectors
C  is updated by postmultiplying it with Z.
C
C  T must be in Schur canonical form (as returned by DHSEQR), that is,
C  block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each
C  2-by-2 diagonal block has its diagonal elements equal and its
C  off-diagonal elements of opposite sign.
C
C  Arguments
C  =========
C
C  COMPQ   (input) CHARACTER*1
C          = 'V':  update the matrix Q of Schur vectors;
C          = 'N':  do not update Q.
C
C  N       (input) INTEGER
C          The order of the matrix T. N >= 0.
C
C  T       (input/output) DOUBLE PRECISION array, dimension (LDT,N)
C          On entry, the upper quasi-triangular matrix T, in Schur
C          Schur canonical form.
C          On exit, the reordered upper quasi-triangular matrix, again
C          in Schur canonical form.
C
C  LDT     (input) INTEGER
C          The leading dimension of the array T. LDT >= max(1,N).
C
C  Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N)
C          On entry, if COMPQ = 'V', the matrix Q of Schur vectors.
C          On exit, if COMPQ = 'V', Q has been postmultiplied by the
C          orthogonal transformation matrix Z which reorders T.
C          If COMPQ = 'N', Q is not referenced.
C
C  LDQ     (input) INTEGER
C          The leading dimension of the array Q.  LDQ >= max(1,N).
C
C  IFST    (input/output) INTEGER
C  ILST    (input/output) INTEGER
C          Specify the reordering of the diagonal blocks of T.
C          The block with row index IFST is moved to row ILST, by a
C          sequence of transpositions between adjacent blocks.
C          On exit, if IFST pointed on entry to the second row of a
C          2-by-2 block, it is changed to point to the first row; ILST
C          always points to the first row of the block in its final
C          position (which may differ from its input value by +1 or -1).
C          1 <= IFST <= N; 1 <= ILST <= N.
C
C  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
C
C  INFO    (output) INTEGER
C          = 0:  successful exit
C          < 0:  if INFO = -i, the i-th argument had an illegal value
C          = 1:  two adjacent blocks were too close to swap (the problem
C                is very ill-conditioned); T may have been partially
C                reordered, and ILST points to the first row of the
C                current position of the block being moved.
C
C  =====================================================================
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
C     .. Executable Statements ..
C
C     Decode and test the input arguments.
C
      info = 0
      wantq = LSAME(compq, 'V')
      result1 = LSAME(compq, 'N')
      IF (.NOT.wantq .AND. (.NOT.result1)) THEN
        info = -1
      ELSE IF (n .LT. 0) THEN
        info = -2
      ELSE
        IF (1 .LT. n) THEN
          max1 = n
        ELSE
          max1 = 1
        END IF
        IF (ldt .LT. max1) THEN
          info = -4
        ELSE
          IF (1 .LT. n) THEN
            max2 = n
          ELSE
            max2 = 1
          END IF
          IF (ldq .LT. 1 .OR. wantq .AND. ldq .LT. max2) THEN
            info = -6
          ELSE IF (ifst .LT. 1 .OR. ifst .GT. n) THEN
            info = -7
          ELSE IF (ilst .LT. 1 .OR. ilst .GT. n) THEN
            info = -8
          END IF
        END IF
      END IF
      IF (info .NE. 0) THEN
        arg1 = -info
        CALL XERBLA('DTREXC', arg1)
        RETURN
      ELSE IF (n .LE. 1) THEN
C
C     Quick return if possible
C
        RETURN
      ELSE
C
C     Determine the first row of specified block
C     and find out it is 1 by 1 or 2 by 2.
C
        IF (ifst .GT. 1) THEN
          IF (t(ifst, ifst-1) .NE. zero) ifst = ifst - 1
        END IF
        nbf = 1
        IF (ifst .LT. n) THEN
          IF (t(ifst+1, ifst) .NE. zero) nbf = 2
        END IF
C
C     Determine the first row of the final block
C     and find out it is 1 by 1 or 2 by 2.
C
        IF (ilst .GT. 1) THEN
          IF (t(ilst, ilst-1) .NE. zero) ilst = ilst - 1
        END IF
        nbl = 1
        IF (ilst .LT. n) THEN
          IF (t(ilst+1, ilst) .NE. zero) nbl = 2
        END IF
C
        IF (ifst .EQ. ilst) THEN
          RETURN
        ELSE
C
          IF (ifst .LT. ilst) THEN
C
C        Update ILST
C
            IF (nbf .EQ. 2 .AND. nbl .EQ. 1) ilst = ilst - 1
            IF (nbf .EQ. 1 .AND. nbl .EQ. 2) ilst = ilst + 1
C
            here = ifst
C
C
C        Swap block with next one below
C
 10         IF (nbf .EQ. 1 .OR. nbf .EQ. 2) THEN
C
C           Current block either 1 by 1 or 2 by 2
C
              nbnext = 1
              IF (here + nbf + 1 .LE. n) THEN
                IF (t(here+nbf+1, here+nbf) .NE. zero) nbnext = 2
              END IF
              CALL DLAEXC_DV(wantq, n, t, td, ldt, q, qd, ldq, here, nbf
     +                       , nbnext, work, workd, info, nbdirs)
              IF (info .NE. 0) THEN
                GOTO 100
              ELSE
                here = here + nbnext
C
C           Test if 2 by 2 block breaks into two 1 by 1 blocks
C
                IF (nbf .EQ. 2) THEN
                  IF (t(here+1, here) .EQ. zero) nbf = 3
                END IF
              END IF
            ELSE
C
C
C           Current block consists of two 1 by 1 blocks each of which
C           must be swapped individually
C
              nbnext = 1
              IF (here + 3 .LE. n) THEN
                IF (t(here+3, here+2) .NE. zero) nbnext = 2
              END IF
              arg1 = here + 1
              CALL DLAEXC_DV(wantq, n, t, td, ldt, q, qd, ldq, arg1, 1,
     +                       nbnext, work, workd, info, nbdirs)
              IF (info .NE. 0) THEN
                GOTO 120
              ELSE IF (nbnext .EQ. 1) THEN
C
C              Swap two 1 by 1 blocks, no problems possible
C
                CALL DLAEXC_DV(wantq, n, t, td, ldt, q, qd, ldq, here, 1
     +                         , nbnext, work, workd, info, nbdirs)
                here = here + 1
              ELSE
C
C              Recompute NBNEXT in case 2 by 2 split
C
                IF (t(here+2, here+1) .EQ. zero) nbnext = 1
                IF (nbnext .EQ. 2) THEN
C
C                 2 by 2 Block did not split
C
                  CALL DLAEXC_DV(wantq, n, t, td, ldt, q, qd, ldq, here
     +                           , 1, nbnext, work, workd, info, nbdirs)
                  IF (info .NE. 0) THEN
                    GOTO 110
                  ELSE
                    here = here + 2
                  END IF
                ELSE
C
C                 2 by 2 Block did split
C
                  CALL DLAEXC_DV(wantq, n, t, td, ldt, q, qd, ldq, here
     +                           , 1, 1, work, workd, info, nbdirs)
                  arg1 = here + 1
                  CALL DLAEXC_DV(wantq, n, t, td, ldt, q, qd, ldq, arg1
     +                           , 1, 1, work, workd, info, nbdirs)
                  here = here + 2
                END IF
              END IF
            END IF
            IF (here .LT. ilst) THEN
              GOTO 10
            ELSE
              GOTO 160
            END IF
 100        ilst = here
            RETURN
 110        ilst = here
            RETURN
 120        ilst = here
            RETURN
          ELSE
C
C
            here = ifst
C
C        Swap block with next one above
C
 20         IF (nbf .EQ. 1 .OR. nbf .EQ. 2) THEN
C
C           Current block either 1 by 1 or 2 by 2
C
              nbnext = 1
              IF (here .GE. 3) THEN
                IF (t(here-1, here-2) .NE. zero) nbnext = 2
              END IF
              arg1 = here - nbnext
              CALL DLAEXC_DV(wantq, n, t, td, ldt, q, qd, ldq, arg1,
     +                       nbnext, nbf, work, workd, info, nbdirs)
              IF (info .NE. 0) THEN
                GOTO 130
              ELSE
                here = here - nbnext
C
C           Test if 2 by 2 block breaks into two 1 by 1 blocks
C
                IF (nbf .EQ. 2) THEN
                  IF (t(here+1, here) .EQ. zero) nbf = 3
                END IF
              END IF
            ELSE
C
C
C           Current block consists of two 1 by 1 blocks each of which
C           must be swapped individually
C
              nbnext = 1
              IF (here .GE. 3) THEN
                IF (t(here-1, here-2) .NE. zero) nbnext = 2
              END IF
              arg1 = here - nbnext
              CALL DLAEXC_DV(wantq, n, t, td, ldt, q, qd, ldq, arg1,
     +                       nbnext, 1, work, workd, info, nbdirs)
              IF (info .NE. 0) THEN
                GOTO 150
              ELSE IF (nbnext .EQ. 1) THEN
C
C              Swap two 1 by 1 blocks, no problems possible
C
                CALL DLAEXC_DV(wantq, n, t, td, ldt, q, qd, ldq, here,
     +                         nbnext, 1, work, workd, info, nbdirs)
                here = here - 1
              ELSE
C
C              Recompute NBNEXT in case 2 by 2 split
C
                IF (t(here, here-1) .EQ. zero) nbnext = 1
                IF (nbnext .EQ. 2) THEN
C
C                 2 by 2 Block did not split
C
                  arg1 = here - 1
                  CALL DLAEXC_DV(wantq, n, t, td, ldt, q, qd, ldq, arg1
     +                           , 2, 1, work, workd, info, nbdirs)
                  IF (info .NE. 0) THEN
                    GOTO 140
                  ELSE
                    here = here - 2
                  END IF
                ELSE
C
C                 2 by 2 Block did split
C
                  CALL DLAEXC_DV(wantq, n, t, td, ldt, q, qd, ldq, here
     +                           , 1, 1, work, workd, info, nbdirs)
                  arg1 = here - 1
                  CALL DLAEXC_DV(wantq, n, t, td, ldt, q, qd, ldq, arg1
     +                           , 1, 1, work, workd, info, nbdirs)
                  here = here - 2
                END IF
              END IF
            END IF
            IF (here .GT. ilst) THEN
              GOTO 20
            ELSE
              GOTO 160
            END IF
 130        ilst = here
            RETURN
 140        ilst = here
            RETURN
 150        ilst = here
            RETURN
          END IF
 160      ilst = here
C
          RETURN
        END IF
      END IF
      END
