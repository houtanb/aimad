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
C  Differentiation of dhseqr in forward (tangent) mode: (multi-directional mode)
C   variations  of output variables: wr h z work wi
C   with respect to input variables: wr h z work wi
      SUBROUTINE DHSEQR_DV(job, compz, n, ilo, ihi, h, hd, ldh, wr, wrd
     +                     , wi, wid, z, zd, ldz, work, workd, lwork,
     +                     info, nbdirs)
      USE DIFFSIZES
      IMPLICIT NONE
C
C     ==== End of DHSEQR ====
C
      CHARACTER compz, job
      INTEGER ihi, ilo, info, lwork, n
      INTEGER ldh, ldz, nbdirs
      DOUBLE PRECISION h(ldh, *), hd(nbdirsmax, ldh, *), wi(*), wid(
     +                 nbdirsmax, *), work(*), workd(nbdirsmax, *), wr(*
     +                 ), wrd(nbdirsmax, *), z(ldz, *), zd(nbdirsmax,
     +                 ldz, *)
      INTEGER nl
      PARAMETER (nl=49)
      INTEGER ntiny
      PARAMETER (ntiny=11)
      DOUBLE PRECISION one, zero
      PARAMETER (one=1.0d0, zero=0.0d0)
      DOUBLE PRECISION hl(nl, nl), hld(nbdirsmax, nl, nl), workl(nl),
     +                 workld(nbdirsmax, nl)
      INTEGER arg1, arg2, i, ii2, ii3, ILAENV, kbot, max1, max2, max3,
     +        max4, max5, max6, max7, min1, nd, nmin
      LOGICAL initz, lquery, LSAME, result1, wantt, wantz
      DOUBLE PRECISION x1, x2
      EXTERNAL ILAENV, XERBLA, LSAME
      INTRINSIC MAX, DBLE, MIN
C
C  -- LAPACK driver routine (version 3.1) --
C     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
C     November 2006
C
C     .. Scalar Arguments ..
C     ..
C     .. Array Arguments ..
C     ..
C     Purpose
C     =======
C
C     DHSEQR computes the eigenvalues of a Hessenberg matrix H
C     and, optionally, the matrices T and Z from the Schur decomposition
C     H = Z T Z**T, where T is an upper quasi-triangular matrix (the
C     Schur form), and Z is the orthogonal matrix of Schur vectors.
C
C     Optionally Z may be postmultiplied into an input orthogonal
C     matrix Q so that this routine can give the Schur factorization
C     of a matrix A which has been reduced to the Hessenberg form H
C     by the orthogonal matrix Q:  A = Q*H*Q**T = (QZ)*T*(QZ)**T.
C
C     Arguments
C     =========
C
C     JOB   (input) CHARACTER*1
C           = 'E':  compute eigenvalues only;
C           = 'S':  compute eigenvalues and the Schur form T.
C
C     COMPZ (input) CHARACTER*1
C           = 'N':  no Schur vectors are computed;
C           = 'I':  Z is initialized to the unit matrix and the matrix Z
C                   of Schur vectors of H is returned;
C           = 'V':  Z must contain an orthogonal matrix Q on entry, and
C                   the product Q*Z is returned.
C
C     N     (input) INTEGER
C           The order of the matrix H.  N .GE. 0.
C
C     ILO   (input) INTEGER
C     IHI   (input) INTEGER
C           It is assumed that H is already upper triangular in rows
C           and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
C           set by a previous call to DGEBAL, and then passed to DGEHRD
C           when the matrix output by DGEBAL is reduced to Hessenberg
C           form. Otherwise ILO and IHI should be set to 1 and N
C           respectively.  If N.GT.0, then 1.LE.ILO.LE.IHI.LE.N.
C           If N = 0, then ILO = 1 and IHI = 0.
C
C     H     (input/output) DOUBLE PRECISION array, dimension (LDH,N)
C           On entry, the upper Hessenberg matrix H.
C           On exit, if INFO = 0 and JOB = 'S', then H contains the
C           upper quasi-triangular matrix T from the Schur decomposition
C           (the Schur form); 2-by-2 diagonal blocks (corresponding to
C           complex conjugate pairs of eigenvalues) are returned in
C           standard form, with H(i,i) = H(i+1,i+1) and
C           H(i+1,i)*H(i,i+1).LT.0. If INFO = 0 and JOB = 'E', the
C           contents of H are unspecified on exit.  (The output value of
C           H when INFO.GT.0 is given under the description of INFO
C           below.)
C
C           Unlike earlier versions of DHSEQR, this subroutine may
C           explicitly H(i,j) = 0 for i.GT.j and j = 1, 2, ... ILO-1
C           or j = IHI+1, IHI+2, ... N.
C
C     LDH   (input) INTEGER
C           The leading dimension of the array H. LDH .GE. max(1,N).
C
C     WR    (output) DOUBLE PRECISION array, dimension (N)
C     WI    (output) DOUBLE PRECISION array, dimension (N)
C           The real and imaginary parts, respectively, of the computed
C           eigenvalues. If two eigenvalues are computed as a complex
C           conjugate pair, they are stored in consecutive elements of
C           WR and WI, say the i-th and (i+1)th, with WI(i) .GT. 0 and
C           WI(i+1) .LT. 0. If JOB = 'S', the eigenvalues are stored in
C           the same order as on the diagonal of the Schur form returned
C           in H, with WR(i) = H(i,i) and, if H(i:i+1,i:i+1) is a 2-by-2
C           diagonal block, WI(i) = sqrt(-H(i+1,i)*H(i,i+1)) and
C           WI(i+1) = -WI(i).
C
C     Z     (input/output) DOUBLE PRECISION array, dimension (LDZ,N)
C           If COMPZ = 'N', Z is not referenced.
C           If COMPZ = 'I', on entry Z need not be set and on exit,
C           if INFO = 0, Z contains the orthogonal matrix Z of the Schur
C           vectors of H.  If COMPZ = 'V', on entry Z must contain an
C           N-by-N matrix Q, which is assumed to be equal to the unit
C           matrix except for the submatrix Z(ILO:IHI,ILO:IHI). On exit,
C           if INFO = 0, Z contains Q*Z.
C           Normally Q is the orthogonal matrix generated by DORGHR
C           after the call to DGEHRD which formed the Hessenberg matrix
C           H. (The output value of Z when INFO.GT.0 is given under
C           the description of INFO below.)
C
C     LDZ   (input) INTEGER
C           The leading dimension of the array Z.  if COMPZ = 'I' or
C           COMPZ = 'V', then LDZ.GE.MAX(1,N).  Otherwize, LDZ.GE.1.
C
C     WORK  (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
C           On exit, if INFO = 0, WORK(1) returns an estimate of
C           the optimal value for LWORK.
C
C     LWORK (input) INTEGER
C           The dimension of the array WORK.  LWORK .GE. max(1,N)
C           is sufficient, but LWORK typically as large as 6*N may
C           be required for optimal performance.  A workspace query
C           to determine the optimal workspace size is recommended.
C
C           If LWORK = -1, then DHSEQR does a workspace query.
C           In this case, DHSEQR checks the input parameters and
C           estimates the optimal workspace size for the given
C           values of N, ILO and IHI.  The estimate is returned
C           in WORK(1).  No error message related to LWORK is
C           issued by XERBLA.  Neither H nor Z are accessed.
C
C
C     INFO  (output) INTEGER
C             =  0:  successful exit
C           .LT. 0:  if INFO = -i, the i-th argument had an illegal
C                    value
C           .GT. 0:  if INFO = i, DHSEQR failed to compute all of
C                the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR
C                and WI contain those eigenvalues which have been
C                successfully computed.  (Failures are rare.)
C
C                If INFO .GT. 0 and JOB = 'E', then on exit, the
C                remaining unconverged eigenvalues are the eigen-
C                values of the upper Hessenberg matrix rows and
C                columns ILO through INFO of the final, output
C                value of H.
C
C                If INFO .GT. 0 and JOB   = 'S', then on exit
C
C           (*)  (initial value of H)*U  = U*(final value of H)
C
C                where U is an orthogonal matrix.  The final
C                value of H is upper Hessenberg and quasi-triangular
C                in rows and columns INFO+1 through IHI.
C
C                If INFO .GT. 0 and COMPZ = 'V', then on exit
C
C                  (final value of Z)  =  (initial value of Z)*U
C
C                where U is the orthogonal matrix in (*) (regard-
C                less of the value of JOB.)
C
C                If INFO .GT. 0 and COMPZ = 'I', then on exit
C                      (final value of Z)  = U
C                where U is the orthogonal matrix in (*) (regard-
C                less of the value of JOB.)
C
C                If INFO .GT. 0 and COMPZ = 'N', then Z is not
C                accessed.
C
C     ================================================================
C             Default values supplied by
C             ILAENV(ISPEC,'DHSEQR',JOB(:1)//COMPZ(:1),N,ILO,IHI,LWORK).
C             It is suggested that these defaults be adjusted in order
C             to attain best performance in each particular
C             computational environment.
C
C            ISPEC=1:  The DLAHQR vs DLAQR0 crossover point.
C                      Default: 75. (Must be at least 11.)
C
C            ISPEC=2:  Recommended deflation window size.
C                      This depends on ILO, IHI and NS.  NS is the
C                      number of simultaneous shifts returned
C                      by ILAENV(ISPEC=4).  (See ISPEC=4 below.)
C                      The default for (IHI-ILO+1).LE.500 is NS.
C                      The default for (IHI-ILO+1).GT.500 is 3*NS/2.
C
C            ISPEC=3:  Nibble crossover point. (See ILAENV for
C                      details.)  Default: 14% of deflation window
C                      size.
C
C            ISPEC=4:  Number of simultaneous shifts, NS, in
C                      a multi-shift QR iteration.
C
C                      If IHI-ILO+1 is ...
C
C                      greater than      ...but less    ... the
C                      or equal to ...      than        default is
C
C                           1               30          NS -   2(+)
C                          30               60          NS -   4(+)
C                          60              150          NS =  10(+)
C                         150              590          NS =  **
C                         590             3000          NS =  64
C                        3000             6000          NS = 128
C                        6000             infinity      NS = 256
C
C                  (+)  By default some or all matrices of this order
C                       are passed to the implicit double shift routine
C                       DLAHQR and NS is ignored.  See ISPEC=1 above
C                       and comments in IPARM for details.
C
C                       The asterisks (**) indicate an ad-hoc
C                       function of N increasing from 10 to 64.
C
C            ISPEC=5:  Select structured matrix multiply.
C                      (See ILAENV for details.) Default: 3.
C
C     ================================================================
C     Based on contributions by
C        Karen Braman and Ralph Byers, Department of Mathematics,
C        University of Kansas, USA
C
C     ================================================================
C     References:
C       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
C       Algorithm Part I: Maintaining Well Focused Shifts, and Level 3
C       Performance, SIAM Journal of Matrix Analysis, volume 23, pages
C       929--947, 2002.
C
C       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
C       Algorithm Part II: Aggressive Early Deflation, SIAM Journal
C       of Matrix Analysis, volume 23, pages 948--973, 2002.
C
C     ================================================================
C     .. Parameters ..
C
C     ==== Matrices of order NTINY or smaller must be processed by
C     .    DLAHQR because of insufficient subdiagonal scratch space.
C     .    (This is a hard limit.) ====
C
C     ==== NL allocates some local workspace to help small matrices
C     .    through a rare DLAHQR failure.  NL .GT. NTINY = 11 is
C     .    required and NL .LE. NMIN = ILAENV(ISPEC=1,...) is recom-
C     .    mended.  (The default value of NMIN is 75.)  Using NL = 49
C     .    allows up to six simultaneous shifts and a 16-by-16
C     .    deflation window.  ====
C
C     ..
C     .. Local Arrays ..
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
C     ==== Decode and check the input parameters. ====
C
      wantt = LSAME(job, 'S')
      initz = LSAME(compz, 'I')
      result1 = LSAME(compz, 'V')
      wantz = initz .OR. result1
      IF (1 .LT. n) THEN
        max1 = n
      ELSE
        max1 = 1
      END IF
      DO nd=1,nbdirs
        workd(nd, 1) = 0.D0
      ENDDO
      work(1) = DBLE(max1)
      lquery = lwork .EQ. -1
C
      info = 0
      result1 = LSAME(job, 'E')
      IF (.NOT.result1 .AND. (.NOT.wantt)) THEN
        info = -1
      ELSE
        result1 = LSAME(compz, 'N')
        IF (.NOT.result1 .AND. (.NOT.wantz)) THEN
          info = -2
        ELSE IF (n .LT. 0) THEN
          info = -3
        ELSE
          IF (1 .LT. n) THEN
            max2 = n
          ELSE
            max2 = 1
          END IF
          IF (ilo .LT. 1 .OR. ilo .GT. max2) THEN
            info = -4
          ELSE
            IF (ilo .GT. n) THEN
              min1 = n
            ELSE
              min1 = ilo
            END IF
            IF (ihi .LT. min1 .OR. ihi .GT. n) THEN
              info = -5
            ELSE
              IF (1 .LT. n) THEN
                max3 = n
              ELSE
                max3 = 1
              END IF
              IF (ldh .LT. max3) THEN
                info = -7
              ELSE
                IF (1 .LT. n) THEN
                  max4 = n
                ELSE
                  max4 = 1
                END IF
                IF (ldz .LT. 1 .OR. wantz .AND. ldz .LT. max4) THEN
                  info = -11
                ELSE
                  IF (1 .LT. n) THEN
                    max5 = n
                  ELSE
                    max5 = 1
                  END IF
                  IF (lwork .LT. max5 .AND. (.NOT.lquery)) info = -13
                END IF
              END IF
            END IF
          END IF
        END IF
      END IF
C
      IF (info .NE. 0) THEN
C
C        ==== Quick return in case of invalid argument. ====
C
        arg1 = -info
        CALL XERBLA('DHSEQR', arg1)
        RETURN
C
      ELSE IF (n .EQ. 0) THEN
C
C        ==== Quick return in case N = 0; nothing to do. ====
C
        RETURN
C
      ELSE IF (lquery) THEN
C
C        ==== Quick return in case of a workspace query ====
C
        CALL DLAQR0_DV(wantt, wantz, n, ilo, ihi, h, hd, ldh, wr, wrd,
     +                 wi, wid, ilo, ihi, z, zd, ldz, work, workd, lwork
     +                 , info, nbdirs)
        IF (1 .LT. n) THEN
          max6 = n
        ELSE
          max6 = 1
        END IF
        x1 = DBLE(max6)
        IF (x1 .LT. work(1)) THEN
          work(1) = work(1)
        ELSE
          DO nd=1,nbdirs
            workd(nd, 1) = 0.D0
          ENDDO
          work(1) = x1
        END IF
        RETURN
C
      ELSE
C
C        ==== copy eigenvalues isolated by DGEBAL ====
C
        DO i=1,ilo-1
          DO nd=1,nbdirs
            wrd(nd, i) = hd(nd, i, i)
            wid(nd, i) = 0.D0
          ENDDO
          wr(i) = h(i, i)
          wi(i) = zero
        ENDDO
        DO i=ihi+1,n
          DO nd=1,nbdirs
            wrd(nd, i) = hd(nd, i, i)
            wid(nd, i) = 0.D0
          ENDDO
          wr(i) = h(i, i)
          wi(i) = zero
        ENDDO
C
C        ==== Initialize Z, if requested ====
C
        IF (initz) CALL DLASET_DV('A', n, n, zero, one, z, zd, ldz,
     +                            nbdirs)
C
C        ==== Quick return if possible ====
C
        IF (ilo .EQ. ihi) THEN
          DO nd=1,nbdirs
            wrd(nd, ilo) = hd(nd, ilo, ilo)
            wid(nd, ilo) = 0.D0
          ENDDO
          wr(ilo) = h(ilo, ilo)
          wi(ilo) = zero
          RETURN
        ELSE
C
C        ==== DLAHQR/DLAQR0 crossover point ====
C
          nmin = ILAENV(12, 'DHSEQR', job(:1)//compz(:1), n, ilo, ihi,
     +      lwork)
          IF (ntiny .LT. nmin) THEN
            nmin = nmin
          ELSE
            nmin = ntiny
          END IF
C
C        ==== DLAQR0 for big matrices; DLAHQR for small ones ====
C
          IF (n .GT. nmin) THEN
            CALL DLAQR0_DV(wantt, wantz, n, ilo, ihi, h, hd, ldh, wr,
     +                     wrd, wi, wid, ilo, ihi, z, zd, ldz, work,
     +                     workd, lwork, info, nbdirs)
          ELSE
C
C           ==== Small matrix ====
C
            CALL DLAHQR_DV(wantt, wantz, n, ilo, ihi, h, hd, ldh, wr,
     +                     wrd, wi, wid, ilo, ihi, z, zd, ldz, info,
     +                     nbdirs)
C
            IF (info .GT. 0) THEN
C
C              ==== A rare DLAHQR failure!  DLAQR0 sometimes succeeds
C              .    when DLAHQR fails. ====
C
              kbot = info
C
              IF (n .GE. nl) THEN
C
C                 ==== Larger matrices have enough subdiagonal scratch
C                 .    space to call DLAQR0 directly. ====
C
                CALL DLAQR0_DV(wantt, wantz, n, ilo, kbot, h, hd, ldh,
     +                         wr, wrd, wi, wid, ilo, ihi, z, zd, ldz,
     +                         work, workd, lwork, info, nbdirs)
C
              ELSE
                DO nd=1,nbdirs
                  DO ii2=1,49
                    DO ii3=1,49
                      hld(nd, ii3, ii2) = 0.D0
                    ENDDO
                  ENDDO
                  DO ii2=1,49
                    workld(nd, ii2) = 0.D0
                  ENDDO
                ENDDO
C
C                 ==== Tiny matrices don't have enough subdiagonal
C                 .    scratch space to benefit from DLAQR0.  Hence,
C                 .    tiny matrices must be copied into a larger
C                 .    array before calling DLAQR0. ====
C
                CALL DLACPY_DV('A', n, n, h, hd, ldh, hl, hld, nl,
     +                         nbdirs)
                DO nd=1,nbdirs
                  hld(nd, n+1, n) = 0.D0
                ENDDO
                hl(n+1, n) = zero
                arg1 = nl - n
                CALL DLASET_DV('A', nl, arg1, zero, zero, hl(1, n+1),
     +                         hld(1, 1, n+1), nl, nbdirs)
                CALL DLAQR0_DV(wantt, wantz, nl, ilo, kbot, hl, hld, nl
     +                         , wr, wrd, wi, wid, ilo, ihi, z, zd, ldz
     +                         , workl, workld, nl, info, nbdirs)
                IF (wantt .OR. info .NE. 0) CALL DLACPY_DV('A', n, n, hl
     +                                                     , hld, nl, h
     +                                                     , hd, ldh,
     +                                                     nbdirs)
              END IF
            END IF
          END IF
C
C        ==== Clear out the trash, if necessary. ====
C
          IF ((wantt .OR. info .NE. 0) .AND. n .GT. 2) THEN
            arg1 = n - 2
            arg2 = n - 2
            CALL DLASET_DV('L', arg1, arg2, zero, zero, h(3, 1), hd(1, 3
     +                     , 1), ldh, nbdirs)
          END IF
          IF (1 .LT. n) THEN
            max7 = n
          ELSE
            max7 = 1
          END IF
          x2 = DBLE(max7)
          IF (x2 .LT. work(1)) THEN
            work(1) = work(1)
          ELSE
            DO nd=1,nbdirs
              workd(nd, 1) = 0.D0
            ENDDO
            work(1) = x2
          END IF
        END IF
      END IF
      END