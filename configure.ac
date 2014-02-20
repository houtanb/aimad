AC_PREREQ([2.61])
AC_INIT([AIM2AD], [1.0-unstable])
AC_CONFIG_SRCDIR([configure.ac])
AM_INIT_AUTOMAKE([-Wall -Werror foreign 1.13])

AC_CANONICAL_HOST

AC_PROG_CC
AC_PROG_FC
AM_PROG_AR
AC_PROG_F77([gfortran g77 f77])
AM_CONDITIONAL([HAVE_FORT], [test "x$F77" != "x"])

AC_CHECK_PROGS([AR], [ar ${host_alias}-ar])
AC_PROG_RANLIB

FFLAGS="$FFLAGS -Wall"
FCFLAGS="$FCFLAGS -Wall"

BUILD_EVERYTHING="yup"

if test "xAR" == "x"; then
  BUILD_EVERYTHING="nope (ar not found)"
fi

if test "xHAVE_FORT" == "x"; then
  BUILD_EVERYTHING="nope (fortran compiler not found)"
fi

AC_MSG_NOTICE([

****************************************
*
* AIMAD
*   Will build proceed correctly: $BUILD_EVERYTHING
*
****************************************
])


AC_OUTPUT([Makefile
           Derivs/AD/libraries/lapack-3.2.1/BLAS/SRC/Makefile
           Derivs/AD/libraries/lapack-3.2.1/SRC/Makefile
           Derivs/AD/libraries/LIBAIM2/Makefile
           Derivs/AD/FirstAD/SW/Makefile
           Derivs/FD/SW/Makefile
])
