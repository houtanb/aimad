
bin_PROGRAMS = orig finite
noinst_LIBRARIES = libaim.a

FILES=loadcof.f90 exact_shift.f90 get_zerorows.f90 get_selector.f90 \
	shiftright.f90 numeric_shiftsvd.f90 constructq.f90 transp.f90 matrixmult.f90 \
	get_stable.f90 build_a.f90 eigensystem.f90 sort.f get_zerocols.f90 \
	reduced_form.f90 copy_w.f90 aim.f90 printmat.f90 printmatz.f90 writemat.f90 \
	selectf.f90  timing.c get_stable_singvals.f90 invert.f90 numeric_shift.f90

orig_LDADD = ../lapack-3.2.1/SRC/liblapack.a ../lapack-3.2.1/BLAS/SRC/libblas.a
finite_LDADD = ../lapack-3.2.1/SRC/liblapack.a ../lapack-3.2.1/BLAS/SRC/libblas.a

orig_SOURCES = \
	call_aim.f90 \
	$(FILES)

finite_SOURCES = \
	finite_diff.f90 \
	$(FILES)

nodist_libaim_a_SOURCES = \
	$(FILES)
