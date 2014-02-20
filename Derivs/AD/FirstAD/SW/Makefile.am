
bin_PROGRAMS = firstad

firstad_LDADD = ../../libraries/LIBAIM2/libaim.a ../../libraries/lapack-3.2.1/SRC/liblapack.a ../../libraries/lapack-3.2.1/BLAS/SRC/libblas.a

firstad_SOURCES=DIFFSIZES.f90 modminimizer.f90 \
	daxpy_dv.f   dgehrd_dv.f  dgetrs_dv.f  dlahr2_dv.f  dlaqr2_dv.f  dlarft_dv.f  dlasy2_dv.f  dswap_dv.f \
	dcopy_dv.f   dgemm_dv.f   dhseqr_dv.f  dlaln2_dv.f  dlaqr3_dv.f  dlarfx_dv.f  dnrm2_dv.f   dtrexc_dv.f \
	ddot_dv.f    dgemv_dv.f   dlacn2_dv.f  dlange_dv.f  dlaqr4_dv.f  dlartg_dv.f  dorg2r_dv.f  dtrmm_dv.f \
	dgebak_dv.f  dger_dv.f    dlacpy_dv.f  dlanv2_dv.f  dlaqr5_dv.f  dlascl_dv.f  dorghr_dv.f  dtrmv_dv.f \
	dgebal_dv.f  dgesv_dv.f   dladiv_dv.f  dlapy2_dv.f  dlarfb_dv.f  dlaset_dv.f  dorgqr_dv.f  dtrsen_dv.f \
	dgeesx_dv.f  dgetf2_dv.f  dlaexc_dv.f  dlaqr0_dv.f  dlarf_dv.f   dlassq_dv.f  drot_dv.f    dtrsm_dv.f \
	dgehd2_dv.f  dgetrf_dv.f  dlahqr_dv.f  dlaqr1_dv.f  dlarfg_dv.f  dlaswp_dv.f  dscal_dv.f   dtrsyl_dv.f \
	aim_dv.f90               eyemat.f90               loaddata.f90                  strmatch.f90 \
	beginfile.f90            eyematnod.f90            matrixmult_dv.f90             svd_analytical_deriv.f90 \
	build_a_dv.f90           get_hmat_dv.f90          moderror.f90                  svd_gary_analytical_1deriv.f90 \
	call_likel_dv.f90        get_selector.f90         svd_gary_analytical_deriv.f90 \
	constructPiMat.f90       get_stable_singvals.f90  numeric_shift_dv.f90          transp_dv.f90 \
	copy_w_dv.f90            get_zerocols.f90         writemat.f90 \
	get_zerorows.f90         objfun_likel.f90              writematint.f90 \
	eigensystem_analytd.f90  inverse_d.f90            pseudoinv.f90                 writematint_total.f90 \
	eigensystem_dv.f90       invert_dv.f90            reduced_form_dv.f90           writemat_total_3d.f90 \
	exact_shift_dv.f90       kalman_dv.f90            selectf.f90                   writemat_total.f90 \
	existstrue.f90           kroneckerProduct.f90     shiftright_dv.f90 \
	eyemat_dv.f90            likelihood_dv.f90        splitmat_dv.f90
