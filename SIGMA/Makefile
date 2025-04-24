##
#
# Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
# This file is part of RGWBS. It is distributed under the GPL v1.
#
##

##
SRC	= alccheck.f90 define_kmap.f90 define_xc.f90 delete_file.f90 \
	  get_date.f90 get_grad_sym.f90 get_lap_sym.f90 input_g.f90 \
	  inverse_map.f90 map_check.f90 mapinverse.f90 mtrxin.f90 \
	  newtonr.f90 paratec_atoms.f90 paratec_xc.f90 parsec_atoms.f90 \
	  parsec_xc.f90 printsigma_nogw.f90 qpt_matching.f90 quicksort.f90 \
	  read_ptgroups.f90 read_scissors.f90 read_sigma.f90 \
	  select1234.f90 spline.f90 \
	  spline_evaluate.f90 splint.f90 stopwatch.f90 sum_rule.f90 \
	  timacc.f90 unfold.f90 vector_check.f90 \
	  write_sigma.f90 linear_solver.f90 printmatrix.f90 \
	  lanczos_spectra.f90 apply_matvec.f90\
      lanczos_spectra_UltraLowMem.f90 apply_matvec_UltraLowMem.f90 \
	  lanczos_spectra_BLOCK.f90 apply_matvec_BLOCK.f90 \
	  polyfit.f90 polyfitcheb.f90 $(AUX_SRC)
PSRCZ = calculate_tdlda.F90z calculate_sigma.F90z \
        calculate_sigma_lanczos.F90z \
        calculate_sigma_lanczos_UltraLowMem.F90z \
        calculate_sigma_lanczos_BLOCK.F90z \
	    calculate_sigma_lanczos_lowcomsym.F90z \
	    charac_group.F90z  contract_en.F90z correlation.F90z  \
	    redistribute_psi_loc.F90z \
	    define_pmap.F90z diag_pol.F90z dipole.F90z eigensolver.F90z eigensolver_new.F90z energy_mtxel.F90z \
	    etotal.F90z fock_exchange.F90z gw_correlation.F90z eigensolver_hipMAGMA.F90z\
	    get_pol.F90z k_integrate.F90z get_pol_isdf.F90z \
	    k_print.F90z kernel.F90z matvec3.F90z model_exchange.F90z model_correlation.F90z multiply_vec.F90z \
	    potential.F90z print_qp.F90z proj_calc.F90z psum.F90z \
	    read_pol.F90z rotate_qp.F90z save_pol.F90z setup_g.F90z setup_s.F90z static.F90z \
	    static_corr.F90z vertex.F90z vxc.F90z gemm_hl.F90z \
	    vxcmtxel.F90z wpol_v.F90z wpol0.F90z write_pol.F90z \
	    k_integrate_isdf.F90z wpol0_isdf.F90z k_integrate_isdf_fastselect.F90z \
	    isdf_parallel_sym_UltraLowMem.F90z lanczos_spectra_lowcomsym.F90z apply_matvec_lowcomsym.F90z
# submodule.F90z
PSRC = die.F90 grid_setup.F90 finalize.F90 input_s.F90 paratec_wfn.F90 \
	   parsec_wfn.F90 printsigma.F90 qp_eigensolver.F90 \
	   save_sigma.F90 sigma.F90 isdf_parallel_sym.F90 cvt_wfn.F90 cvt_parallel_sym1.F90 \
	   cvt_parallel_sym2.F90 cvt.F90 isdf_parallel_sym_lessmemory2.F90 int_psum.F90
# CSRC	= hadamard_prod_mat_sym.cpp

OBJ	= $(SRC:.f90=.o)
POBJ	= $(PSRC:.F90=.o)
POBJZ	= $(PSRCZ:.F90z=.o)
#COBJ	= ${CSRC:.cpp=.o}

MODULE	= $(AUX_MOD) myconstants.o typedefs.o esdf_key_module.o esdf_module.o \
	  xc_module.o mpi_module.o fft_module.o fd_module.o psp_module.o 
# submodule.o

INITIAL	= init_tmp.o

VPATH	= ../shared:../SIGMA:$(AUX_PATH)

#sigma: module obj initial pobjz pobj cobj
#	$(F90) $(OPTS) $(OBJ) $(POBJ) $(POBJZ) $(COBJ) \
#	$(MODULE) $(INITIAL) $(LIBFFT) $(LIBLAPACK) $(LIBSCALAPACK) $(LIBXC) \
#	$(LIBHDF5) $(LIBMAGMA) $(LIBHIP) $(LIBHIP) -o ../bin/sigma$(EXT) 

sigma: module obj initial pobjz pobj 
	$(F90) $(OBJ) $(POBJ) $(POBJZ) \
	$(MODULE) $(INITIAL) $(LIBFFT) $(LIBLAPACK) $(LIBSCALAPACK) $(LIBXC) \
	$(LIBHDF5) $(LIBMAGMA) $(LIBHIP) $(LIBHIP) $(OPTS) -o ../bin/sigma$(EXT) 

.SUFFIXES: .f90 .F90 .F90z $(SUFFIXES)

obj: module
	make FOPTS="$(OPTS)" $(OBJ)

pobj: module
	make FOPTS="$(OPTS)" $(POBJ)

pobjz: module
	make FOPTS="$(OPTS)" $(POBJZ)

#cobj:
#	make $(COBJ)
	
module:
	make FOPTS="$(OPTS)" $(MODULE)

initial:
	sed "s/MACH/$(MACH)/g" ../shared/header.F90 \
	| sed "s/CPPOPT/$(CPPOPT)/g" \
	| sed "s/OPTS/$(OPTS)/g" \
	| sed "s/GIT_COMMIT/$(GIT_COMMIT)/g" \
	| sed "s/BUILD_DATE/$(BUILD_DATE)/g" > init_tmp.F90
	make FOPTS="$(OPTS)" $(INITIAL)

$(POBJ): $(MODULE)
$(POBJZ): $(MODULE)

.f90.o:
	$(FCPP) $< > $*.tmp.f90
	$(F90s) -c $(OPTS) $(OPTS2) $*.tmp.f90 -o $@

.F90.o:
	$(FCPP) $(CPPOPT) $< > $*.tmp.f90
	$(F90) -c $(OPTS) $(OPTS2) $*.tmp.f90 -o $@

.F90z.o:
	$(FCPP) $(CPPOPT) $< > $*.tmp.f90
	$(FCPP) $(CPPOPT) -DCPLX $< >> $*.tmp.f90
	$(F90) -c $(OPTS) $(OPTS2) $*.tmp.f90 -o $@

#.cpp.o:
#	cp $< $*.cpp
#	$(HIPCC) -c $*.cpp -o $@

