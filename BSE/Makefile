##
#
# Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
# This file is part of RGWBS. It is distributed under the GPL v1.
#
##

##
SRC	= alccheck.f90 define_kmap.f90 define_xc.f90 \
	  delete_file.f90 get_date.f90 \
	  get_grad_sym.f90 get_lap_sym.f90 input_b.f90 input_g.f90 \
	  map_check.f90 mapinverse.f90 mtrxin.f90 paratec_atoms.f90 \
	  paratec_xc.f90 parsec_atoms.f90 parsec_xc.f90 qpt_matching.f90 \
	  quicksort.f90 read_ptgroups.f90 read_scissors.f90 select1234.f90 \
	  stopwatch.f90 sum_rule.f90 \
	  timacc.f90 unfold.f90 $(AUX_SRC)
PSRCZ	= charac_group.F90z calculate_bse.F90z calculate_tdlda.F90z \
	  define_pmap.F90z diag_bse.F90z diag_bse_mix.F90z diag_pol.F90z dipole.F90z \
	  direct_b.F90z direct_mix.F90z direct_s.F90z eigensolver.F90z eigensolver_new.F90z\
	  exchange_b.F90z get_pol.F90z eigensolver_hipMAGMA.F90z \
	  group_reduce_bse.F90z k_integrate.F90z k_print.F90z kernel.F90z \
          k_integrate_isdf.F90z k_integrate_isdf_fastselect.F90z \
	  matvec3.F90z multiply_vec.F90z potential.F90z \
	  proj_calc.F90z psum.F90z read_bse.F90z read_pol.F90z \
	  save_pol.F90z setup_g.F90z setup_b.F90z write_bse.F90z write_pol.F90z
PSRC	= die.F90 finalize.F90 grid_setup.F90 paratec_wfn.F90 \
	  parsec_wfn.F90 quasiparticle.F90 \
	  bsesolv.F90

OBJ	= $(SRC:.f90=.o)
POBJ	= $(PSRC:.F90=.o)
POBJZ	= $(PSRCZ:.F90z=.o)

MODULE	= $(AUX_MOD) myconstants.o typedefs.o esdf_key_module.o esdf_module.o \
	  xc_module.o mpi_module.o fft_module.o fd_module.o psp_module.o

INITIAL	= init_tmp.o

VPATH	= ../shared:../BSE:$(AUX_PATH)

bsesolv: module obj initial pobjz pobj
	$(F90) $(OPTS) -o ../bin/bsesolv$(EXT) $(OBJ) $(POBJ) $(POBJZ) \
	$(MODULE) $(INITIAL) $(LIBFFT) $(LIBLAPACK) $(LIBSCALAPACK) $(LIBXC) $(LIBMAGMA)

.SUFFIXES: .f90 .F90 .F90z $(SUFFIXES)

obj: module
	make FOPTS="$(OPTS)" $(OBJ)

pobj: module
	make FOPTS="$(OPTS)" $(POBJ)

pobjz: module
	make FOPTS="$(OPTS)" $(POBJZ)

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
	$(F90s) -c $(OPTS) $*.tmp.f90 -o $@

.F90.o:
	$(FCPP) $(CPPOPT) $< > $*.tmp.f90
	$(F90) -c $(OPTS) $*.tmp.f90 -o $@

.F90z.o:
	$(FCPP) $(CPPOPT) $< > $*.tmp.f90
	$(FCPP) $(CPPOPT) -DCPLX $< >> $*.tmp.f90
	$(F90) -c $(OPTS) $*.tmp.f90 -o $@
