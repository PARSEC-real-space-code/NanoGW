##
#
# Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
# This file is part of RGWBS. It is distributed under the GPL v1.
#
##

##
## Common files
##
TARG	= absp chkpt_bin_asc proj_pol rotate_pol

tools: $(TARG)

VPATH	= ../shared:../utils:$(AUX_PATH)

.SUFFIXES: .f90 $(SUFFIXES)

.f90.o:
	$(FCPP) $< > $*.tmp.f90
	$(F90s) -c $(OPTS) $*.tmp.f90 -o $@

absp: absp.o
	$(F90s) $(OPTS) -o ../bin/absp absp.o

chkpt_bin_asc: chkpt_bin_asc.o
	$(F90s) $(OPTS) -o ../bin/chkpt_bin_asc chkpt_bin_asc.o

proj_pol: proj_pol.o quicksort.o
	$(F90s) $(OPTS) -o ../bin/proj_pol proj_pol.o quicksort.o

rotate_pol: rotate_pol.o get_date.o quicksort.o
	$(F90s) $(OPTS) -o ../bin/rotate_pol rotate_pol.o get_date.o \
	quicksort.o

