#####################
#
# Author:       Murilo Tiago, Univ. of Minnesota (2004)
#               mtiago@msi.umn.edu
#
# Copyright (C) 2009 Murilo L. Tiago, http://users.ices.utexas.edu/~mtiago
# This file is part of RGWBS. It is distributed under the GPL v1.
#
#####################
#  Make sure the MACH paramter reflects the correct platform.
#  Existing parameter files for different platforms are found under
#  subdirectory ./config
#####################

#MACH = sugon_kunshan_dcu
MACH = perlmutter_nvdia
#MACH = perlmutter_intel
#MACH = linux_nv
#MACH  = linux_nv_openblas
#MACH = linux_nv_mkl

# Version date: must have exactly 28 characters!
VERSIONDATE = "Wed Aug 23 17:07:49 EST 2023"

#####################
# Typically, no changes are needed beyond this point
#####################

MACH_CONF = ./config/make.$(MACH).h

include $(MACH_CONF)

default:
	@echo 
	@echo "Using compilation parameters for machine $(MACH)"
	@echo " if this is not your machine, edit ./Makefile file"
	@echo
	@echo " Pre-compilation options: $(CPPOPT)"
	@echo " executable extension: $(EXT)"
	@echo
	@echo "  COMPILATION OPTIONS:"
	@echo " make tdlda     - code for calculation of TDLDA excited states"
	@echo " make sigma     - code for self-energy in GW approximation"
	@echo " make bsesolv   - code for solution of Bethe-Salpeter equation"
	@echo " make w_blip  - code for BLIP interpolation" 
	@echo " make util      - compile pre-processing and post-processing tools"
	@echo " make all       - tdlda+sigma+bsesolv+w_blip+util"
	@echo " make clean     - delete all generated files"
	@echo " make allclean  - also delete binaries and backup files"
	@echo
	@echo " Note: after failure, use make clean to create a clean slate."
	@echo

DEPTH	= .

UTILS	= absp chkpt_bin_asc proj_pol rotate_pol

#
all: tdlda sigma bsesolv w_blip util

tdlda: create_tmp
	@ cd tmp_$(MACH) ; \
	/bin/rm -f mach_tmp ; echo "MACH  = " $(MACH) > mach_tmp ; \
	echo "VERSIONDATE = " $(VERSIONDATE)  >> mach_tmp ; \
	/bin/rm -f Makefile ; cat ../$(MACH_CONF) mach_tmp \
	../TDLDA/Makefile > makefile ; $(MAKE) -f makefile

sigma: create_tmp
	@ cd tmp_$(MACH) ; \
	/bin/rm -f mach_tmp ; echo "MACH  = " $(MACH) > mach_tmp ; \
	echo "VERSIONDATE = " $(VERSIONDATE) >> mach_tmp ; \
	/bin/rm -f Makefile ; cat ../$(MACH_CONF) mach_tmp \
	../SIGMA/Makefile > makefile ; $(MAKE) -f makefile

bsesolv: create_tmp
	@ cd tmp_$(MACH) ; \
	/bin/rm -f mach_tmp ; echo "MACH  = " $(MACH) > mach_tmp ; \
	echo "VERSIONDATE = " $(VERSIONDATE) >> mach_tmp ; \
	/bin/rm -f Makefile ; cat ../$(MACH_CONF) mach_tmp \
	../BSE/Makefile > makefile ; $(MAKE) -f makefile

w_blip: create_tmp
	@ cd tmp_$(MACH) ; \
	/bin/rm -f mach_tmp ; echo "MACH  = " $(MACH) > mach_tmp ; \
	echo "VERSIONDATE = " $(VERSIONDATE) >> mach_tmp ; \
	/bin/rm -f Makefile ; cat ../$(MACH_CONF) mach_tmp \
	../W_BLIP/Makefile > makefile ; $(MAKE) -f makefile

util: $(UTILS)

$(UTILS): create_tmp
	@ cd tmp_$(MACH) ; /bin/rm -f Makefile ; \
	cat ../$(MACH_CONF) ../utils/Makefile > \
	makefile ; $(MAKE) -f makefile -e $@

create_tmp:
	@if [ ! -d tmp_$(MACH) ] ; then mkdir tmp_$(MACH) ; fi
	@if [ ! -d bin ] ; then mkdir bin ; fi

clean:
	@/bin/rm -fR tmp_*

allclean: clean
	/bin/rm -f bin/* SIGMA/*~ BSE/*~ TDLDA/*~ W_BLIP/*~ \
	shared/*~ utils/*~ config/*~

cleanall: allclean

