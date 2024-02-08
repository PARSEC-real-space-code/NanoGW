#!/bin/bash
#
# Script that performs a series of self-consistent GW iterations using
# a scissors operator.
#
# input files in work directory :
#    parsec_dft.dat  : DFT wave-functions and eigenvalues
#    rgwbs.in_0      : parameter file in zero-th iteration
#    scissors.pl     : appropriate script to build scissors operator
#    occup.in        : orbital occupancies, optional
#
# Output of n-th iteration is kept under subdirectory iteration_n.
#
# Copyright (C) 2009 Murilo Tiago, Univ. of Texas, Austin, TX, USA
# mtiago@ices.utexas.edu
#
# First version written by Murilo Tiago, Oak Ridge National Laboratory, February 2008.
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 1, or (at your option)
#    any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston MA  02110-1301 USA
#

# Bindings for the executables:
sigma_run='~/bin/sigma.ser'

# number of iterations
niter=5

# starting iteration
nstart=0

###############################################################
#  No modifications needed below this point.
############################################################### 
echo 'Start calculation of ' $niter ' iterations starting from ' $nstart

#
# start the queue of jobs
#
iter=$nstart

while (( iter<=$niter ))
do

    let "iterm1 = iter - 1"

    echo
    echo ------------------------------------------------------------------

# delete any junk directory and reinitialize work directory
    if [[ -d iteration_$iter ]]
    then
        rm -fR iteration_$iter
    fi
    echo running iteration $iter
    mkdir iteration_$iter
    cd iteration_$iter
    if [[ -f ../occup.in ]]
    then
        ln -s ../occup.in .
    fi
    if [[ -f ../parsec.in ]]
    then
        ln -s ../parsec.in .
    fi

# input files: parsec_dft.dat, rgwbs.in_0
# output files: parsec_qp.dat, sigma_mtxel_qp.dat
        ln -s ../parsec_dft.dat parsec.dat
        cp -p ../rgwbs.in_0 rgwbs.in
        cat rgwbs.in ../scissors > tt
        rm -f rgwbs.in
        mv tt rgwbs.in

        date
        echo ------------------------------------------------------------------

        $sigma_run | tee sigma.out

        if [[ ! -f hmat_qp ]]
        then
           exit
        fi

        ../scissors.pl < hmat_qp > ../scissors

        echo
        echo ------------------------------------------------------------------
        echo

    cd ../

    let "iter = iter + 1"

    if [[ -f stop_scgw ]]
    then
       exit
    fi

done
