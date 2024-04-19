#!/usr/bin/perl
#
# Builds a scissors operator from the width of occupied levels and
# the HOMO-LUMO gap. Input is a file with the same format as hmat_qp.
#
# usage: ./make_scissors.pl < [in_file ] > [out_file]
#
# Copyright (C) 2009 Murilo Tiago, Univ. of Texas, Austin, TX, USA
# mtiago@ices.utexas.edu
#
# First version written by Murilo Tiago, Oak Ridge National Laboratory, March 2008.
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

# System-dependent parameters: 
$up_B = 1;               # lowest occupied level, spin 1
$up_HOMO =  3;           # highest occupied level, spin 1
$up_LUMO =  4;           # lowest unoccupied level, spin 1
$up_MAX = 1000;          # maximum number of orbitals, spin 1
$down_B = 1;             # lowest occupied level, spin 2
$down_HOMO =  1;         # highest occupied level, spin 2
$down_LUMO =  2;         # lowest unoccupied level, spin 2
$down_MAX = 1000;        # maximum number of orbitals, spin 2

############################################################

$input = shift ;
$nspin = 1;

while (<>) {
    ($j1,$j2,$eqp,$jsp,$elda) = split(' ');
#
# Search for lowest occupied level, spin 1
#
    if ($j1 == $up_B && $j2 == $up_B && $jsp == 1) {
       $up_delta = $eqp - $elda;
       $up_bottom = $elda;
#print "match B $up_delta $up_bottom \n";
    }
#
# Search for highest occupied level, spin 1
#
    if ($j1 == $up_HOMO && $j2 == $up_HOMO && $jsp == 1) {
       $up_0= $elda;
       $up_const = $eqp - $elda;
       $up_top = $elda;
#print "match VBM $up_const $up_top \n";
    }
#
# Search for lowest unoccupied level, spin 1
#
    if ($j1 == $up_LUMO && $j2 == $up_LUMO && $jsp == 1) {
       $up_c= $elda;
       $up_const_c = $eqp - $elda;
#print "match CBM $up_c $up_const_c \n";
    }
#
# Search for lowest occupied level, spin 2
#
    if ($j1 == $down_B && $j2 == $down_B && $jsp == 2) {
       $nspin = 2;
       $down_delta = $eqp - $elda;
       $down_bottom = $elda;
#print "match B $down_delta $down_bottom \n";
    }
#
# Search for highest occupied level, spin 2
#
    if ($j1 == $down_HOMO && $j2 == $down_HOMO && $jsp == 2) {
       $nspin = 2;
       $down_0= $elda;
       $down_const = $eqp - $elda;
       $down_top = $elda;
#print "match VBM $down_const $down_top \n";
    }
#
# Search for lowest unoccupied level, spin 2
#
    if ($j1 == $down_LUMO && $j2 == $down_LUMO && $jsp == 2) {
       $nspin = 2;
       $down_c= $elda;
       $down_const_c = $eqp - $elda;
#print "match CBM $down_c $down_const_c \n";
    }
}

#
# Calculate scissors slopes and print them out.
#
if ($up_bottom == $up_top) {
    $up_slope = 0
}
else {
    $up_slope = ($up_delta - $up_const)/ ($up_bottom - $up_top);
}

if ($down_bottom == $down_top) {
    $down_slope = 0
} else {
    $down_slope = ($down_delta - $down_const)/ ($down_bottom - $down_top);
}

#
# Print parameters out.
#
print "begin scissors \n";
print "1 $up_B $up_HOMO $up_const $up_0 $up_slope \n";
print "1 $up_LUMO $up_MAX $up_const_c $up_c 0.0 \n";
if ( $nspin == 2) {
    print "2 $down_B $down_HOMO $down_const $down_0 $down_slope \n";
    print "2 $down_LUMO $down_MAX $down_const_c $down_c 0.0 \n";
}
print "end scissors \n";
