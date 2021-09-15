#!/bin/csh
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#FUNCTION: make-rank-crank-modp-data - 
#                    
#
#SYNTAX: make-rank-crank-modp-data t ln smod sres
#
#PARAMETERS: 
#       t - prime
#       ln - last n
#       smod - positive integer#
#       sres - residue  (0 ... smod-1)
#
#DESCRIPTION:
# Computes  N(k, t, smod*n + sres)  mod t for k=0, .., (t-1)/2 and n=0 .. ln
# Computes  M(k, t, smod*n + sres)  mod t for k=0, .., (t-1)/2 and n=0 .. ln
#  Writes  (nr(j),j=0,t1) to file ranksave
#  Writes  (ncr(j),j=0,t1) to file cranksave
#  Writes  p(n) mod t to file ptnsave
#  Writes  t, t1, smod, sres, lj+1  to file mdim
#
#EXAMPLE:
#  make-rank-crank-modp-data 19 80000 289 277
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
set nargs = $#argv
if ($nargs == 0) then
 echo "SYNTAX: make-rank-crank-modp-data t lastn smod sres"
 echo "EXAMPLE: make-rank-crank-modp-data 19 100000 361 346"
 exit -9
endif
set t = $argv[1]
set ln = $argv[2]
set smod = $argv[3]
set sres = $argv[4]
echo "make-rank-crank-modp-data $argv"
spin 5
#############################################################################
echo "STEP 1:  Create init file inputdata2"
echo "$t" > inputdata2
echo "$ln" >> inputdata2
echo "$smod" >> inputdata2
echo "$sres" >> inputdata2
#############################################################################
echo "STEP 2:  Run fortran program makegfdat"
     spin 3
     makegfdat
     spin 3
#############################################################################
echo "STEP 3:  Convert fortran output to maple input"
     sed -f sedsp2tab ranksave > ranksave-"$t"-"$ln"-"$smod"-"$sres"
     sed -f sedsp2tab cranksave > cranksave-"$t"-"$ln"-"$smod"-"$sres"
     sed -f sedsp2tab ptnsave > ptnsave-"$t"-"$ln"-"$smod"-"$sres"
echo "SEE new files: "
     ls -dl  ranksave-"$t"-"$ln"-"$smod"-"$sres"
     ls -dl cranksave-"$t"-"$ln"-"$smod"-"$sres"
     ls -dl ptnsave-"$t"-"$ln"-"$smod"-"$sres"

