#!/bin/bash
#export DYLD_LIBRARY_PATH=/Applications/MATLAB_R2014b.app/sys/os/maci64
########################################################
cd '/Users/davidmorante/Desktop/GALILEO-DEIMOS'
n=18
tfinal=200
dir=Test3-Galileo-nocons
initialGuess=/Users/davidmorante/Desktop/MOLTO-OR-step2/Initialguess/solution_nocons.sol
########################################################
echo $n
echo $tfinal
echo $dir
echo initialGuess
########################################################
matlab -nodisplay  -nosplash -r " mkdir(['./','$dir']); fid=  fopen( ['./','$dir','/ampl_guess_coll.dat',], 'wt' ); fid2=  fopen( ['./','$dir','/model_param.dat',], 'wt' ); fid3=  fopen( ['./','$dir','/ampl_main.mod',], 'wt' ); n= $n; tfinal= $tfinal; dir='$dir';initialGuess='$initialGuess'; run ./Auxiliar/generateGuess.m; exit"

/Users/davidmorante/Desktop/AMPL/solvers/ampl "$dir/ampl_main.mod"

matlab -nodisplay  -r "dir='$dir' ; load([dir,'/output.out']);run Auxiliar/ampl_plot_results.m; exit"
#########################################################