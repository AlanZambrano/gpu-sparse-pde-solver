#!/usr/bin/gnuplot 

reset
set terminal pngcairo
set output "fig_time_vs_N.png"

fname = "compute-6-10.dat"

set xlabel "size (N)"

set logscale x
set logscale y

set key top left

plot fname using 1:3 title "time (ms)" with linespoints linewidth 3 linecolor rgb "red", \
fname using 1:2 title "used memory (MiB)" with linespoints linewidth 3 linecolor rgb "blue"
