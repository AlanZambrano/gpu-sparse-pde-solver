reset
set terminal pngcairo
set output "fig_time_mem_vs_N.png"

fname = "compute-6-10.dat"

f(x) = a*x**n
fit f(x) fname using 1:3 via a,n

#set xrange [900:100000]
#set yrange  [10000:2002000]

#set xrange [1:50000]
#set yrange  [500:2000000]

set xrange [10:200]
set yrange [100:10000000]


set logscale x
set logscale y

set tics font ",15"
set xlabel "Nodos en grilla (N)" font ",15"
#set ylabel "tiempo (ms)" font ",15"
set title "Tiempo y memoria vs tamaño del sistema" font ",18"

set key right bottom

plot fname u 1:3 title "tiempo (ms)"  w lp lw 3 lc rgb "red", \
fname u 1:2 title "used memory (MiB)" w lp lw 3 lc rgb "blue", \
f(x) title sprintf("%3.3g x^{(%3.2g)}",a,n) lw 3 lc rgb "black"
