reset

set key bottom left

set xrange [0:500]
set ytics nomirror 
set yrange  [-0.01:1.02]
set y2tics nomirror 
set y2range [-0.01:1.02] 

set tics font ",15"
set xlabel "Dominio (x)" font ",15"
set ylabel "Densidades poblacionales" font ",15"
set y2label "Difusión efectiva" font ",15" 


## animation

set samples 5000
system("rm ./seis_curvas/curvas_*.png")
set terminal pngcairo

do for [i = 0:49] {
  set output sprintf('./seis_curvas/curvas_%03d.png', i)
  plot \
    "thick_blue_red_green_e_10.dat" using 1 index i title "ovejas"           with lines axes x1y1 lw 3 lc rgb "blue" , \
    "thick_blue_red_green_e_10.dat" using 2 index i title "pumas"            with lines axes x1y1 lw 3 lc rgb "red"  , \
    "thick_blue_red_green_e_10.dat" using 3 index i title "D_eff"            with lines axes x1y2 lw 3 lc rgb "green", \
    "dotted_blue_0_green.dat"  using 1 index i title "solas" with lines axes x1y1 lw 0.5 lc rgb "blue" , \
    "dotted_0_red_1.dat"       using 2 index i title "solos" with lines axes x1y1 lw 0.5 lc rgb "red"  , \
    "dotted_blue_0_green.dat"  using 3 index i title "sola"  with lines axes x1y2 lw 0.5 lc rgb "green"
}

system("convert -delay 10 -loop 0 ./seis_curvas/curvas_*.png ./seis_curvas/seis_curvas.gif")







