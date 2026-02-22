reset
fname = "test_con_D.dat"


set xrange [0:500]
set ytics nomirror ## modificacion ##
set yrange  [-0.1:1.1]
set y2tics nomirror ## modificacion ##
set y2range [-0.1:1.1] ## modificacion ##

set tics font ",15"
set xlabel "Dominio (x)" font ",15"
set ylabel "Densidad poblacional de oveja" font ",15"
set y2label "Difusión efectiva" font ",15" ## modificacion ##


## animation

set samples 5000
system("rm ./plots_cond_D/test_con_D_*.png")
set terminal pngcairo

do for[i = 0 : 49] {
set output sprintf('./plots_cond_D/test_con_D_%03d.png', i)
plot fname u 1 index i title "rho = 0.5" w l axes x1y1 lw 3 lc rgb "blue", fname u 2 index i title "D_eff" w l axes x1y2 lw 3 lc rgb "green"

}

system("convert -delay 25 -loop 0 ./plots_cond_D/test_con_D_*.png ./plots_cond_D/test_con_D.gif")

#system("rm test_con_D_*.png")






