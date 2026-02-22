reset


#set autoscale xfixmin
#set autoscale xfixmax
set xrange [0:500]
#set ytics nomirror ## modificacion ##
set yrange  [-0.1:1.1]
#set y2tics nomirror ## modificacion ##
#set y2range [-0.1:1.1] ## modificacion ##

set tics font ",15"
set xlabel "Dominio (x)" font ",15"
set ylabel "Densidad poblacional de oveja" font ",15"
#set y2label "Difusión efectiva" font ",15" ## modificacion ##

# con columna de x
	#plot fname u 1:2 title "ovejas" w l lw 3 lc rgb #	"blue", #fname u 1:3 title "pumas" w l lw 3 lc rgb 	#"red"


# sin columna de x:
	#plot fname u 1 index 1 title "ovejas" w l lw 3 lc rgb "blue" #, fname u 2 title "pumas" w l lw 3 lc rgb "red"

## animation

set terminal gif animate
set output "todas_rho.gif"
set samples 5000
do for[i = 0 : 99] {
plot "test_rho_01_zona2.dat" u 1 index i title "rho = 0.1" w l axes x1y1 lw 3 lc rgb "black",\
 "test_rho_05_zona2.dat" u 1 index i title "rho = 0.5" w l axes x1y1 lw 3 lc rgb "blue",\
"test_rho_097_zona2.dat" u 1 index i title "rho = 0.97" w l axes x1y1 lw 3 lc rgb "red" 

pause 0.05
}

