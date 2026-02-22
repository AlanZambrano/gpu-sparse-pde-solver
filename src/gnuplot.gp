N = 10 # tamanio del cada eje
K = 20 # cantidad de fotogramas

set terminal pngcairo
fname = "oveja_puma_difusividad.dat"

set xrange[0:N]
set yrange[0:N]
set cbrange[0:1]

set tics out
set size ratio 1

set samples 5000

# Ovejas:
set title "Densidad Poblacional de Ovejas"
set xlabel "X(nodos)"
set ylabel "Y(nodos)"

do for [i = 0:K] {
	set output sprintf('../figures/ovejas_%03d.png', i)

	plot fname using 1:2:3 index i with image title "ovejas"
}

system("convert -delay 10 -loop 0 ../figures/ovejas_*.png ../figures/ovejas.gif")

system("rm ../figures/ovejas_*.png")


# Pumas:
set title "Densidad Poblacional de Pumas"
set xlabel "X(nodos)"
set ylabel "Y(nodos)"

do for [i = 0:K] {
	set output sprintf('../figures/pumas_%03d.png', i)

	plot fname using 1:2:4 index i with image title "pumas"
}

system("convert -delay 10 -loop 0 ../figures/pumas_*.png ../figures/pumas.gif")

system("rm ../figures/pumas_*.png")


# Difusividad de ovejas:
set title "Difusividad de ovejas"
set xlabel "X(nodos)"
set ylabel "Y(nodos)"

do for [i = 0:K] {
	set output sprintf('../figures/difusividad_%03d.png', i)

	plot fname using 1:2:5 index i with image title "difusividad"
}

system("convert -delay 10 -loop 0 ../figures/difusividad_*.png ../figures/difusividad.gif")

system("rm ../figures/difusividad_*.png")










