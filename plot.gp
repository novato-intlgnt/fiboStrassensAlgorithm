# Configuración general
set terminal pngcairo size 1024,768 font "Arial,12"
set grid
set key top left
set xlabel "Tamaño de la matriz (n x n)"
set ylabel "Tiempo (microsegundos)"

# Usar escala logarítmica es FUNDAMENTAL para ver la diferencia en algoritmos O(n^3)
set logscale xy

set logscale xy
set output './proofs/graph_loglog.png'
set title "Comparación de Algoritmos (Escala Log-Log)"

plot "results.dat" using 1:2 with linespoints lw 2 title "Ord-mult", \
     "results.dat" using 1:3 with linespoints lw 2 title "DC-mult", \
     "results.dat" using 1:4 with linespoints lw 2 title "Str-mult"


unset logscale             
set output './proofs/graph_normal.png'
set title "Comparación de Algoritmos (Escala Lineal)"

# Reajustar rangos si es necesario, o dejar que gnuplot lo haga solo
plot "results.dat" using 1:2 with linespoints lw 2 title "Ord-mult", \
     "results.dat" using 1:3 with linespoints lw 2 title "DC-mult", \
     "results.dat" using 1:4 with linespoints lw 2 title "Str-mult"

set output
