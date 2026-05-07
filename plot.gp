set terminal pngcairo size 800,600
set output './proofs/graph.png'

set title "Análisis Empírico: Multiplicacion de Matrices"
set xlabel "Tamaño de entrada (n)"
set ylabel "Tiempo (microsegundos)"
set key top left
set grid

plot "results.dat" using 1:2 with linespoints title "ordSort", \
     "results.dat" using 1:3 with linespoints title "DC-Sort", \
     "results.dat" using 1:4 with linespoints title "ST-Sort"

