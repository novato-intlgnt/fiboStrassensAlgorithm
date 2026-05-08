set terminal pngcairo size 800,600
set output './proofs/intersection_zoom.png'

set title "Punto de Intersección: Iterativo vs Strassen"
set xlabel "n"
set ylabel "Tiempo (microsegundos)"

# Quitamos el logscale para ver la magnitud real del cruce
unset logscale

# Ajustamos el rango a donde ocurre la competencia (ej. de 1024 a 4096)
set xrange [1024:4096]
set yrange [500000:60000000] # Ajusta según tus datos de avg

plot "results.dat" using 1:2 with linespoints lw 2 title "Iterativo", \
     "results.dat" using 1:6 with linespoints lw 2 title "Strassen"
