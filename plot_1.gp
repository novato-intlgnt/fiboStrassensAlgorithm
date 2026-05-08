# Configuración general
set terminal pngcairo size 1024,768 font "Arial,12"
set grid
set key top left
set xlabel "Tamaño de la matriz (n x n)"
set ylabel "Tiempo (microseg)"

# --- Gráfico 1: Líneas simples ---
set output './proofs/graph-nlg.png'
set title "Comparación de Algoritmos"

plot "results_range.dat" using 1:2 with linespoints lw 2 title "Ord-mult", \
     "results_range.dat" using 1:3 with linespoints lw 2 title "DC-mult", \
     "results_range.dat" using 1:4 with linespoints lw 2 title "Str-mult"

