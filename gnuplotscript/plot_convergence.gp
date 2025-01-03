set terminal png size 800,600
set output "../assets/convergence_history.png"

set title "Convergence History"
set xlabel "Iteration"
set ylabel "Residual Norm"
set ytics 0.2
set yrange [0:1]
set grid
set style line 1 linecolor rgb "#0060ad" linewidth 2

plot "../data/convergence_history.dat" using 1:2 with lines linestyle 1 title "Residual Norm"