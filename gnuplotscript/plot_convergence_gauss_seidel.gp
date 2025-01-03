# Gauss Seidel Convergence Graph
set terminal png size 800,600
set output "../assets/convergence_history_gauss_seidel.png"

set title "Convergence History - Gauss Seidel"
set xlabel "Iteration"
set ylabel "Residual Norm"
set ytics 0.2
set yrange [0:1]
set grid
set style line 2 linecolor rgb "#60ad00" linewidth 2

plot "../data/convergence_history_gauss_seidel.dat" using 1:2 with lines linestyle 2 title "Residual Norm (Jacobi)"