# gplot.plt
set title "Question 1: Cholesky Decomp"
set nokey
set grid
set xlabel "x"
set ylabel "y"
m="plot.dat"
n="regression.dat"
set terminal png size 1600, 900
set output "plot.png"
plot m u 1:2, n u 1:2 w l
