# gplot.plt
set title "Question 1: Householder QR Decomp"
set nokey
set grid
set xlabel "x"
set ylabel "y"
m="plot.dat"
n="qr.dat"
set terminal png size 1600, 900
set output "qr.png"
plot m u 1:2, n u 1:2 w l
