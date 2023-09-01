set title 'Mach number vs Time'
set xlabel 'Time'
set ylabel 'Mach'

set format x "%g"
set format y "%1.4f"

set grid

plot 'gnuplot/log.dat' notitle with lines

#pause -1 "Hit any key to continue"

#refresh
#replot
while (1) {
    replot
    pause 1
}
