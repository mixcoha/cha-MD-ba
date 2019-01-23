#!/bin/bash

RUNDIR=$PWD
cd $RUNDIR/$1/npt/cell

cd analysis

## ******   RMSD CALCULATION   **** ###
cd $2


gnuplot << EOF
set terminal postscript eps color enhanced
set output "$2.eps"

set style line 1 lt -1 lw 3 lc rgb "red"
set style line 2 lt -1 lw 3 lc rgb "green"
set style line 3 lt -1 lw 3 lc rgb "blue"
set style line 4 lt -1 lw 3 lc rgb "magenta"
set style line 5 lt -1 lw 3 lc rgb "cyan"
set style line 6 lt -1 lw 3 lc rgb "yellow"
set style line 7 lt -1 lw 3 lc rgb "black"
set style line 8 lt -1 lw 3 lc rgb "gold"
set style line 9 lt -1 lw 3 lc rgb "orange"
set style line 10 lt -1 lw 3 lc rgb "grey"
set style line 11 lt -1 lw 3 pi -3 pt 7 ps 0.5
set style line 12 lt 2  lw 3 lc 9
set encoding iso_8859_1

#
#PLOT1 All results
#


set title "$1 $2"
set format x"%g"   # reset xformat
set format y"%g"
set xlabel "Time (ns)" font "Helvetica,14" # set xlabel
set ylabel "$2 (nm)" font "Helvetica,14"

set xrange [GPVAL_X_MIN:GPVAL_X_MAX]
set yrange [GPVAL_Y_MIN:GPVAL_Y_MAX]
set key ins
set key horiz
set key default
set key outside left bottom
set key noinvert samplen 0.2 spacing 2 width 3 height 0
set key outside horizontal



f(x) = mean_y
fit f(x) 'gyrate.txt' u 1:2 via mean_y

stddev_y = sqrt(FIT_WSSR / (FIT_NDF + 1 ))

# Plotting the range of standard deviation with a shaded background
set label 1 gprintf("~m{.4-} = %g", mean_y) at 400, mean_y-stddev_y+0.2
set label 2 gprintf("{/Symbol s} = %g", stddev_y) at 400, mean_y-stddev_y+0.3
plot mean_y-stddev_y with filledcurves y1=mean_y lt 1 lc rgb "#EAEAF4" title "~m{.4-} + {/Symbol s}" , \
mean_y+stddev_y with filledcurves y1=mean_y lt 1 lc rgb "#EAEAF4" title "~m{.4-} - {/Symbol s}", \
mean_y w l ls 3 title "~m{.4-}", 'gyrate.txt' u 1:2 w l ls 10 t '$2'

EOF

