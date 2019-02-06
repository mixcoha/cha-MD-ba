#!/bin/bash

cd $1/R1/prod/cell
RUNDIR=$PWD

mv tprout.tpr topol.tpr

xtc_file=$RUNDIR/cell.xtc
tpr_file=$RUNDIR/topol.tpr

echo -e "q \n" | gmx make_ndx -f topol.tpr

ndx_file=$RUNDIR/index.ndx

mkdir analysis
cd analysis

## ******   RMSD CALCULATION   **** ###

mkdir rms
cd rms

echo 4 4 | gmx rms -f $xtc_file -s $tpr_file -dist
echo 4 13 | gmx rms -f $xtc_file -s $tpr_file -dist -o rmsd_liginprot.xvg 
echo 13 13 | gmx rms -f $xtc_file -s $tpr_file -dist -o rmsd_lig.xvg 

### DESELECCIONAR ESTAS OPCIONES PARA EL CÁLCULO DE LA MATRIZ DE RMSD
#echo 4 4 | gmx rms -f $xtc_file -s $tpr_file -dist -m
#gmx xpm2ps -f rmsd.xpm

sed -e '1,17d' rmsd.xvg | awk '{ print $1/1000, $2 }' > rmsd.txt
sed -e '1,17d' rmsd_lig.xvg | awk '{ print $1/1000, $2 }' > rmsd-lig.txt
sed -e '1,17d' rmsd_liginprot.xvg | awk '{ print $1/1000, $2 }' > rmsd-liginprot.txt

##### GRAFICA RMSD PROTEINA ####
gnuplot << EOF
set terminal postscript eps color enhanced
set output "rmsd.eps"

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

#set tmarg 1
#set bmarg 1
#set lmarg 10
#set rmarg 3
set title "$1 RMSD"
set format x"%g"   # reset xformat
set format y"%g"
set xlabel "Time (ns)" font "Helvetica,14" # set xlabel
set ylabel "RMSD (nm)" font "Helvetica,14"
set xtics 50
set mxtics 5
set ytics 0.1 nomirror
#set autoscale xy  #
set xrange [GPVAL_X_MIN:GPVAL_X_MAX]
set yrange [GPVAL_Y_MIN:GPVAL_Y_MAX]
set key ins
set key horiz
set key default
set key outside left bottom
set key noinvert samplen 0.2 spacing 2 width 3 height 0
set key outside horizontal

f(x) = mean_y
fit f(x) 'rmsd.txt' u 1:2 via mean_y

stddev_y = sqrt(FIT_WSSR / (FIT_NDF + 1 ))

# Plotting the range of standard deviation with a shaded background
set label 1 gprintf("~m{.4-} = %g", mean_y) at 200, mean_y-stddev_y-0.1
set label 2 gprintf("{/Symbol s} = %g", stddev_y) at 200, mean_y-stddev_y-0.2
plot mean_y-stddev_y with filledcurves y1=mean_y lt 1 lc rgb "#EAEAF4" title "~m{.4-} +/- {/Symbol s}" , \
mean_y+stddev_y with filledcurves y1=mean_y lt 1 lc rgb "#EAEAF4" notitle, \
mean_y w l ls 3 title "~m{.4-}", 'rmsd.txt' u 1:2 w l ls 1 t 'RMSD'

EOF


##### GRAFICA RMSD LIGANDO ####
gnuplot << EOF
set terminal postscript eps color enhanced
set output "rmsd-lig.eps"

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

#set tmarg 1
#set bmarg 1
#set lmarg 10
#set rmarg 3
set title "$1 RMSD"
set format x"%g"   # reset xformat
set format y"%g"
set xlabel "Time (ns)" font "Helvetica,14" # set xlabel
set ylabel "RMSD (nm)" font "Helvetica,14"
set xtics 50
set mxtics 5
set ytics 0.1 nomirror
#set autoscale xy  #
set xrange [GPVAL_X_MIN:GPVAL_X_MAX]
set yrange [GPVAL_Y_MIN:GPVAL_Y_MAX]
set key ins
set key horiz
set key default
set key outside left bottom
set key noinvert samplen 0.2 spacing 2 width 3 height 0
set key outside horizontal

f(x) = mean_y
fit f(x) 'rmsd-lig.txt' u 1:2 via mean_y

stddev_y = sqrt(FIT_WSSR / (FIT_NDF + 1 ))

# Plotting the range of standard deviation with a shaded background
set label 1 gprintf("~m{.4-} = %g", mean_y) at 200, mean_y-stddev_y-0.1
set label 2 gprintf("{/Symbol s} = %g", stddev_y) at 200, mean_y-stddev_y-0.2
plot mean_y-stddev_y with filledcurves y1=mean_y lt 1 lc rgb "#EAEAF4" title "~m{.4-} +/- {/Symbol s}" , \
mean_y+stddev_y with filledcurves y1=mean_y lt 1 lc rgb "#EAEAF4" notitle, \
mean_y w l ls 3 title "~m{.4-}", 'rmsd-lig.txt' u 1:2 w l ls 1 t 'RMSD', 'rmsd-liginprot.txt' u 1:2 w l ls 2 t 'RMSD Lig in MOR1'

EOF


cd $RUNDIR/analysis

## ******   RoG CALCULATION   **** ###

mkdir rog
cd rog

echo 1 | gmx gyrate -f $xtc_file -s $tpr_file

sed -e '1,17d' gyrate.xvg | awk '{ print $1/1000, $2 }' > gyrate.txt

gnuplot << EOF
set terminal postscript eps color enhanced
set output "gyrate.eps"

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


#set tmarg 1
#set bmarg 1
#set lmarg 10
#set rmarg 3
set title "$1 Radius of gyration"
set format x"%g"   # reset xformat
set format y"%g"
set xlabel "Time (ns)" font "Helvetica,14" # set xlabel
set ylabel "RoG (nm)" font "Helvetica,14"
set xtics 50
set mxtics 5
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
set label 1 gprintf("~m{.4-} = %g", mean_y) at 400, mean_y-stddev_y-0.01
set label 2 gprintf("{/Symbol s} = %g", stddev_y) at 400, mean_y-stddev_y-0.02
plot mean_y-stddev_y with filledcurves y1=mean_y lt 1 lc rgb "#EAEAF4" title "~m{.4-} + {/Symbol s}" , \
mean_y+stddev_y with filledcurves y1=mean_y lt 1 lc rgb "#EAEAF4" title "~m{.4-} - {/Symbol s}", \
mean_y w l ls 3 title "~m{.4-}", 'gyrate.txt' u 1:2 w l ls 1 title "RoG"

EOF

cd $RUNDIR/analysis

## ***** RMSF ***** ##

mkdir rmsf
cd rmsf

echo 3 | gmx rmsf -f $xtc_file -s $tpr_file -n $ndx_file -oq -ox -od

ed -n '17,$p' rmsf.xvg | awk '{print $2}'  > p1
sed -n 's/ATOM/&/p' xaver.pdb | awk '{print $6}'  > p2

paste p2 p1 > rmsf.txt

gnuplot << EOF

##########################################
set terminal postscript eps color enhanced
set output "rmsf.eps"
##########################################

set style data histogram
set style histogram clustered gap 1

set ylabel "RMSF" font ",18" offset character 2.5,0,0

set key under nobox
set style histogram clustered gap 1 title offset 2,0.25
set style fill solid noborder
set boxwidth 0.95
unset xtics
set xtics nomirror rotate by -45 scale 0 font ",8"
set ytics font ",8"
set grid y
set auto y
set xrange [0:400]

plot newhistogram, 'rmsf.txt' using 2:xtic(50) t 'RMSF' lc rgbcolor "blue" lt 1

EOF

cd $RUNDIR/analysis

## ******   Hydrogen Bonds CALCULATION   **** ###

mkdir hbonds
cd hbonds

echo 1 1 | gmx hbond -f $xtc_file -s $tpr_file -n $RUNDIR/index.ndx -num prot-prot.xvg -tu ns
echo 1 13 | gmx hbond -f $xtc_file -s $tpr_file -n $RUNDIR/index.ndx -num prot-lig.xvg -tu ns

sed -e '1,24d' prot-prot.xvg | awk '{ print $1, $2 }' > hbond-prot.txt
sed -e '1,24d' prot-lig.xvg | awk '{ print $1, $2 }' > hbond-PL.txt



gnuplot << EOF

set terminal postscript eps color enhanced
set output "hbonds.eps"


set style line 1 lt -1 lw 2 lc rgb "black"
set style line 2 lt -1 lw 2 lc rgb "red"
set style line 3 lt -1 lw 2 lc rgb "gold"
set style line 4 lt -1 lw 2 lc rgb "blue"
set style line 5 lt -1 lw 2 lc rgb "green"
set style line 6 lt -1 lw 2 lc rgb "cyan"
set style line 7 lt -1 lw 2 lc rgb "orange

set xlabel "Time [ns]"
set ylabel "Total hydrogen bond"
set title "Total hydrogen bonds intra Protein"
set autoscale ymax
set autoscale xmax
#set xrange [ 0 : 250 ]
#set yrange [ 0 : 60 ]
set key default
set key samplen .2
set key center top
set key ins vert
set key horiz
set format x "%g"

f(x) = mean_y
fit f(x) 'hbond-prot.txt' u 1:2 via mean_y

stddev_y = sqrt(FIT_WSSR / (FIT_NDF + 1 ))

# Plotting the range of standard deviation with a shaded background
set label 1 gprintf("~m{.4-} = %g", mean_y) at 400, mean_y-stddev_y-0.01
set label 2 gprintf("{/Symbol s} = %g", stddev_y) at 400, mean_y-stddev_y-0.02
plot mean_y-stddev_y with filledcurves y1=mean_y lt 1 lc rgb "#EAEAF4" title "~m{.4-} + {/Symbol s}" , \
mean_y+stddev_y with filledcurves y1=mean_y lt 1 lc rgb "#EAEAF4" title "~m{.4-} - {/Symbol s}", \
mean_y w l ls 3 title "~m{.4-}", 'hbond-prot.txt' u 1:2 w l ls 1 title "Prot-Prot"
EOF

##### Gráfica del los puentes de hidrógeno de la prot con el lig

gnuplot << EOF

set terminal postscript eps color enhanced
set output "hbonds_PL.eps"


set style line 1 lt -1 lw 2 lc rgb "black"
set style line 2 lt -1 lw 2 lc rgb "red"
set style line 3 lt -1 lw 2 lc rgb "gold"
set style line 4 lt -1 lw 2 lc rgb "blue"
set style line 5 lt -1 lw 2 lc rgb "green"
set style line 6 lt -1 lw 2 lc rgb "cyan"
set style line 7 lt -1 lw 2 lc rgb "orange

set xlabel "Time [ns]"
set ylabel "Total hydrogen bond"
set title "Total hydrogen bonds Protein-Ligand"
set autoscale ymax
set autoscale xmax
#set xrange [ 0 : 250 ]
#set yrange [ 0 : 60 ]
set key default
set key samplen .2
set key center top
set key ins vert
set key horiz
set format x "%g"

f(x) = mean_y
fit f(x) 'hbond-PL.txt' u 1:2 via mean_y

stddev_y = sqrt(FIT_WSSR / (FIT_NDF + 1 ))

# Plotting the range of standard deviation with a shaded background
set label 1 gprintf("~m{.4-} = %g", mean_y) at 400, mean_y-stddev_y-0.01
set label 2 gprintf("{/Symbol s} = %g", stddev_y) at 400, mean_y-stddev_y-0.02
plot mean_y-stddev_y with filledcurves y1=mean_y lt 1 lc rgb "#EAEAF4" title "~m{.4-} + {/Symbol s}" , \
mean_y+stddev_y with filledcurves y1=mean_y lt 1 lc rgb "#EAEAF4" title "~m{.4-} - {/Symbol s}", \
mean_y w l ls 3 title "~m{.4-}", 'hbond-PL.txt' u 1:2 w l ls 1 title "Prot-Lig"
EOF

cd $RUNDIR/analysis

##### ***** Distance matrices consisting of the smallest distance **** #####

mkdir mdmat
cd mdmat

echo 1 | gmx mdmat -f $xtc_file -s $tpr_file -n $RUNDIR/index.ndx -mean -no
gmx xpm2ps -f dm.xpm

cd ..

##### ***** Clusters  **** #####
mkdir cluster
cd cluster

echo 4 2 | gmx cluster -f $xtc_file -s $tpr_file -n $RUNDIR/index.ndx -tu ns -o -g -sz -dist -cl -method gromos -cutoff 0.35
##### DESELECCIONAR ESTA OPCIÓN Y COMENTAR LA ANTERIOR CUANDO ESTÉ LA MATRIZ DE RMSD CALCULADA
#echo 4 2 | gmx cluster -f $xtc_file -s $tpr_file -n $RUNDIR/index.ndx -dm $RUNDIR/analysis/rms/rmsd.xpm -tu ns -o -g -sz -dist -cl -method gromos -cutoff 0.35

gmx xpm2ps -f rmsd-clust.xpm

sed -n '18,$p' clust-size.xvg > clust-size.txt

gnuplot << EOF
set terminal postscript eps color enhanced
set output "cluster.eps"

set border 3 front linetype -1 linewidth 1.000
set boxwidth 0.95 absolute
set style fill solid 1.00 noborder
set grid nopolar
set grid noxtics nomxtics ytics nomytics noztics nomztics \
nox2tics nomx2tics noy2tics nomy2tics nocbtics nomcbtics
set grid layerdefault   linetype 0 linewidth 1.000,  linetype 0 linewidth 1.000
#set key bmargin center horizontal Left reverse noenhanced autotitles columnhead nobox
set key default
set key samplen .5
set key center top
set key ins vert
set key horiz
set key right top
set style histogram clustered gap 1.5 title  offset character 2, 0.25, 0
set datafile missing '-'
set style data histograms
set xtics border in scale 0,0 nomirror rotate by -90  offset character 0, 0, 0 autojustify
set xtics  norangelimit font ",8"
set xtics   ()
set ytics border in scale 0,0 mirror norotate  offset character 0, 0, 0 autojustify
set ztics border in scale 0,0 nomirror norotate  offset character 0, 0, 0 autojustify
set cbtics border in scale 0,0 mirror norotate  offset character 0, 0, 0 autojustify
set rtics axis in scale 0,0 nomirror norotate  offset character 0, 0, 0 autojustify
#set title "RMSF"
set xlabel "Cluster ID"
set xlabel  offset character 0, -2, 0 font "" textcolor lt -1 norotate
set ylabel "Structures"
set autoscale ymax
set autoscale xmax

plot 'clust-size.txt' using 2:xtic(1) t 'Cluster'

EOF

cd $RUNDIR/analysis

##### ***** PCA  **** #####
mkdir covar
cd covar

echo 3 3 | gmx covar -f $xtc_file -s $tpr_file -n $ndx_file -o -v -av -xpm -xpma

gmx xpm2ps -f covara.xpm

echo 3 3 | gmx anaeig -v eigenvec.trr -f $xtc_file -s $tpr_file -n $RUNDIR/index.ndx -eig eigenval.xvg -proj -rmsf -3d -filt filter1-3.pdb -first 1 -last 3 -skip 100

echo 3 3 | gmx anaeig -f $xtc_file -s $tpr_file -n $RUNDIR/index.ndx -filt filter1-2.pdb -first 1 -last 2 -skip 100

cd ..
## ***** Electrostatical potential across the box ***** ##

mkdir Vbox
cd Vbox

echo 1 | gmx potential -f $xtc_file -s $tpr_file -n $ndx_file


cd $RUNDIR

