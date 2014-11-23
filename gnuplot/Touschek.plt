ps = 0; eps = 0;

if (!ps) set terminal x11;
if (ps && !eps) set terminal postscript enhanced color solid;
if (ps && eps) set terminal postscript eps enhanced color solid;

set grid;

if (ps) set output "Touschek.ps"
set title "Momentum Aperture";
set xlabel "s [m]";
set ylabel "{/Symbol d} [%]";
plot "mom_aper.out" using 2:3 notitle with fsteps ls 2, \
     "mom_aper.out" using 2:4 notitle with fsteps ls 2;
if (!ps) pause -1;
