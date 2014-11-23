ps = 0; eps = 0; phys_app = 0;

if (!ps) set terminal x11;
if (ps && !eps) \
  set terminal postscript enhanced color solid lw 2 "Times-Roman" 20;
if (ps && eps) \
  set terminal postscript eps enhanced color solid lw 2 "Times-Roman" 20;

set grid;

set style line 1 lt 1 lw 1 lc rgb "blue";
set style line 2 lt 1 lw 1 lc rgb "green";
set style line 3 lt 1 lw 1 lc rgb "red";

# draw projection of mechanical aperture
Ax = 17.5; Ay = 12.5;
beta_max_y = 25.5; beta_inj_y =  3.1;

if (phys_app) \
  x_hat = Ax; y_hat = Ay*sqrt(beta_inj_y/beta_max_y); \
  set arrow from -x_hat, 0.0 to -x_hat, y_hat nohead \
  lt 1 lw 1 lc rgb "black"; \
  set arrow from -x_hat, y_hat to  x_hat, y_hat nohead \
  lt 1 lw 1 lc rgb "black"; \
  set arrow from  x_hat, y_hat to  x_hat,   0.0 nohead \
  lt 1 lw 1 lc rgb "black";

if (ps) set output "dynap_err_1.ps"
set title "Dynamic Aperture\n";
set xlabel "x [mm]"; set ylabel "y [mm]";
plot "DA_bare_0.00.out" using  1:2 title "bare" with linespoints ls 1, \
     "DA_real_0.00.out" using  1:2 notitle with points ls 3;
if (!ps) pause -1;

unset arrow;

if (ps && !eps) \
  set terminal postscript portrait enhanced color solid lw 2 "Times-Roman" 20;

if (ps) set output "dynap_err_2.ps"

set multiplot;

set size 1.0, 0.5; set origin 0.0, 0.5;
set title "Horizontal Momentum Aperture\n";
set xlabel "{/Symbol d} [%]"; set ylabel "x^ [mm]";
set yrange [0:];
plot "DA_bare.out" using 1:5 title "bare" with linespoints ls 2, \
     "DA_real.out" using 1:11:13 title "w errors" with errorbars ls 1, \
     "DA_real.out" using 1:11 notitle with lines ls 1;

set origin 0.0, 0.0;
set title "Vertical Momentum Aperture\n";
set xlabel "{/Symbol d} [%]"; set ylabel "y^ [mm]";
set yrange [0:];
plot "DA_bare.out" using 1:6 title "bare" with linespoints ls 2, \
     "DA_real.out" using 1:14:16 title "w errors" with errorbars ls 3, \
     "DA_real.out" using 1:14 notitle with lines ls 3;

unset multiplot;
if (!ps) pause -1;

if (ps) set output "dynap_err_3.ps"

set multiplot;

set size 1.0, 0.5; set origin 0.0, 0.5;
set title "Horizontal Momentum Acceptance\n";
set xlabel "{/Symbol d} [%]"; set ylabel "A_x [mm{/Symbol \327}mrad]";
set yrange [0:];
plot "DA_bare.out" using 1:3 title "bare" with linespoints ls 2, \
     "DA_real.out" using 1:5:7 title "w errors" with errorbars ls 1, \
     "DA_real.out" using 1:5 notitle with lines ls 1;

set origin 0.0, 0.0;
set title "Vertical Momentum Acceptance\n";
set xlabel "{/Symbol d} [%]"; set ylabel "A_y [mm{/Symbol \327}mrad]";
set yrange [0:];
plot "DA_bare.out" using 1:4 title "bare" with linespoints ls 2, \
     "DA_real.out" using 1:8:10 title "w errors" with errorbars ls 3, \
     "DA_real.out" using 1:8 notitle with lines ls 3;

unset multiplot;
if (!ps) pause -1;
