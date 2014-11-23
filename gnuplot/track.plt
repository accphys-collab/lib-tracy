ps = 0; eps = 0; J_phi = 0;

if (!ps) set terminal x11;
if (ps && !eps) set terminal postscript enhanced color solid;
if (ps && eps) set terminal postscript eps enhanced color solid;

set grid;

set style line 1 lt 1 lw 1 lc rgb "blue";
set style line 2 lt 1 lw 1 lc rgb "green";
set style line 3 lt 1 lw 1 lc rgb "red";

if (ps) set output "track_1.ps"
if (!J_phi) \
  set title "Horizontal Phase Space"; \
  set xlabel "x [mm]"; set ylabel "p_x [mrad]"; \
  plot "track.out" using 2:3 notitle with points ls 1;
if (J_phi) \
  set title "Horizontal Action-Angle Variables"; \
  set xlabel "{/Symbol f}_x [rad]"; set ylabel "J_x [mm.mrad]"; \
  set yrange [0:]; \
  plot "dJ.out" using 3:2 notitle with points ls 1;
if (!ps) pause -1;

if (ps) set output "track_2.ps"
if (!J_phi) \
  set title "Vertical Phase Space"; \
  set xlabel "y [mm]"; set ylabel "p_y [mrad]"; \
  plot "track.out" using 4:5 notitle with points ls 3;
if (J_phi) \
  set title "Vertical Action-Angle Variables"; \
  set xlabel "{/Symbol f}_y [rad]"; set ylabel "J_y [mm.mrad]"; \
  set yrange [0:]; \
  plot "dJ.out" using 5:4 notitle with points ls 3;
if (!ps) pause -1;

if (ps) set output "track_3.ps"
set title "Longitudinal Phase Space";
set xlabel "{/Symbol f} [{/Symbol \260}]";
set ylabel "{/Symbol d} [%]";
set yrange [*:*];
plot "track.out" using 7:6 notitle with points ls 2;
if (!ps) pause -1;
