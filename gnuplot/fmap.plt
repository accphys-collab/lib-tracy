ps = 0; eps = 0;

if (!ps) set terminal x11;
if (ps && eps) set terminal postscript eps enhanced color solid;
if (ps  && !eps) set terminal postscript enhanced color solid;

N = 1; N_x = 33; N_y = 16;
#N = 3; N_x = 11; N_y = 5;
#N = 15; N_x =  2; N_y = 1;

set grid;

file_name = "`echo $TRACY_LIB`/gnuplot/jet.dat";
# Load 64-color palette for Jet
set palette model RGB file file_name using ($1/255):($2/255):($3/255);

set cbrange [-10:-2];

set noztics; unset clabel;

set cntrparam level 100;
set view map;
# To set y-axis to left side and avoid compression of color box.
unset pm3d;

if (ps) set output "fmap_1.ps"

set multiplot;

set size 1.0, 0.5; set origin 0.0, 0.5;
set title "Tune Shift";
set xlabel "{/Symbol n}_x"; set ylabel "{/Symbol n}_y";
splot "fmap.out" using \
      ((abs($3-int($3)) > 1e-6)? N*(N_x+$3) : NaN): \
      ((abs($4-int($4)) > 1e-6)? N*(N_y+$4) : NaN):7 \
      notitle w points pt 13 lt palette z;

set pm3d at b map;
#set contour;
unset colorbox;

set origin 0.0, 0.0;
set title "Diffusion Map";
set xlabel "A_x"; set ylabel "A_y";
set autoscale x; set autoscale y; 
splot "fmap.out" using 1:2:(($7 != -2.0)? $7 : NaN) notitle lt palette z;

unset multiplot;
if (!ps) pause(-1);

if (ps) set output "fmap_2.ps"

set multiplot;

unset contour;
unset pm3d;
set colorbox;

set size 1.0, 0.5; set origin 0.0, 0.5;
set title "Tune Shift";
set xlabel "{/Symbol n}_x"; set ylabel "{/Symbol n}_y";
splot "fmapdp.out" using \
      ((abs($3-int($3)) > 1e-6)? N*(N_x+$3) : NaN): \
      ((abs($4-int($4)) > 1e-6)? N*(N_y+$4) : NaN):7 \
      notitle w points pt 13 lt palette z;

set pm3d at b map;
#set contour;
unset colorbox;

set origin 0.0, 0.0;
set title "Diffusion Map";
set xlabel "{/Symbol d} [%]"; set ylabel "A_x";
set autoscale x; set autoscale y; 
splot "fmapdp.out" using 1:2:(($7 != -2.0)? $7 : NaN) notitle lt palette z;

unset multiplot;
if (!ps) pause(-1);
