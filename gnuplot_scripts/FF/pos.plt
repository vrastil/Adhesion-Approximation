
cd 'output/run'
set terminal pngcairo size 1280,960

A = 187826
n = 1
z_in = 200
z_out = 10

Omega_0 = 1
h = 0.67
q(x) = x / (Omega_0*h)

P(x) = A*(x**n)
T(x) = log(1+2.34*q(x))/(2.34*q(x))*(1 + 3.89*q(x) + (16.1*q(x))**2 + (5.4*q(x))**3 + (6.71*q(x))**4)**(-1./4.);

P_0(x) = P(x)*T(x)*T(x)
P_i(x) = P_0(x)/(z_in+1)**2
P_f(x) = P_0(x)/(z_out+1)**2

f_0(x)= 0
f_1(x)= -0.2

set xlabel "x [Mpc/h]"
set ylabel "y [Mpc/h]"

set output 'cut_FF_z200.png'
plot 'par_cut_FF_z200.dat' u 1:2 pt 7 ps 0.4 t 'FF, z = 200', 'track_par_pos_FF_z200.dat' pt 7 ps 0.6 lc rgb "red"
set output 'cut_FF_z65.png'
plot 'par_cut_FF_z65.dat' u 1:2 pt 7 ps 0.4 t 'FF, z = 65', 'track_par_pos_FF_z65.dat' w l lt 1 lw 4 lc rgb "red"
set output 'cut_FF_z39.png'
plot 'par_cut_FF_z39.dat' u 1:2 pt 7 ps 0.4 t 'FF, z = 39', 'track_par_pos_FF_z39.dat' w l lt 1 lw 4 lc rgb "red"
set output 'cut_FF_z27.png'
plot 'par_cut_FF_z27.dat' u 1:2 pt 7 ps 0.4 t 'FF, z = 27', 'track_par_pos_FF_z27.dat' w l lt 1 lw 4 lc rgb "red"
set output 'cut_FF_z21.png'
plot 'par_cut_FF_z21.dat' u 1:2 pt 7 ps 0.4 t 'FF, z = 21', 'track_par_pos_FF_z21.dat' w l lt 1 lw 4 lc rgb "red"
set output 'cut_FF_z17.png'
plot 'par_cut_FF_z17.dat' u 1:2 pt 7 ps 0.4 t 'FF, z = 17', 'track_par_pos_FF_z17.dat' w l lt 1 lw 4 lc rgb "red"
set output 'cut_FF_z14.png'
plot 'par_cut_FF_z14.dat' u 1:2 pt 7 ps 0.4 t 'FF, z = 14', 'track_par_pos_FF_z14.dat' w l lt 1 lw 4 lc rgb "red"
set output 'cut_FF_z12.png'
plot 'par_cut_FF_z12.dat' u 1:2 pt 7 ps 0.4 t 'FF, z = 12', 'track_par_pos_FF_z12.dat' w l lt 1 lw 4 lc rgb "red"
set output 'cut_FF_z10.png'
plot 'par_cut_FF_z10.dat' u 1:2 pt 7 ps 0.4 t 'FF, z = 10', 'track_par_pos_FF_z10.dat' w l lt 1 lw 4 lc rgb "red"