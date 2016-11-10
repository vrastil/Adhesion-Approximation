
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

set output 'cut_ZA_z200.png'
plot 'par_cut_ZA_z200.dat' u 1:2 pt 7 ps 0.4 t 'ZA, z = 200', 'track_par_pos_ZA_z200.dat' pt 7 ps 0.6 lc rgb "red"
set output 'cut_ZA_z65.png'
plot 'par_cut_ZA_z65.dat' u 1:2 pt 7 ps 0.4 t 'ZA, z = 65', 'track_par_pos_ZA_z65.dat' w l lt 1 lw 4 lc rgb "red"
set output 'cut_ZA_z39.png'
plot 'par_cut_ZA_z39.dat' u 1:2 pt 7 ps 0.4 t 'ZA, z = 39', 'track_par_pos_ZA_z39.dat' w l lt 1 lw 4 lc rgb "red"
set output 'cut_ZA_z27.png'
plot 'par_cut_ZA_z27.dat' u 1:2 pt 7 ps 0.4 t 'ZA, z = 27', 'track_par_pos_ZA_z27.dat' w l lt 1 lw 4 lc rgb "red"
set output 'cut_ZA_z21.png'
plot 'par_cut_ZA_z21.dat' u 1:2 pt 7 ps 0.4 t 'ZA, z = 21', 'track_par_pos_ZA_z21.dat' w l lt 1 lw 4 lc rgb "red"
set output 'cut_ZA_z17.png'
plot 'par_cut_ZA_z17.dat' u 1:2 pt 7 ps 0.4 t 'ZA, z = 17', 'track_par_pos_ZA_z17.dat' w l lt 1 lw 4 lc rgb "red"
set output 'cut_ZA_z14.png'
plot 'par_cut_ZA_z14.dat' u 1:2 pt 7 ps 0.4 t 'ZA, z = 14', 'track_par_pos_ZA_z14.dat' w l lt 1 lw 4 lc rgb "red"
set output 'cut_ZA_z12.png'
plot 'par_cut_ZA_z12.dat' u 1:2 pt 7 ps 0.4 t 'ZA, z = 12', 'track_par_pos_ZA_z12.dat' w l lt 1 lw 4 lc rgb "red"
set output 'cut_ZA_z10.png'
plot 'par_cut_ZA_z10.dat' u 1:2 pt 7 ps 0.4 t 'ZA, z = 10', 'track_par_pos_ZA_z10.dat' w l lt 1 lw 4 lc rgb "red"

set terminal pngcairo size 640,480

set xlabel "k [h/Mpc]"
set ylabel "P(k) [(Mpc/h)^3]"
set output 'pwr_spec_128p_128m_1024_ZA.png'
set autoscale
set logscale xy
set format y "10^{%L}"

plot 'pwr_spec_ZA_z200.dat' u 1:2 t 'ZA, z = 200' pt 1 ps 0.5, \
'pwr_spec_ZA_z39.dat' u 1:2 t 'ZA, z = 39' pt 1 ps 0.5, \
'pwr_spec_ZA_z21.dat' u 1:2 t 'ZA, z = 21' pt 1 ps 0.5, \
'pwr_spec_ZA_z14.dat' u 1:2 t 'ZA, z = 14' pt 1 ps 0.5, \
'pwr_spec_ZA_z10.dat' u 1:2 t 'ZA, z = 10' pt 1 ps 0.5, \
P_f(x) t 'P_f(k)', P_i(x) t 'P_i(k)'
#'pwr_spec_ZA_b100.dat' u 1:2 t 'ZA, b = 1,0' pt 1 ps 0.5, \

set ylabel "[P(k)-P_{lim}(k)]/P_{lin}(k)"
set output 'pwr_spec_diff_128p_128m_1024b_ZA.png'
unset logscale y
unset format y
set yrange [-0.4:0.1]


plot 'pwr_spec_diff_ZA_z200.dat' u 1:2 t 'ZA, z = 200' pt 1 ps 0.5, \
'pwr_spec_diff_ZA_z39.dat' u 1:2 t 'ZA, z = 39' pt 1 ps 0.5, \
'pwr_spec_diff_ZA_z21.dat' u 1:2 t 'ZA, z = 21' pt 1 ps 0.5, \
'pwr_spec_diff_ZA_z14.dat' u 1:2 t 'ZA, z = 14' pt 1 ps 0.5, \
'pwr_spec_diff_ZA_z10.dat' u 1:2 t 'ZA, z = 10' pt 1 ps 0.5, \
f_0(x), f_1(x)
#'pwr_spec_diff_ZA_b100.dat' u 1:2 t 'ZA, b = 1,0' pt 1 ps 0.5, \



unset logscale xy
set autoscale
set xlabel "x [Mpc/h]"
set ylabel "y [Mpc/h]"
set pm3d map
set cbrange [-0.5:3.5]
set output 'dens_map_ZA_z200.png'
splot 'rho_map_ZA_z200.dat' u 1:2:3 with pm3d t 'ZA, density, z = 200'
set output 'dens_map_ZA_z65.png'
splot 'rho_map_ZA_z65.dat' u 1:2:3 with pm3d t 'ZA, density, z = 65'
set output 'dens_map_ZA_z39.png'
splot 'rho_map_ZA_z39.dat' u 1:2:3 with pm3d t 'ZA, density, z = 39'
set output 'dens_map_ZA_z27.png'
splot 'rho_map_ZA_z27.dat' u 1:2:3 with pm3d t 'ZA, density, z = 27'
set output 'dens_map_ZA_z21.png'
splot 'rho_map_ZA_z21.dat' u 1:2:3 with pm3d t 'ZA, density, z = 21'
set output 'dens_map_ZA_z17.png'
splot 'rho_map_ZA_z17.dat' u 1:2:3 with pm3d t 'ZA, density, z = 17'
set output 'dens_map_ZA_z14.png'
splot 'rho_map_ZA_z14.dat' u 1:2:3 with pm3d t 'ZA, density, z = 14'
set output 'dens_map_ZA_z12.png'
splot 'rho_map_ZA_z12.dat' u 1:2:3 with pm3d t 'ZA, density, z = 12'
set output 'dens_map_ZA_z10.png'
splot 'rho_map_ZA_z10.dat' u 1:2:3 with pm3d t 'ZA, density, z = 10'