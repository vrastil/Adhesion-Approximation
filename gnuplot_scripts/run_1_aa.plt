
cd 'output/run_1'
set terminal pngcairo

A = 187826
n = 1
z_in = 200

k_G = 0.5
k2_G = k_G**2

a = 6.4
b = 3.0
c = 1.7
v = 1.13

Omega_0 = 1
h = 0.67
q(x) = x / (Omega_0*h)

P(x) = A*(x**n)
T(x) = log(1+2.34*q(x))/(2.34*q(x))*(1 + 3.89*q(x) + (16.1*q(x))**2 + (5.4*q(x))**3 + (6.71*q(x))**4)**(-1./4.);
T_k_G(x) = exp(-x**2/(2*k2_G))

P_0(x) = P(x)*T(x)*T(x)
P_i(x) = P_0(x)/(z_in+1)**2
P_ti(x) = P_0(x)*T_k_G(x)/(z_in+1)**2

set xlabel "x [Mpc/h]"
set ylabel "y [Mpc/h]"

set output 'cut_z200.png'
plot 'par_cut_z200.dat' u 1:2 w dots

set output 'cut_z63.png'
plot 'par_cut_z63.dat' u 1:2 w dots

set output 'cut_z19.png'
plot 'par_cut_z19.dat' u 1:2 w dots

set output 'cut_z5.png'
plot 'par_cut_z5.dat' u 1:2 w dots

set output 'cut_z1.png'
plot 'par_cut_z1.dat' u 1:2 w dots

set xlabel "k [h/Mpc]"
set ylabel "P(k) [(Mpc/h)^3]"

set output 'pwr_spec_AA.png'
set autoscale
set logscale xy
set format y "10^{%L}"

plot 'pwr_spec_z200.dat' u 1:2 t 'z = 200' pt 1 ps 0.5, \
'pwr_spec_z63.dat' u 1:2 t 'z = 63' pt 1 ps 0.5, \
'pwr_spec_z19.dat' u 1:2 t 'z = 19' pt 1 ps 0.5, \
'pwr_spec_z5.dat' u 1:2 t 'z = 5' pt 1 ps 0.5, \
'pwr_spec_z1.dat' u 1:2 t 'z = 1' pt 1 ps 0.5, \
P_0(x) t 'P_0(k)', P_i(x) t 'P_i(k)', P_ti(x) t 'P^t_i(k)'
