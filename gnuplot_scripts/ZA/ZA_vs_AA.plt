
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

set output 'cut_AAd_b0.png'
plot 'par_cut_AAd_b0.dat' u 1:2 w points pt 7 ps 0.5

set output 'cut_AAd_b2.png'
plot 'par_cut_AAd_b2.dat' u 1:2 w points pt 7 ps 0.5

set output 'cut_AAd_b4.png'
plot 'par_cut_AAd_b4.dat' u 1:2 w points pt 7 ps 0.5

set output 'cut_AAd_b6.png'
plot 'par_cut_AAd_b6.dat' u 1:2 w points pt 7 ps 0.5

set output 'cut_AAd_b8.png'
plot 'par_cut_AAd_b8.dat' u 1:2 w points pt 7 ps 0.5

set output 'cut_AAd_b10.png'
plot 'par_cut_AAd_b10.dat' u 1:2 w points pt 7 ps 0.5

set output 'cut_ZAq_b0.png'
plot 'par_cut_ZAq_b0.dat' u 1:2 w points pt 7 ps 0.5

set output 'cut_ZAq_b2.png'
plot 'par_cut_ZAq_b2.dat' u 1:2 w points pt 7 ps 0.5

set output 'cut_ZAq_b4.png'
plot 'par_cut_ZAq_b4.dat' u 1:2 w points pt 7 ps 0.5

set output 'cut_ZAq_b6.png'
plot 'par_cut_ZAq_b6.dat' u 1:2 w points pt 7 ps 0.5

set output 'cut_ZAq_b8.png'
plot 'par_cut_ZAq_b8.dat' u 1:2 w points pt 7 ps 0.5

set output 'cut_ZAq_b10.png'
plot 'par_cut_ZAq_b10.dat' u 1:2 w points pt 7 ps 0.5

set terminal pngcairo size 800, 600

set xlabel "k [h/Mpc]"
set ylabel "P(k) [(Mpc/h)^3]"

set output 'pwr_spec_q_vs_k_diff.png'
set autoscale
set logscale xy
set format y "10^{%L}"

plot 'pwr_spec_AAd_b10.dat' u 1:2 t 'AAd, b = 10' pt 1 ps 0.5, \
'pwr_spec_ZAq_b10.dat' u 1:2 t 'ZAq, b = 10' pt 1 ps 0.5, \
'pwr_spec_AAd_b8.dat' u 1:2 t 'AAd, b = 8' pt 1 ps 0.5, \
'pwr_spec_ZAq_b8.dat' u 1:2 t 'ZAq, b = 8' pt 1 ps 0.5, \
'pwr_spec_AAd_b6.dat' u 1:2 t 'AAd, b = 6' pt 1 ps 0.5, \
'pwr_spec_ZAq_b6.dat' u 1:2 t 'ZAq, b = 6' pt 1 ps 0.5, \
'pwr_spec_AAd_b4.dat' u 1:2 t 'AAd, b = 4' pt 1 ps 0.5, \
'pwr_spec_ZAq_b4.dat' u 1:2 t 'ZAq, b = 4' pt 1 ps 0.5, \
'pwr_spec_AAd_b2.dat' u 1:2 t 'AAd, b = 2' pt 1 ps 0.5, \
'pwr_spec_ZAq_b2.dat' u 1:2 t 'ZAq, b = 2' pt 1 ps 0.5, \
'pwr_spec_AAd_b0.dat' u 1:2 t 'AAd, b = 0' pt 1 ps 0.5, \
'pwr_spec_ZAq_b0.dat' u 1:2 t 'ZAq, b = 0' pt 1 ps 0.5, \
P_0(x) t 'P_0(k)', P_i(x) t 'P_i(k)', P_ti(x) t 'P^t_i(k)'
