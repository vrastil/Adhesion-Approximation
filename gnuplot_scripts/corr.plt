
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

P_0(x) = P(x)*T(x)*T(x)
P_i(x) = P_0(x)/(z_in+1)**2

set terminal pngcairo size 800, 600

set xlabel "k [h/Mpc]"
set ylabel "P(k) [(Mpc/h)^3]"

set output 'pwr_spec_w_vs_k_diff.png'
set autoscale
set logscale xy
set format y "10^{%L}"

plot 'pwr_spec_ZA_b0.dat' u 1:2 t 'ZA, b = 0' pt 1 ps 0.5, \
'pwr_spec_ZA_b4.dat' u 1:2 t 'ZA, b = 4' pt 1 ps 0.5, \
'pwr_spec_ZA_b10.dat' u 1:2 t 'ZA, b = 10' pt 1 ps 0.5, \
'pwr_spec_ZAopt0_b0.dat' u 1:2 t 'ZAopt0, b = 0' pt 1 ps 0.5, \
'pwr_spec_ZAopt0_b4.dat' u 1:2 t 'ZAopt0, b = 4' pt 1 ps 0.5, \
'pwr_spec_ZAopt0_b10.dat' u 1:2 t 'ZAopt0, b = 10' pt 1 ps 0.5, \
'pwr_spec_ZAopt1_b0.dat' u 1:2 t 'ZAopt1, b = 0' pt 1 ps 0.5, \
'pwr_spec_ZAopt1_b4.dat' u 1:2 t 'ZAopt1, b = 4' pt 1 ps 0.5, \
'pwr_spec_ZAopt1_b10.dat' u 1:2 t 'ZAopt1, b = 10' pt 1 ps 0.5, \
'pwr_spec_ZAopt2_b0.dat' u 1:2 t 'ZAopt2, b = 0' pt 1 ps 0.5, \
'pwr_spec_ZAopt2_b4.dat' u 1:2 t 'ZAopt2, b = 4' pt 1 ps 0.5, \
'pwr_spec_ZAopt2_b10.dat' u 1:2 t 'ZAopt2, b = 10' pt 1 ps 0.5, \
'pwr_spec_ZAopt3_b0.dat' u 1:2 t 'ZAopt3, b = 0' pt 1 ps 0.5, \
'pwr_spec_ZAopt3_b4.dat' u 1:2 t 'ZAopt3, b = 4' pt 1 ps 0.5, \
'pwr_spec_ZAopt3_b10.dat' u 1:2 t 'ZAopt3, b = 10' pt 1 ps 0.5, \
P_i(x) t 'P_i(k)', P_0(x) t 'P_0(k)'