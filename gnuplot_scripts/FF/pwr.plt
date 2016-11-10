
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


set terminal pngcairo size 640,480

set xlabel "k [h/Mpc]"
set ylabel "P(k) [(Mpc/h)^3]"
set output 'pwr_spec_128p_128m_1024_FF.png'
set autoscale
set logscale xy
set format y "10^{%L}"

plot 'pwr_spec_FF_z200.dat' u 1:2 t 'FF, z = 200' pt 1 ps 0.5, \
'pwr_spec_FF_z39.dat' u 1:2 t 'FF, z = 39' pt 1 ps 0.5, \
'pwr_spec_FF_z21.dat' u 1:2 t 'FF, z = 21' pt 1 ps 0.5, \
'pwr_spec_FF_z14.dat' u 1:2 t 'FF, z = 14' pt 1 ps 0.5, \
'pwr_spec_FF_z10.dat' u 1:2 t 'FF, z = 10' pt 1 ps 0.5, \
P_f(x) t 'P_f(k)', P_i(x) t 'P_i(k)'
#'pwr_spec_FF_b100.dat' u 1:2 t 'FF, b = 1,0' pt 1 ps 0.5, \