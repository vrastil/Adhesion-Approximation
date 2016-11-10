
cd 'output/run_1'
set terminal pngcairo size 1280,960

A = 187826
n = 1
z_in = 200

Omega_0 = 1
h = 0.67
q(x) = x / (Omega_0*h)

P(x) = A*(x**n)
T(x) = log(1+2.34*q(x))/(2.34*q(x))*(1 + 3.89*q(x) + (16.1*q(x))**2 + (5.4*q(x))**3 + (6.71*q(x))**4)**(-1./4.);

P_0(x) = P(x)*T(x)*T(x)
P_i(x) = P_0(x)/(z_in+1)**2

f_0(x)= 0
f_1(x)= -0.2

set xlabel "x [Mpc/h]"
set ylabel "y [Mpc/h]"

set output 'cut.png'
plot 'par_cut_ZA_b100.dat' u 1:2 pt 7 ps 0.3

set terminal pngcairo size 640,480

set xlabel "k [h/Mpc]"
set ylabel "P(k) [(Mpc/h)^3]"
#set output 'SW_pwr_spec_b64.png'
set autoscale
set logscale xy
set format y "10^{%L}"

#plot 'pwr_spec_ZA_z50.dat' u 1:2 t 'SW, z = 50, b = 128' pt 1 ps 0.5