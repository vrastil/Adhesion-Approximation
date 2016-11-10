
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

set output 'cut_FF_b0.0.png'
plot 'par_cut_FF_z100.dat' u 1:2 pt 7 ps 0.3
set output 'cut_FF_b0.4.png'
plot 'par_cut_FF_b40.dat' u 1:2 pt 7 ps 0.3
set output 'cut_FF_b0.6.png'
plot 'par_cut_FF_b80.dat' u 1:2 pt 7 ps 0.3
set output 'cut_FF_b1.0.png'
plot 'par_cut_FF_b100.dat' u 1:2 pt 7 ps 0.3

set terminal pngcairo size 640,480

set xlabel "k [h/Mpc]"
set ylabel "P(k) [(Mpc/h)^3]"
set output 'pwr_spec_128p_128m_4096b_FF.png'
set autoscale
set logscale xy
set format y "10^{%L}"

plot 'pwr_spec_FF_z100.dat' u 1:2 t 'FF, z = 200' pt 1 ps 0.5, \
'pwr_spec_FF_b40.dat' u 1:2 t 'FF, b = 0,4' pt 1 ps 0.5, \
'pwr_spec_FF_b80.dat' u 1:2 t 'FF, b = 0,8' pt 1 ps 0.5, \
'pwr_spec_FF_b100.dat' u 1:2 t 'FF, b = 1,0' pt 1 ps 0.5
# P_0(x) t 'P_0(k)', P_i(x) t 'P_i(k)'

set ylabel "[P(k)-P_{lim}(k)]/P_{lin}(k)"
set output 'pwr_spec_diff_128p_128m_4096b_FF.png'
unset logscale y
unset format y
set yrange [-1.8:0.2]


plot 'pwr_spec_diff_FF_z100.dat' u 1:2 t 'FF, z = 200' pt 1 ps 0.5, \
'pwr_spec_diff_FF_b40.dat' u 1:2 t 'FF, b = 0,4' pt 1 ps 0.5, \
'pwr_spec_diff_FF_b80.dat' u 1:2 t 'FF, b = 0,8' pt 1 ps 0.5, \
'pwr_spec_diff_FF_b100.dat' u 1:2 t 'FF, b = 1,0' pt 1 ps 0.5, \
f_0(x), f_1(x)