
cd 'output/run_1'
set terminal png

A = 0.0333428
n = 1

k_G = 0.69
k2_G = k_G**2

a = 6.4
b = 3.0
c = 1.7
v = 1.13

P_i(x) = A*(x**n)

T(x) = (1+(a*x+(b*x)**(3/2.0)+(c*x)**(2.0))**v)**(-2.0/v)

T_k_G(x) = exp(-x**2/(2*k2_G))
P(x) = P_i(x)*T(x)
P_G(x) = P(x)*T_k_G(x)

set xlabel "x [Mpc/h]"
set ylabel "y [Mpc/h]"

set output 'cut_ZA.png'
plot 'par_cut_z200.dat' u 1:2 w dots


set xlabel "k [h/Mpc]"
set ylabel "P(k) [(Mpc/h)^3]"
set output 'pwr_spec_ZA.png'
set autoscale
set logscale xy
set format y "10^{%L}"

plot 'pwr_spec_ZA.dat' u 1:2 t 'z = 200' pt 1 ps 0.5, P(x) t 'Initial P(k)'