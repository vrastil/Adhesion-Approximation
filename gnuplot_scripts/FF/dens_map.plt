
cd 'output/run'
set terminal pngcairo size 960,720

unset logscale xy
set autoscale
set xlabel "x [Mpc/h]"
set ylabel "y [Mpc/h]"
set pm3d map
set cbrange [-0.5:3.5]
set output 'dens_map_FF_z200.png'
splot 'rho_map_FF_z200.dat' u 1:2:3 with pm3d t 'FF, density, z = 200'
set output 'dens_map_FF_z65.png'
splot 'rho_map_FF_z65.dat' u 1:2:3 with pm3d t 'FF, density, z = 65'
set output 'dens_map_FF_z39.png'
splot 'rho_map_FF_z39.dat' u 1:2:3 with pm3d t 'FF, density, z = 39'
set output 'dens_map_FF_z27.png'
splot 'rho_map_FF_z27.dat' u 1:2:3 with pm3d t 'FF, density, z = 27'
set output 'dens_map_FF_z21.png'
splot 'rho_map_FF_z21.dat' u 1:2:3 with pm3d t 'FF, density, z = 21'
set output 'dens_map_FF_z17.png'
splot 'rho_map_FF_z17.dat' u 1:2:3 with pm3d t 'FF, density, z = 17'
set output 'dens_map_FF_z14.png'
splot 'rho_map_FF_z14.dat' u 1:2:3 with pm3d t 'FF, density, z = 14'
set output 'dens_map_FF_z12.png'
splot 'rho_map_FF_z12.dat' u 1:2:3 with pm3d t 'FF, density, z = 12'
set output 'dens_map_FF_z10.png'
splot 'rho_map_FF_z10.dat' u 1:2:3 with pm3d t 'FF, density, z = 10'