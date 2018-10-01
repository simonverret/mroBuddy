w = 600
h = 600

set term wxt persist size w,h
bind "w" "unset output; exit gnuplot"
# set terminal epslatex size w cm,h cm color dashlength 0.3
# set output 'plot.tex'
# set termopt enhanced 

set linestyle 1 lc rgb 'red' lw 0.5 dt 1 ps 0.3 pt 7

set title "$t=0.190$ eV, $\\hbar/\\tau=t/25$, $t_z=-0.07t$"

x="(column('theta'))"
set xlabel 'theta'
set xrange [0:pi]
set xtics ('0' 0,'45' pi/4,'90' pi/2,'135' pi/3, '180' pi)

y="(column('sigma_ref')/column('sigma_zz'))"
set ylabel '\rho_zz' offset -1,0
#set yrange [0.8:1.2]
set ytics 0.1

unset key
# plot \
# for [ii=240:247]'FS.dat' i ii u 1:2 w l  lw 1.5 lc rgb 'black',\
# 'amro.dat' u @x:@y w l lw 1 lc rgb 'red' dt 3

plot \
'amro.dat' u @x:@y w l lw 2 lc rgb 'blue',\

#'fastAMRO.dat' u @x:@y w l lw 1 lc rgb 'red'

pause -1