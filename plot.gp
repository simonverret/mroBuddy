w = 600
h = 600

set term wxt persist size w,h
# set terminal epslatex size w cm,h cm color dashlength 0.3
# set output 'plot.tex'
# set termopt enhanced 

set linestyle 1 lc rgb 'red' lw 0.5 dt 1 ps 0.3 pt 7

set title "$t=0.190$ eV, $\\hbar/\\tau=t/25$, $t_z=-0.07t$"

x="(column('kx'))"
set xlabel 'kx'
set xrange [-3.1416:3.1416]
set xtics .5

y="(column('ky'))"
set ylabel 'ky' offset -1,0
set yrange [-3.1416:3.1416]
set ytics .5

z="(column('kz'))"
set zlabel 'kz' offset -1,0
set zrange [-6.283:6.283]
set ztics 1.

set ticslevel 0
unset key

# plot \
# for [ii=240:247]'FS.dat' i ii u 1:2 w l  lw 1.5 lc rgb 'black',\
# 'amro.dat' u @x:@y w l lw 1 lc rgb 'red' dt 3

splot \
'FS.dat' u @x:@y:@z w lp pt 7 ps 0.1  lw 0.5 lc rgb 'black',\
'amro.dat' u @x:@y:@z w l lw 2 lc rgb 'red'

pause mouse close
exit