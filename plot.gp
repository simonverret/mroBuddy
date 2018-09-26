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
# set xrange [-3.1416:3.1416]
# set xtics .5

y="(column('ky'))"
set ylabel 'ky' offset -1,0
# set yrange [-3.1416:3.1416]
# set ytics .5

z="(column('kz'))"
set zlabel 'kz' offset -1,0
# set zrange [-6.283:6.283]
# set ztics 1.

set ticslevel 0
#set key top left samplen 2 spacing 0.8 maxrows 2

splot \
'amro.dat' u @x:@y:@z w lp ls 1 t "test"


pause mouse close
exit