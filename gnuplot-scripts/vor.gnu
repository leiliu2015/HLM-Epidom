d(i) = i==1 ? 'A-23' : (i==2 ? 'I-14' : 'R-11')
f(i) = sprintf("../%s/gnm.t0-0.2016zhuang418.", d(i))

set style line 1 lt 1 lc rgb '#EF3B2C' # red
set style line 2 lt 1 lc rgb '#969696' # grey
set style line 3 lt 1 lc rgb '#084594' # blue
set tics back scale 0.75
unset key

unset xlabel
set xrange [0:3]
set xtics ('A-23'0.5, 'I-14'1.5, 'R-11'2.5) rotate by -90 nomirror
set xtics out scale 0.4 offset 0, +0.2

mps = 1
mlw = 3
set border 3
set boxwidth 0.5 absolute
set style fill solid 0.2
set term wxt size 500,600
set multiplot layout 2, 2

#
apx = 'irg'
set ylabel 'L/R@^{3}_{g} [a^{-3}]'; set yrange [0:5.95]; set ytics 0,2,10 nomirror; set mytics 2
plot for [i=1:3] f(i).apx u (i-0.5):12:11:15:14 w candlesticks ls i lw mlw notitle, \
     for [i=1:3] f(i).apx u (i-0.5):13:13:13:13 w candlesticks ls i lw mlw notitle, \
     for [i=1:3] f(i).apx u (i-0.5):8 w p ls i lw mlw pt 2 ps mps notitle

#
apx = 'asp'
set ylabel 'Asp.'; set yrange [0:0.85]; set ytics 0,0.3,2.0 nomirror; set mytics 3
plot for [i=1:3] f(i).apx u (i-0.5):12:11:15:14 w candlesticks ls i lw mlw notitle, \
     for [i=1:3] f(i).apx u (i-0.5):13:13:13:13 w candlesticks ls i lw mlw notitle, \
     for [i=1:3] f(i).apx u (i-0.5):8 w p ls i lw mlw pt 2 ps mps notitle

#
apx = 'rho'
set ylabel 'Density [a^{-3}]'; set yrange [0.3:0.49]; set ytics 0,0.1,1.0 nomirror; set mytics 2
plot for [i=1:3] f(i).apx u (i-0.5):12:11:15:14 w candlesticks ls i lw mlw notitle, \
     for [i=1:3] f(i).apx u (i-0.5):13:13:13:13 w candlesticks ls i lw mlw notitle, \
     for [i=1:3] f(i).apx u (i-0.5):8 w p ls i lw mlw pt 2 ps mps notitle

#
apx = 'srf'
set ylabel 'S/S_{0}'; set yrange [1.5:2.65]; set ytics 0,0.5,3.0 nomirror; set mytics 5
plot for [i=1:3] f(i).apx u (i-0.5):12:11:15:14 w candlesticks ls i lw mlw notitle, \
     for [i=1:3] f(i).apx u (i-0.5):13:13:13:13 w candlesticks ls i lw mlw notitle, \
     for [i=1:3] f(i).apx u (i-0.5):8 w p ls i lw mlw pt 2 ps mps notitle

#
unset multiplot

