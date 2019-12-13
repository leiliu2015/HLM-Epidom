fa = '../A-23/gnm.t0-0.2016zhuang418.rgs'
fb = '../I-14/gnm.t0-0.2016zhuang418.rgs'
fc = '../R-11/gnm.t0-0.2016zhuang418.rgs'

mps = 1.5
mlw = 1.5
mlwf= 2.5
fq  = 1
set style line 1 lt 1 lc rgb '#FB9A99' lw mlw pt 6  ps 1.1*mps # light red
set style line 2 lt 1 lc rgb '#969696' lw mlw pt 4  ps 1.0*mps # light gray
set style line 3 lt 1 lc rgb '#A6CEE3' lw mlw pt 12 ps 1.4*mps # light blue
set style line 4 lt 1 lc rgb '#CB181D' lw mlw pt 6  ps 1.1*mps # red
set style line 5 lt 1 lc rgb '#000000' lw mlw pt 4  ps 1.0*mps # black
set style line 6 lt 1 lc rgb '#2171B5' lw mlw pt 12 ps 1.4*mps # blue
set tics front scale 0.75
set border 3
set key left top samplen 0.5 spacing 1.2
set size square

rs = 0.005
set xrange [0.1:0.45]; set xtics add(0.2,0.4) nomirror
set yrange [2:4.2]; set ytics add(2,3,4) nomirror

set xlabel 's [Mb]'
set ylabel 'r_{g}(s) [a]'
set grid xtics lt 0 lc rgb 'gray'
set grid ytics lt 0 lc rgb 'gray'
set logscale

ga(x) = (x>=0.1) && (x<0.5) ? (10**(-0.226)) * ((x/rs)**0.429) : NaN
gb(x) = (x>=0.1) && (x<0.5) ? (10**(-0.026)) * ((x/rs)**0.260) : NaN
gc(x) = (x>=0.1) && (x<0.5) ? (10**(-0.144)) * ((x/rs)**0.346) : NaN


plot fa every fq u ($1*rs):2 w p ls 1 t 'A-23', \
     fb every fq u ($1*rs):2 w p ls 2 t 'I-14', \
     fc every fq u ($1*rs):2 w p ls 3 t 'R-11', \
     ga(x) w l ls 4 lw mlwf t 's^{0.43}', \
     gb(x) w l ls 5 lw mlwf t 's^{0.26}', \
     gc(x) w l ls 6 lw mlwf t 's^{0.35}'

