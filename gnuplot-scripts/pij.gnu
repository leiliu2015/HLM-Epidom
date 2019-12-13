# Hi-C results of A-23/I-14/R-11
fa = '../Hi-C/Kc167-5kb.chr2R.19600000-20100000.cm'
fb = '../Hi-C/Kc167-5kb.chr2R.15700000-16200000.cm'
fc = '../Hi-C/Kc167-5kb.chr3R.02450000-02950000.cm'
# HLM results of A-23/I-14/R-11
ga = '../A-23/gnm.t0-0.rc2.2.cm'
gb = '../I-14/gnm.t0-0.rc2.2.cm'
gc = '../R-11/gnm.t0-0.rc2.2.cm'

set size square
set tics front out scale 0.75

sft= 0.5
rs = 0.005
xmx= 100*rs
set xrange [0:xmx]
set yrange [0:xmx]
set xtics 0,0.2,1 offset 0,0.2 nomirror
set mxtics 2
set ytics 0,0.2,1 offset 0.4,0 nomirror
set mytics 2

set cbrange [-2.0:-0.5]; set cbtics -3,1,0 offset -0.8,0
set cbtics add('10^{-1}'-1,'10^{-2}'-2)
set cbtics in

set label 1 'Hi-C' at first 0.03, 0.47 front
set label 2 'HLM' at first 0.40, 0.05 front

set term wxt 0 size 900, 360
set multiplot layout 1,3

# A-23
set title 'A-23' offset 0, -0.5
set xtics add('19.6'0,'19.8'0.2,'20'0.4)
set ytics add('19.6'0,'19.8'0.2,'20'0.4)
set ytics add(''0,''0.2,''0.4)
set palette defined ( 0 '#FFFFFF',\
                      1 '#FFF5F0',\
    	    	      2 '#FEE0D2',\
		      3 '#FCBBA1',\
		      4 '#FC9272',\
		      5 '#FB6A4A',\
		      6 '#EF3B2C',\
		      7 '#CB181D',\
		      8 '#99000D' )
plot fa matrix u (($2+sft)*rs):(($1+sft)*rs):($3>0 && $2<$1 ? log10($3) : NaN) w image notitle, \
     ga matrix u (($2+sft)*rs):(($1+sft)*rs):($3>0 && $2>$1 ? log10($3) : NaN) w image notitle

# I-14
set title 'I-14' offset 0, -0.5
set xtics add('15.7'0,'15.9'0.2,'16.1'0.4)
set ytics add('15.7'0,'15.9'0.2,'16.1'0.4)
set ytics add(''0,''0.2,''0.4)
set palette defined ( 0 '#FFFFFF',\
    	    	      1 '#F0F0F0',\
		      2 '#D9D9D9',\
		      3 '#BDBDBD',\
		      4 '#969696',\
		      5 '#737373',\
		      6 '#525252',\
		      7 '#252525' )
plot fb matrix u (($2+sft)*rs):(($1+sft)*rs):($3>0 && $2<$1 ? log10($3) : NaN) w image notitle, \
     gb matrix u (($2+sft)*rs):(($1+sft)*rs):($3>0 && $2>$1 ? log10($3) : NaN) w image notitle

# R-11
set title 'R-11' offset 0, -0.5
set xtics add('2.45'0,'2.65'0.2,'2.85'0.4)
set ytics add('2.45'0,'2.65'0.2,'2.85'0.4)
set ytics add(''0,''0.2,''0.4)
set palette defined ( 0 '#FFFFFF',\
                      1 '#F7FBFF',\
    	    	      2 '#DEEBF7',\
		      3 '#C6DBEF',\
		      4 '#9ECAE1',\
		      5 '#6BAED6',\
		      6 '#4292C6',\
		      7 '#2171B5',\
		      8 '#084594' )
plot fc matrix u (($2+sft)*rs):(($1+sft)*rs):($3>0 && $2<$1 ? log10($3) : NaN) w image notitle, \
     gc matrix u (($2+sft)*rs):(($1+sft)*rs):($3>0 && $2>$1 ? log10($3) : NaN) w image notitle

#
unset multiplot

