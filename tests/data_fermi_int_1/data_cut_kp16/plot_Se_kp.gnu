set term post eps enhanced color 24
set output "Se_vs_cut.eps"
#set xrange [0.4:1.6]
#set yrange [0.0:2.0e+3]
#set mxtics 2
#set ytics 500
set mytics 2
set ylabel '|S| ({/Symbol m}V/K)'
set xlabel 'cut'
#set arrow from 0,0.01 to 0,10000 nohead
set style line 1 lt 1 lc 1 lw 3
set style line 2 lt 1 lc 2 lw 3
set style line 3 lt 1 lc 3 lw 3
set style line 4 lt 2 lc 4 lw 3
set style line 5 lt 1 lc 5 lw 3
set style line 6 lt 2 lc 7 lw 3
#set pointsize 1.0
set key top left
plot "sig_vs_cut.out" u 1:(-1e6*$3) w lp ls 1 notitle
