set term post eps enhanced color 24
set output "sig_vs_cut.eps"
#set xrange [0.4:1.6]
set yrange [350:850]
set ytics 100
set mytics 2
set ylabel '{/Symbol s} ({/Symbol W} cm)^{-1}'
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
plot "sig_vs_cut.out" u 1:2 w lp ls 1 notitle
