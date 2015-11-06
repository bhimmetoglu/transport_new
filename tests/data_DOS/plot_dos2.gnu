set term post eps enhanced color 24
set output "DOS2-32.eps"
set xrange [11:15]
#set yrange [0:10]
set xtics 1
set mxtics 2
set ytics 2
set mytics 2
set ylabel 'DOS (eV^{-1})'
set xlabel 'E (eV)'
set style line 1 lt 1 lc 1 lw 3
set style line 2 lt 1 lc 2 lw 3
set style line 3 lt 2 lc 3 lw 3
set style line 4 lt 2 lc 4 lw 3
set style line 5 lt 1 lc 5 lw 3
set style line 6 lt 2 lc 7 lw 3
#set pointsize 1.0
set key top left
plot "dos.32-1.4.out" u 1:2 w l ls 1 t 'adaptive,a=1.4',\
     "dos.32-1.0.out" u 1:2 w l ls 3 t 'adaptive,a=1.0',\
     "dos.32-0.8.out" u 1:2 w l ls 6 t 'adaptive,a=0.8'
