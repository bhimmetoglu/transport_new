set term post eps enhanced color 24
set output "DOS1.eps"
set xrange [4:15]
set yrange [0:14]
set xtics 2
set mxtics 2
set ytics 2
set mytics 2
set ylabel 'DOS (eV^{-1})'
set xlabel 'E (eV)'
#set arrow from 0,0.01 to 0,10000 nohead
set style line 1 lt 1 lc 1 lw 3
set style line 2 lt 1 lc 2 lw 3
set style line 3 lt 1 lc 3 lw 3
set style line 4 lt 2 lc 4 lw 3
set style line 5 lt 1 lc 5 lw 3
set style line 6 lt 2 lc 7 lw 3
#set pointsize 1.0
#set key top left
plot "dos.32-1.4.out" u 1:2 w l ls 1 t 'adaptive,a=1.4',\
     "dos.0.2.out" u 1:2 w l ls 6 t '{/Symbol s}=0.2 eV'
