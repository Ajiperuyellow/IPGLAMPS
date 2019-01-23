reset

! cat /proc/cpuinfo | awk -F":" '$1 ~ "model name" { print $2 }'|head -1 > name.txt
CPUNAME="`cat name.txt`"
print CPUNAME

set terminal push
set terminal postscript eps enhanced color 'Helvetica' 24 size 6.0,6.0 dl 2
set output "AAA.eps"


set yrange [0:]
set xrange [0:30]

set xtics 0,5,20
set mxtics 5

set label 1 CPUNAME at graph 0.1, graph 0.9
set label 2 "red: SSE" at graph 0.85, graph 0.2 center
set label 3 "green: trunk" at graph 0.85, graph 0.15 center

set xlabel "run"

set ylabel "cpu time/loop (TSC)"

set key right center

#plot "aaa.rep.setBetaCM" u ($5) w l lt 1 lc 1 lw 6 t "setBetaCM",\
#     "aaa.rep.setBetaCMScalar" u ($5) w l lt 1 lc 2 lw 6 t "setBetaCMScalar",\
#     "aaa.rep.boost2" u ($5) w l lt 1 lc 1 lw 6 t "boost2",\
#     "aaa.rep.boostScalar2" u ($5) w l lt 1 lc 2 lw 6 t "boostScalar2",\
#     "aaa.rep.boost" u ($5) w l lt 1 lc 1 lw 2 t "boost",\
#     "aaa.rep.boostScalar" u ($5) w l lt 1 lc 2 lw 2 t "boostScalar",\
#     "aaa.rep.total" u (-$5) w l lt 1 lc 3 lw 2 t "total"

plot "aaa.rep.all" u ($1+0.5):($3) w lp lt 1 lc 1 lw 6 t "setBetaCM",\
     "aaa.rep.all" u ($1+0.5):($4) w lp lt 1 lc 2 lw 6 t "",\
     "aaa.rep.all" u ($1+0.5):($5) w lp lt 2 lc 1 lw 6 t "boost2",\
     "aaa.rep.all" u ($1+0.5):($6) w lp lt 2 lc 2 lw 6 t "",\
     "aaa.rep.all" u ($1+0.5):($7) w lp lt 3 lc 1 lw 2 t "boost",\
     "aaa.rep.all" u ($1+0.5):($8) w lp lt 3 lc 2 lw 2 t "",\
     "aaa.rep.all" u ($1+0.5):(-$2) w lp lt 4 lc 3 lw 2 t "total"


set output
set terminal pop