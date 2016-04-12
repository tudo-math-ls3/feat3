#####################
# Defect / Iters
#####################
set xlabel "Iterartion"
set ylabel "Defect"
set style data linespoints
set key autotitle columnheader
set terminal postscript eps enhanced color font "Helvetica,16"

files = system("ls -1 defect-iters*.dat")
do for [file in files]{
  print file
  stats file nooutput
  set logscale y
  set output file.".eps"
  plot for [IDX=0:STATS_blocks-1] file i IDX u 1:4 with yerrorlines
  unset logscale
}

unset style
unset logscale


#####################
# time / ranks
#####################
reset

set xlabel "Rank"
set ylabel "Time [s]"
set xtics rotate out
set style data histogram
set style histogram errorbars lw 1
set style fill solid border
set style fill solid 0.3
set bars front
set terminal postscript eps enhanced color font "Helvetica,16"

files = system("ls -1 time-ranks*.dat")
do for [file in files]{
  print file
  set output file.".eps"
  plot for [COL=1:4:3] file using COL+1:COL+2:COL+3:xticlabels(1) title columnheader
}



#####################
# Global time / ranks
#####################
reset

set key outside
set xlabel "Rank"
set ylabel "Time [s]"
set xtics rotate out
set style data histogram
set style histogram rowstacked
set boxwidth 0.6 relative
set style fill solid border
set style fill solid 0.3
set bars front
set terminal postscript eps enhanced color font "Helvetica,16"

files = system("ls -1 global-time-ranks*.dat")
do for [file in files]{
  print file
  set output file.".eps"
  plot for [COL=2:6] file using COL:xticlabels(1) title columnheader
}



#####################
# Global size / ranks
#####################
reset

set key outside
set xlabel "Rank"
set ylabel "Size [MByte]"
set xtics rotate out
set style data histogram
set style histogram rowstacked
set boxwidth 0.6 relative
set style fill solid border
set style fill solid 0.3
set bars front
set terminal postscript eps enhanced color font "Helvetica,16"

files = system("ls -1 global-size-ranks*.dat")
do for [file in files]{
  print file
  set output file.".eps"
  plot for [COL=2:4] file using COL:xticlabels(1) title columnheader
}


#####################
# iters / ranks
#####################
reset

set key outside
set yrange [0:]
set xlabel "Rank"
set ylabel "Iterations"
set xtics rotate out
set style data histogram
set style histogram errorbars lw 1
set style fill solid border
set style fill solid 0.3
set bars front
set terminal postscript eps enhanced color font "Helvetica,16"

files = system("ls -1 iters-ranks*.dat")
do for [file in files]{
  print file
  set output file.".eps"
  plot file using 2:3:4:xticlabels(1) title columnheader
}


#####################
# time / levels
#####################
reset

set key outside
set xlabel "Level"
set ylabel "Time [s]"
set style data histogram
set style histogram errorbars lw 1
set style fill solid border
set style fill solid 0.3
set bars front
set key autotitle columnheader
set terminal postscript eps enhanced color font "Helvetica,16"

files = system("ls -1 time-levels*.dat")
do for [file in files]{
  print file
  set output file.".eps"
  plot for [COL=1:4:3] file using COL+1:COL+2:COL+3:xticlabels(1)
}
