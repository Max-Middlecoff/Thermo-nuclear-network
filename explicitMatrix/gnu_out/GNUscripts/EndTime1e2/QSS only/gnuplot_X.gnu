set macro  # Enable macro definition

# Some macro definitions

# Colors: hex #RRGGBB
label_color = "#000000"
tic_color = "#000000"
title_color = "#000000"
myblue = "#2196f3"
myred = "#f44336"
mygreen = "#2e7d32"
mypurple = "#9c27b0"
myviolet = "#673ab7"
mybrown = "#795548"
myorange = "#ff9800"

# Width and height of postscript figure in inches
width = 8.5
height = 4.5

# x-axis resolution
set samples 2000

# Line styles.  
# For lines: plot x with lines ls 1
# For points: plot x with points ls 1

set style line 1 lc rgb 'black' pt 5   # fill square
set style line 2 lc rgb myred pt 7   # circle
set style line 3 lc rgb 'blue' pt 9   # triangle
set style line 4 lc rgb mygreen pt 7   # circle
set style line 5 lc rgb mypurple pt 7   # circle
set style line 6 lc rgb myviolet  pt 7   # circle
set style line 7 lc rgb mybrown  pt 7   # circle
set style line 8 lc rgb myorange  pt 7   # circle
set style line 9 lc rgb myblue pt 7   # circle
set style line 10 lc rgb 'green' pt 6   # open circle
set style line 11 lc rgb 'red' pt 4   # open circle
set style line 12 lc rgb 'gray' pt 4   # open circle
set style line 13 lc rgb 'gold' pt 4   # open circle
set style line 14 lc rgb 'dark-grey' pt 4   # open circle
set style line 15 lc rgb 'dark-cyan' pt 4   # open circle
set style line 16 lc rgb 'purple' pt 4   # open circle
set style line 17 lc rgb 'orchid' pt 4   # open circle
set style line 18 lc rgb 'brown' pt 4   # open circle
set style line 19 lc rgb 'dark-green' pt 4   # open circle
set style line 20 lc rgb 'magenta' pt 4   # open circle

#set xtics rotate        # Rotates x tic numbers by 90 degrees
#set ytics rotate        # Rotates y tic numbers by 90 degrees

# Set tic labeling with color
set xtics textcolor rgb tic_color
set ytics textcolor rgb tic_color

set bmargin 4  # Bottom margin

# Set screen display to same aspect ratio as postscript plot
set size ratio height/width

set xlabel 'Log time (s)' textcolor rgb tic_color font "Arial,32"
set ylabel 'Log X' textcolor rgb tic_color font "Arial,32"

# Uncomment following to set log or log-log plots
#set logscale x
#set logscale y

set pointsize 1.5    # Size of the plotted points

set key top outside   # Move legend to outside top
#unset key            # Don't show legend

#set timestamp       # Date/time

ds="QSS T9=5 rho=1e8"
ds = ds.": mass fraction"
set title ds textcolor rgb title_color

file1 = "gnufile.data"


# -------- Axis ranges and ticmarks -----------

xlow = -18
xup = 2
xtics = 1     # Space between major x ticmarks
minxtics = 5  # Number minor x tics

ylow = -14
yup = 0
ytics = 1.0      # Space between major y ticmarks
minytics = 5  # Number minor y tics

set xrange [xlow : xup]
set xtics  xlow, xtics, xup
set mxtics minxtics   # minor x tics per major tic

set yrange[ylow : yup]
set ytics ylow, ytics, yup
set mytics minytics   # minor y tics per major tic

#set grid   # set x-y grid at major ticmarks

# -------- Axis ranges and ticmarks -----------


# Edit the following plot commands to correspond to data
# read in from data file

plot file1 using 1:8 with lines ls 2 title "4He"
replot file1 using 1:9 with lines ls 3 title "12C"
replot file1 using 1:10 with lines ls 1 title "16O"
replot file1 using 1:11 with lines ls 4 title "20Ne"
replot file1 using 1:12 with lines ls 5 title "24Mg"
replot file1 using 1:13 with lines ls 6 title "28Si"
replot file1 using 1:14 with lines ls 7 title "32S"
replot file1 using 1:15 with lines ls 8 title "36Ar"
replot file1 using 1:16 with lines ls 9 title "40Ca"
replot file1 using 1:17 with lines ls 10 title "44Ti"
replot file1 using 1:18 with lines ls 11 title "48Cr"
replot file1 using 1:19 with lines ls 12 title "52Fe"
replot file1 using 1:20 with lines ls 13 title "56Ni"
replot file1 using 1:21 with lines ls 14 title "60Zn"
replot file1 using 1:22 with lines ls 1 title "64Ge"
replot file1 using 1:23 with lines ls 2 title "68Se"
#replot file1 using 1:24 with lines ls 3 title "17"

# Reset font sizes for .eps and .png output2

set title ds textcolor rgb title_color font "Arial,22"
set key top right font "Arial,22"
set xlabel 'Log time (s)' textcolor rgb tic_color font "Arial,28"
set ylabel 'Log X' textcolor rgb tic_color font "Arial,28"

# Plot to postscript file

set out "gnuplot_X.eps"    # Output file
set terminal postscript eps size width, height enhanced color solid lw 2 "Arial" 24
replot               # Plot to postscript file

# Plot to PNG file

set out "gnuplot_X.png"
# Assume 72 pixels/inch and make bitmap twice as large for display resolution
set terminal png transparent size 2*width*72, 2*height*72 lw 2
replot

quit
