set output plot1.png
set term png size 1920, 1080
set size 0.5, 1
set contour
unset surface
set view map
splot '-' using 1:2:3 with lines
0 0 0
0 1 0
0 2 0

1 0 2
1 1 2
1 2 2

2 0 1
2 1 1
2 2 1
