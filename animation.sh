#!/bin/sh

rm *.png -f

gnuplot plot.gplt 

convert -delay 10 -loop 0 result/*.png animation.gif

firefox animation.gif

