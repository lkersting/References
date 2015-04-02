reset
clear
set terminal pngcairo nocrop enhanced font 'Verdana,10' size 1000,500 
set output 'Sphere1.png'
set title "Flux Ratio for cold Hydrogen sphere of radius 0.0005 cm" 
set xlabel "Energy (MeV)" 
set ylabel "Flux Ratio" 
set arrow from 0,1 to 0.01,1 nohead
plot "Sphere1.tsv" title "" with errorbars
clear
set terminal pngcairo nocrop enhanced font 'Verdana,10' size 1000,500 
set output 'Sphere2.png'
set title "Flux Ratio for cold Hydrogen sphere of radius 0.001 cm" 
set xlabel "Energy (MeV)" 
set ylabel "Flux Ratio" 
set arrow from 0,1 to 0.01,1 nohead
plot "Sphere2.tsv" title "" with errorbars
clear
set terminal pngcairo nocrop enhanced font 'Verdana,10' size 1000,500 
set output 'Sphere3.png'
set title "Flux Ratio for cold Hydrogen sphere of radius 0.0015 cm" 
set xlabel "Energy (MeV)" 
set ylabel "Flux Ratio" 
set arrow from 0,1 to 0.01,1 nohead
plot "Sphere3.tsv" title "" with errorbars
clear
set terminal pngcairo nocrop enhanced font 'Verdana,10' size 1000,500 
set output 'Sphere4.png'
set title "Flux Ratio for cold Hydrogen sphere of radius 0.002 cm" 
set xlabel "Energy (MeV)" 
set ylabel "Flux Ratio" 
set arrow from 0,1 to 0.01,1 nohead
plot "Sphere4.tsv" title "" with errorbars
clear
set terminal pngcairo nocrop enhanced font 'Verdana,10' size 1000,500 
set output 'Sphere5.png'
set title "Flux Ratio for cold Hydrogen sphere of radius 0.0025 cm" 
set xlabel "Energy (MeV)" 
set ylabel "Flux Ratio" 
set arrow from 0,1 to 0.01,1 nohead
plot "Sphere5.tsv" title "" with errorbars
clear

