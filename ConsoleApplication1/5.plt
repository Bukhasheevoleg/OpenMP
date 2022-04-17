set xrange[-0.1:1.1]
set xlabel 'x' 
set ylabel 'Функции' 
plot 'C:/Users/1/source/repos/ConsoleApplication1/ConsoleApplication1/V_component_Exact_Sollution.txt' with lines lt rgb 'red' lw 4 title 'точное решение V','C:/Users/1/source/repos/ConsoleApplication1/ConsoleApplication1/V_component.txt' with lp lt rgb 'black' dashtype 7 lw 2 pt 7 ps 1 title 'приближенное решение V' 
