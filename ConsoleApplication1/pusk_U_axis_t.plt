set xrange[-0.5:5.5]
set xlabel 't' 
set ylabel '�������' 
plot 'C:/Users/1/source/repos/ConsoleApplication1/ConsoleApplication1/U_exact_axis_t.txt' with lines lt rgb 'blue' lw 4 title '������ ������� U','C:/Users/1/source/repos/ConsoleApplication1/ConsoleApplication1/U_axis_t.txt' with lp lt rgb 'orange' dashtype 7 lw 2 pt 7 ps 1 title '������������ ������� U' 
