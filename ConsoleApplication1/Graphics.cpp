#include"Graphics.h"

void Errors(BurgersSystem& iburgersSystem, vector<double**>& iresult)
{
	ofstream out7;
	ofstream out9;
	ofstream out11;
	ofstream out12;
	out7.open("Error_U.txt");
	out9.open("Error_V.txt");
	out11.open("U_component_3D.plt");
	out12.open("V_component_3D.plt");

	double x = 0;
	double y = 1;

	double relative_error_u = 0;
	double relative_error_v = 0;
	double abs_error_u = 0;
	double abs_error_v = 0;
	double tmp;

	//Точка в которой максимум погрешности
	double x_U, y_U;
	double x_V, y_V;

	for (int i = 0; i < n + 1; i++) {
		for (int j = 0; j < n + 1; j++) {
			Transform(i, j, x, y);
			out7 << x << " ";
			out7 << y << " ";

			if (i == 0 || i == n || j == 0 || j == n)
				out7 << 0 << endl;
			else {
				tmp = fabs(iburgersSystem.Exact_Solution_Burgers(x, y, time_max)[0] - iresult[0][i][j]);
				out7 << tmp << endl;
				if ((tmp / fabs(iburgersSystem.Exact_Solution_Burgers(x, y, time_max)[0])) * 100 > relative_error_u)
					relative_error_u = 100 * (tmp / fabs(iburgersSystem.Exact_Solution_Burgers(x, y, time_max)[0]));
				if (tmp > abs_error_u) {
					abs_error_u = tmp;
					x_U = x;
					y_U = y;
				}
			}

		}
	}
	x = 0;
	y = 1;
	for (int i = 0; i < n + 1; i++) {
		for (int j = 0; j < n + 1; j++) {
			Transform(i, j, x, y);
			out9 << x << " ";
			out9 << y << " ";

			if (i == 0 || i == n || j == 0 || j == n)
				out9 << 0 << endl;
			else {
				tmp = fabs(iburgersSystem.Exact_Solution_Burgers(x, y, time_max)[1] - iresult[1][i][j]);
				out9 << tmp << endl;
				if (tmp / fabs(iburgersSystem.Exact_Solution_Burgers(x, y, time_max)[1]) * 100 > relative_error_v)
					relative_error_v = 100 * tmp / fabs(iburgersSystem.Exact_Solution_Burgers(x, y, time_max)[1]);
				if (tmp > abs_error_v) {
					abs_error_v = tmp;
					x_V = x;
					y_V = y;
				}
			}

		}
	}

	cout << "Relative Error U : " << relative_error_u << endl;
	cout << "Abs Error U : " << abs_error_u << endl;
	cout << "Point : " << x_U << " " << y_U << endl;

	cout << "Relative Error V : " << relative_error_v << endl;
	cout << "Abs Error V : " << abs_error_v << endl;
	cout << "Point : " << x_V << " " << y_V << endl;


	out11 << "set ticslevel 0" << endl;
	out12 << "set ticslevel 0" << endl;
	out11 << "set dgrid3d 300,300" << endl;
	out12 << "set dgrid3d 300,300 " << endl;

	out11 << "set pm3d at bs" << endl;
	out11 << "set hidden3d" << endl;
	out11 << "set xlabel 'X' " << endl;
	out11 << "set ylabel 'Y' " << endl;

	out12 << "set pm3d at bs" << endl;
	out12 << "set hidden3d" << endl;
	out12 << "set xlabel 'X' " << endl;
	out12 << "set ylabel 'Y' " << endl;

	out11 << "splot 'C:/Users/1/source/repos/ConsoleApplication1/ConsoleApplication1/Error_U.txt' with lines title 'Абсолютная погрешность U' " << endl;
	out12 << "splot 'C:/Users/1/source/repos/ConsoleApplication1/ConsoleApplication1/Error_V.txt' with lines title 'Абсолютная погрешность V' " << endl;

	out7.close();
	out9.close();
	out11.close();
	out12.close();
	system("start U_component_3D.plt");
	system("start V_component_3D.plt");
}

void ExactAndNumericalSolutions(BurgersSystem& iburgersSystem, vector<double**>& iresult)
{

	ofstream out1;
	ofstream out2;
	ofstream out3;
	ofstream out4;
	ofstream out5;
	ofstream out6;

	out1.open("odin.txt");
	out2.open("2.plt");
	out5.open("5.plt");
	out3.open("tri.txt");
	out4.open("V_component.txt");
	out6.open("V_component_Exact_Sollution.txt");

	double x = 0;
	const int test = n / 2;
	//U_component
	for (int j = 0; j < n + 1; j++)
	{
		if (j == 0 || j == n) {
			out3 << x << " " << iburgersSystem.Exact_Solution_Burgers(x, side_of_the_square - test * step_space, time_max)[0] << endl;
		}
		else {
			out3 << x << " " << iresult[0][test][j] << endl; ;
		}
		x = x + step_space;
	}

	x = 0;
	for (int j = 0; j <= 10; j++) {
		out1 << x << " " << iburgersSystem.Exact_Solution_Burgers(x, side_of_the_square - test * step_space, time_max)[0] << endl;
		x += side_of_the_square / 10;
	}

	//V_component
	x = 0;
	for (int j = 0; j < n + 1; j++)
	{
		if (j == 0 || j == n) {
			out4 << x << " " << iburgersSystem.Exact_Solution_Burgers(x, side_of_the_square - test * step_space, time_max)[1] << endl;
		}
		else {
			out4 << x << " " << iresult[1][test][j] << endl; ;
		}
		x = x + step_space;
	}

	x = 0;
	for (int j = 0; j <= 10; j++) {
		out6 << x << " " << iburgersSystem.Exact_Solution_Burgers(x, side_of_the_square - test * step_space, time_max)[1] << endl;
		x = x + side_of_the_square / 10;
	}

	cout << "y = " << side_of_the_square - test * step_space << endl;
	out2 << "set xrange[" << 0 - fabs(0.1 * side_of_the_square) << ":" << side_of_the_square + fabs(0.1 * side_of_the_square) << "]" << endl;
	out5 << "set xrange[" << 0 - fabs(0.1 * side_of_the_square) << ":" << side_of_the_square + fabs(0.1 * side_of_the_square) << "]" << endl;

	out2 << "set xlabel 'x' " << endl;
	out2 << "set ylabel 'Функции' " << endl;
	out5 << "set xlabel 'x' " << endl;
	out5 << "set ylabel 'Функции' " << endl;

	out2 << "plot 'C:/Users/1/source/repos/ConsoleApplication1/ConsoleApplication1/odin.txt' with lines lt rgb 'red' lw 4 title 'точное решение U','C:/Users/1/source/repos/ConsoleApplication1/ConsoleApplication1/tri.txt' with lp lt rgb 'black' dashtype 7 lw 2 pt 7 ps 1 title 'приближенное решение U' " << endl;
	out5 << "plot 'C:/Users/1/source/repos/ConsoleApplication1/ConsoleApplication1/V_component_Exact_Sollution.txt' with lines lt rgb 'red' lw 4 title 'точное решение V','C:/Users/1/source/repos/ConsoleApplication1/ConsoleApplication1/V_component.txt' with lp lt rgb 'black' dashtype 7 lw 2 pt 7 ps 1 title 'приближенное решение V' " << endl;


	out1.close();
	out2.close();
	out3.close();
	out4.close();
	out5.close();
	out6.close();
	system("start 2.plt");
	system("start 5.plt");

}

void BurgersSystem::ExactAndNumericalSolutionsAxisT(const vector<vector<double>> isolution_axis_t)
{
	ofstream out1;
	ofstream out2;
	ofstream out3;
	ofstream out4;
	ofstream out5;
	ofstream out6;

	out1.open("U_axis_t.txt");
	out2.open("U_exact_axis_t.txt");
	out3.open("V_axis_t.txt");
	out4.open("V_exact_axis_t.txt");
	out5.open("pusk_U_axis_t.plt");
	out6.open("pusk_V_axis_t.plt");

	//U_component
	double t = 0;
	for (int j = 0; j <= iterations; j++)
	{
		out1 << t << " " << isolution_axis_t[0][j] << endl;
		t = t + step_time;
	}
	t = 0;
	for (int j = 0; j <= 200; j++) {
		out2 << t << " " << Exact_Solution_Burgers(side_of_the_square / 2, side_of_the_square / 2, t)[0] << endl;
		t += time_max / 200;
	}
	//V_component
	t = 0;
	for (int j = 0; j <= iterations; j++)
	{
		out3 << t << " " << isolution_axis_t[1][j] << endl;
		t = t + step_time;
	}
	t = 0;
	for (int j = 0; j <= 200; j++) {
		out4 << t << " " << Exact_Solution_Burgers(side_of_the_square / 2, side_of_the_square / 2, t)[1] << endl;
		t += time_max / 200;
	}

	out5 << "set xrange[" << 0 - step_time * 2 << ":" << time_max + step_time * 2 << "]" << endl;
	out6 << "set xrange[" << 0 - step_time * 2 << ":" << time_max + step_time * 2 << "]" << endl;

	out5 << "set xlabel 't' " << endl;
	out5 << "set ylabel 'Функции' " << endl;
	out6 << "set xlabel 't' " << endl;
	out6 << "set ylabel 'Функции' " << endl;

	out5 << "plot 'C:/Users/1/source/repos/ConsoleApplication1/ConsoleApplication1/U_exact_axis_t.txt' with lines lt rgb 'blue' lw 4 title 'точное решение U','C:/Users/1/source/repos/ConsoleApplication1/ConsoleApplication1/U_axis_t.txt' with lp lt rgb 'orange' dashtype 7 lw 2 pt 7 ps 1 title 'приближенное решение U' " << endl;
	out6 << "plot 'C:/Users/1/source/repos/ConsoleApplication1/ConsoleApplication1/V_exact_axis_t.txt' with lines lt rgb 'blue' lw 4 title 'точное решение V','C:/Users/1/source/repos/ConsoleApplication1/ConsoleApplication1/V_axis_t.txt' with lp lt rgb 'orange' dashtype 7 lw 2 pt 7 ps 1 title 'приближенное решение V' " << endl;

	out1.close();
	out2.close();
	out3.close();
	out4.close();
	out5.close();
	out6.close();
	system("pusk_U_axis_t.plt");
	system("pusk_V_axis_t.plt");
}

vector<double> Normalization(const vector<double> vec)
{
	vector<double> vector(2);
	vector[0] = vec[0] / sqrt(pow(vec[0], 2) + pow(vec[1], 2));
	vector[1] = vec[1] / sqrt(pow(vec[0], 2) + pow(vec[1], 2));
	return vector;
}

void FieldVelocity(BurgersSystem& iburgersSystem, vector<double**>& iresult) 
{
	ofstream out1;
	ofstream out2;
	ofstream out3;
	ofstream out4;

	out1.open("Exact_field_velocity.txt");
	out2.open("Field_velocity.txt");
	out3.open("fv1_pusk.plt");
	out4.open("fv2_pusk.plt");

	double x, y;
	double delta = side_of_the_square * 0.05;
	vector<double> tmp(2);
	for (int i = 0; i < n + 1; i++) {
		for (int j = 0; j < n + 1; j++) {
			Transform(i, j, x, y);
			out1 << x << " ";
			out1 << y << " ";
			out1 << Normalization(iburgersSystem.Exact_Solution_Burgers(x, y, time_max))[0] * delta << " ";
			out1 << Normalization(iburgersSystem.Exact_Solution_Burgers(x, y, time_max))[1] * delta << endl;

			out2 << x << " ";
			out2 << y << " ";
			tmp[0] = iresult[0][i][j];
			tmp[1] = iresult[1][i][j];
			out2 << Normalization(tmp)[0] * delta << " ";
			out2 << Normalization(tmp)[1] * delta << endl;
		}
	}

	out3 << "set xlabel 'u' " << endl;
	out3 << "set ylabel 'v' " << endl;
	out4 << "set xlabel 'u' " << endl;
	out4 << "set ylabel 'v' " << endl;

	out3 << "set xrange[" << 0 - fabs(0.1 * side_of_the_square) << ":" << side_of_the_square + fabs(0.1 * side_of_the_square) << "]" << endl;
	out4 << "set xrange[" << 0 - fabs(0.1 * side_of_the_square) << ":" << side_of_the_square + fabs(0.1 * side_of_the_square) << "]" << endl;

	out3 << "plot 'C:/Users/1/source/repos/ConsoleApplication1/ConsoleApplication1/Exact_field_velocity.txt' with vectors head filled lt 2 title 'Exact' " << endl;
	out4 << "plot 'C:/Users/1/source/repos/ConsoleApplication1/ConsoleApplication1/Field_velocity.txt' with vectors head filled lt 2 title 'Приближ' " << endl;


	out1.close();
	out2.close();
	out3.close();
	out4.close();
	system("fv1_pusk.plt");
	system("fv2_pusk.plt");
}

void Compare(BurgersSystem& burgersSystem, const vector<double**>& result)
{
	int i, j;

	cout << endl << endl;
	cout << "U component : " << endl;

	Transform2(0.1, 0.1, i, j);
	cout << "Numerical : " << result[0][i][j] << "     " << "Exact : " << burgersSystem.Exact_Solution_Burgers(0.1, 0.1, 2)[0] << endl;

	Transform2(0.5, 0.1, i, j);
	cout << "Numerical : " << result[0][i][j] << "     " << "Exact : " << burgersSystem.Exact_Solution_Burgers(0.5, 0.1, 2)[0] << endl;

	Transform2(0.9, 0.1, i, j);
	cout << "Numerical : " << result[0][i][j] << "     " << "Exact : " << burgersSystem.Exact_Solution_Burgers(0.9, 0.1, 2)[0] << endl;

	Transform2(0.3, 0.3, i, j);
	cout << "Numerical : " << result[0][i][j] << "     " << "Exact : " << burgersSystem.Exact_Solution_Burgers(0.3, 0.3, 2)[0] << endl;

	Transform2(0.7, 0.3, i, j);
	cout << "Numerical : " << result[0][i][j] << "     " << "Exact : " << burgersSystem.Exact_Solution_Burgers(0.7, 0.3, 2)[0] << endl;

	Transform2(0.1, 0.5, i, j);
	cout << "Numerical : " << result[0][i][j] << "     " << "Exact : " << burgersSystem.Exact_Solution_Burgers(0.1, 0.5, 2)[0] << endl;

	Transform2(0.5, 0.5, i, j);
	cout << "Numerical : " << result[0][i][j] << "     " << "Exact : " << burgersSystem.Exact_Solution_Burgers(0.5, 0.5, 2)[0] << endl;

	Transform2(0.9, 0.5, i, j);
	cout << "Numerical : " << result[0][i][j] << "     " << "Exact : " << burgersSystem.Exact_Solution_Burgers(0.9, 0.5, 2)[0] << endl;

	Transform2(0.3, 0.7, i, j);
	cout << "Numerical : " << result[0][i][j] << "     " << "Exact : " << burgersSystem.Exact_Solution_Burgers(0.3, 0.7, 2)[0] << endl;

	Transform2(0.7, 0.7, i, j);
	cout << "Numerical : " << result[0][i][j] << "     " << "Exact : " << burgersSystem.Exact_Solution_Burgers(0.7, 0.7, 2)[0] << endl;

	Transform2(0.1, 0.9, i, j);
	cout << "Numerical : " << result[0][i][j] << "     " << "Exact : " << burgersSystem.Exact_Solution_Burgers(0.1, 0.9, 2)[0] << endl;

	Transform2(0.5, 0.9, i, j);
	cout << "Numerical : " << result[0][i][j] << "     " << "Exact : " << burgersSystem.Exact_Solution_Burgers(0.5, 0.9, 2)[0] << endl;

	Transform2(0.9, 0.9, i, j);
	cout << "Numerical : " << result[0][i][j] << "     " << "Exact : " << burgersSystem.Exact_Solution_Burgers(0.9, 0.9, 2)[0] << endl;



	cout << endl << endl;
	cout << "V component : " << endl;

	Transform2(0.1, 0.1, i, j);
	cout << "Numerical : " << result[1][i][j] << "     " << "Exact : " << burgersSystem.Exact_Solution_Burgers(0.1, 0.1, 2)[1] << endl;

	Transform2(0.5, 0.1, i, j);
	cout << "Numerical : " << result[1][i][j] << "     " << "Exact : " << burgersSystem.Exact_Solution_Burgers(0.5, 0.1, 2)[1] << endl;

	Transform2(0.9, 0.1, i, j);
	cout << "Numerical : " << result[1][i][j] << "     " << "Exact : " << burgersSystem.Exact_Solution_Burgers(0.9, 0.1, 2)[1] << endl;

	Transform2(0.3, 0.3, i, j);
	cout << "Numerical : " << result[1][i][j] << "     " << "Exact : " << burgersSystem.Exact_Solution_Burgers(0.3, 0.3, 2)[1] << endl;

	Transform2(0.7, 0.3, i, j);
	cout << "Numerical : " << result[1][i][j] << "     " << "Exact : " << burgersSystem.Exact_Solution_Burgers(0.7, 0.3, 2)[1] << endl;

	Transform2(0.1, 0.5, i, j);
	cout << "Numerical : " << result[1][i][j] << "     " << "Exact : " << burgersSystem.Exact_Solution_Burgers(0.1, 0.5, 2)[1] << endl;

	Transform2(0.5, 0.5, i, j);
	cout << "Numerical : " << result[1][i][j] << "     " << "Exact : " << burgersSystem.Exact_Solution_Burgers(0.5, 0.5, 2)[1] << endl;

	Transform2(0.9, 0.5, i, j);
	cout << "Numerical : " << result[1][i][j] << "     " << "Exact : " << burgersSystem.Exact_Solution_Burgers(0.9, 0.5, 2)[1] << endl;

	Transform2(0.3, 0.7, i, j);
	cout << "Numerical : " << result[1][i][j] << "     " << "Exact : " << burgersSystem.Exact_Solution_Burgers(0.3, 0.7, 2)[1] << endl;

	Transform2(0.7, 0.7, i, j);
	cout << "Numerical : " << result[1][i][j] << "     " << "Exact : " << burgersSystem.Exact_Solution_Burgers(0.7, 0.7, 2)[1] << endl;

	Transform2(0.1, 0.9, i, j);
	cout << "Numerical : " << result[1][i][j] << "     " << "Exact : " << burgersSystem.Exact_Solution_Burgers(0.1, 0.9, 2)[1] << endl;

	Transform2(0.5, 0.9, i, j);
	cout << "Numerical : " << result[1][i][j] << "     " << "Exact : " << burgersSystem.Exact_Solution_Burgers(0.5, 0.9, 2)[1] << endl;

	Transform2(0.9, 0.9, i, j);
	cout << "Numerical : " << result[1][i][j] << "     " << "Exact : " << burgersSystem.Exact_Solution_Burgers(0.9, 0.9, 2)[1] << endl;
}