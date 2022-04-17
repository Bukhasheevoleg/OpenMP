#include"Input.h"

namespace EquationData {

	const double side_of_the_square = 1;
	const int n = 10;
	const double number_of_particles = 1.E5;
	const double time_max = 5;//5
	const double step_space = side_of_the_square / n;
	const int iterations = 20;//20
	const double step_time = time_max / iterations;
	

	void Transform(const int a1, const int b1, double& x, double& y) 
	{
		x = step_space * b1;
		y = side_of_the_square - step_space * a1;
	}

	void Transform2(const double x, const double y, int& a, int& b)
	{
		if (x > side_of_the_square || y > side_of_the_square)
			throw;
		a = (side_of_the_square - y) / step_space;
		b = x / step_space;
	}

	bool BelongsSegment(const double point, const double start_segment, const double end_segment)
	{
		if (point >= start_segment && point <= end_segment)
			return true;
		else
			return false;
	}

	void create(double**& A) 
	{
		for (int i = 0; i < n + 1; i++) {
			A[i] = new double[n + 1];
		}
	}

	void Delete(double**& A)
	{
		for (int i = 0; i < n + 1; i++) {
			delete[] A[i];
		}
		delete[] A;
	}

	void zero(double**& A)
	{
		for (int i = 0; i < n + 1; i++) {
			for (int j = 0; j < n + 1; j++) {
				A[i][j] = 0;
			}
		}

	}

	Equation::Equation(const std::vector<double**>& ivelocity)
	{
		this->velocity = ivelocity;
	}

	double Equation::coff_C(double x, double y) 
	{
		double c;
		c = 0;
		return c;
	}

	double Equation::coff_D(double x, double y)
	{
		double c;
		c = 100;
		return c;
	}

	double Equation::coff_S(double x, double y) {
		double c;
		c = 100;
		return c;
	}

	vector<double> Equation::Velocity(const int a1, const int b1) 
	{
		vector<double> vec(2);

		vec[0] = velocity[0][a1][b1];
		vec[1] = velocity[1][a1][b1];

		double x, y;
		Transform(a1, b1, x, y);
		if (vec[0] < (2 * coff_D(x, y)) / step_space && vec[1] < (2 * coff_D(x, y)) / step_space)
			return vec;
		else
			throw "Warning";
	}

	vector<double> Equation::Exact_Solution(double x, double y, double t) 
	{
		double c, tmp;

		vector<double> vec(2);

		//Problem 1
	/*	tmp = (-4 * x + 4 * y - t) / (32 * coff_D(x, y));
		vec[0] = 3.0 / 4 - 1.0 / ((1 + exp(tmp)) * 4);
		vec[1] = 3.0 / 4 + 1.0 / ((1 + exp(tmp)) * 4);*/

		//Problem 2
		vec[0] = (x + y - 2 * x * t) / (1 - 2 * pow(t, 2));
		vec[1] = (x - y - 2 * y * t) / (1 - 2 * pow(t, 2));


		return vec;
	}

	vector<double> Equation::Boundary_Condition(const int a1, const int b1, double t) 
	{
		double x, y, tmp;
		vector<double> vec(2);
		Transform(a1, b1, x, y);


		////Robbin boundary conditions
		//double tmp1, tmp2;
		//if (x == 1) {
		//	tmp1 = (1 - 2 * t) / (1 - 2 * pow(t, 2));
		//	tmp2 = 1 / (1 - 2 * pow(t, 2));
		//}
		//else
		//	if (x == 0) {
		//		tmp1 = -(1 - 2 * t) / (1 - 2 * pow(t, 2));
		//		tmp2 = -1 / (1 - 2 * pow(t, 2));
		//	}
		//	else 
		//		if (y == 1) {
		//			tmp1 = 1 / (1 - 2 * pow(t, 2));
		//			tmp2 = (-1 - 2 * t) / (1 - 2 * pow(t, 2));
		//		}
		//		else
		//			if (y == 0) {
		//				tmp1 = - 1 / (1 - 2 * pow(t, 2));
		//				tmp2 = - (-1 - 2 * t) / (1 - 2 * pow(t, 2));
		//			}

		//vec[0] = coff_D(x, y) * tmp1 + coff_S(x, y) * Exact_Solution(x, y, t)[0];
		//vec[1] = coff_D(x, y) * tmp2 + coff_S(x, y) * Exact_Solution(x, y, t)[1];


		vec = Exact_Solution(x, y, t);
		return vec;
	}

	vector<double> Equation::Internal_Condition(const int a1, const int b1, const double t)
	{
		double f, x, y;
		Transform(a1, b1, x, y);
		vector<double> vec(2);

		vec[0] = 0;
		vec[1] = 0;

		return vec;
	}

	vector<double**> Equation::SolveEquations(const double time) 
	{
		//grid
		double** A = new double* [n + 1];
		create(A);
		zero(A);
		//веса узлов
		double** B = new double* [n + 1];
		create(B);
		zero(B);

		//учет траекторий посетивших узел
		double** number_of_trajectories = new double* [n + 1];
		create(number_of_trajectories);
		zero(number_of_trajectories);
		//счетчик дл€ нулевых граничных условий
		double** C = new double* [n + 1];
		create(C);
		zero(C);


		//ѕараллелю
		double** B2 = new double* [n + 1];
		create(B2);
		zero(B2);
		double** C2 = new double* [n + 1];
		create(C2);
		zero(C2);


		int a, b;//строка и столбец узла
		double d;//f(Po,i)
		double Q;

		double internal_V;
		double boundary_V;

		for (int k = 0; k < number_of_particles; k++) {

			//обнуление траектории
			zero(A);
			Generating_a_node(n, a, b);

			//нулевые граничные услови€
			d = Internal_Condition(a, b, time)[0];
			internal_V = Internal_Condition(a, b, time)[1];

			/*if (component == U_component)
				d = Internal_Condition(a, b, time)[0];
			else
				d = Internal_Condition(a, b, time)[1];*/


			A[a][b] = A[a][b] + 1;

			walk(A, a, b);

			//0ые граничные услови€
			for (int i = 0; i < n + 1; i++) {
				for (int j = 0; j < n + 1; j++) {
					C[i][j] = C[i][j] + A[i][j] * d;
				}
			}

			for (int i = 0; i < n + 1; i++) {
				for (int j = 0; j < n + 1; j++) {
					C2[i][j] = C2[i][j] + A[i][j] * internal_V;
				}
			}

			//ненулевые граничные услови€
			if (a == n || b == n || a == 0 || b == 0) {


				Q = Boundary_Condition(a, b, time)[0];
				boundary_V = Boundary_Condition(a, b, time)[1];


				for (int i = 0; i < n + 1; i++) {
					for (int j = 0; j < n + 1; j++) {
						B[i][j] = B[i][j] + sgn(A[i][j]) * Q;
					}
				}

				for (int i = 0; i < n + 1; i++) {
					for (int j = 0; j < n + 1; j++) {
						B2[i][j] = B2[i][j] + sgn(A[i][j]) * boundary_V;
					}
				}

			}


			//учет траекторий посетивших узел
			for (int i = 0; i < n + 1; i++) {
				for (int j = 0; j < n + 1; j++) {
					number_of_trajectories[i][j] += sgn(A[i][j]);
				}
			}

		}


		double** u_function = new double* [n + 1];
		double** v_function = new double* [n + 1];
		create(u_function);
		create(v_function);
		vector<double**> result = { u_function, v_function };

		double G = ((n - 1) * (n - 1) * pow(step_space, 2)) / (number_of_particles * 4);

		for (int i = 0; i < n + 1; i++) {
			for (int j = 0; j < n + 1; j++) {
				if (i == 0 || j == 0 || i == n || j == n) {
					u_function[i][j] = Boundary_Condition(i, j, time)[0];
					v_function[i][j] = Boundary_Condition(i, j, time)[1];
				}
				else {
					u_function[i][j] = B[i][j] / number_of_trajectories[i][j] + C[i][j] * G;
					v_function[i][j] = B2[i][j] / number_of_trajectories[i][j] + C2[i][j] * G;
				}
			}
		}

		Delete(A);
		Delete(B);
		Delete(B2);
		Delete(C);
		Delete(C2);
		Delete(number_of_trajectories);
		return result;
	}

	Direction Equation::DirectionMovement(const double rand, const int a1, const int b1)
	{
		Direction direction;
		double x, y;
		Transform(a1, b1, x, y);
		vector<double> velocity = Velocity(a1, b1);
		double pr_right = (2 * coff_D(x, y) - step_space * velocity[0]) / (8 * coff_D(x, y));
		double pr_left = (2 * coff_D(x, y) + step_space * velocity[0]) / (8 * coff_D(x, y));
		double pr_up = (2 * coff_D(x, y) - step_space * velocity[1]) / (8 * coff_D(x, y));
		double pr_down = (2 * coff_D(x, y) + step_space * velocity[1]) / (8 * coff_D(x, y));
		vector<double>  probabilities = { pr_right, pr_left, pr_up, pr_down };

		if (BelongsSegment(rand, 0, pr_right))
			return direction = Right;
		if (BelongsSegment(rand, pr_right, pr_right + pr_left))
			return direction = Left;
		if (BelongsSegment(rand, pr_right + pr_left, pr_right + pr_left + pr_up))
			return direction = Up;
		if (BelongsSegment(rand, pr_right + pr_left + pr_up, 1))
			return direction = Down;
	}

	double Equation::Pr_Survive(const int a1, const int b1) 
	{
		double x, y;
		Transform(a1, b1, x, y);
		return  4.0 / (4 + coff_C(x, y) * pow(step_space, 2) / coff_D(x, y));
	}

	/*double Equation::Pr_Jump(const int a1, const int b1) 
	{
		double x, y;
		Transform(a1, b1, x, y);
		return coff_D(x, y) / (coff_D(x, y) + coff_S(x, y) * step_space);
	}*/

	void Equation::Generating_a_node(int n, int& a, int& b) 
	{

		a = rand() % (n - 1) + 1;
		b = rand() % (n - 1) + 1;
	}

	double Equation::sgn(const double x) 
	{
		if (x > 0)return 1;
		if (x == 0)return 0;
	}

	bool Equation::Particle_Inside_The_Area(const int a, const int b) 
	{
		if (a > 0 && a < n && b > 0 && b < n)
			return true;
		else
			return false;
	}

	void Equation::walk(double**& A, int& a, int& b) 
	{
		double direction_rand;
		Direction direction;
		int j = 1;
		double WEIGHT = Pr_Survive(a, b);
		double pr, pr_jump;

		while (a != 0 && b != 0 && a != n && b != n) {

			pr = (double)(rand()) / RAND_MAX;
			if (pr <= WEIGHT) {
				direction_rand = (double)(rand()) / RAND_MAX;
				direction = DirectionMovement(direction_rand, a, b);

				switch (direction) {

				case Down: {
					if (Particle_Inside_The_Area(a + 1, b))
					{
						A[a + 1][b] += 1;
						a++;
					}
					else {
						a++;
						break;
					}
				}
								 break;

				case Up: {
					if (Particle_Inside_The_Area(a - 1, b))
					{
						A[a - 1][b] += 1;
						a--;
					}
					else {
						a--;
						break;
					}
				}
							 break;

				case Right: {
					if (Particle_Inside_The_Area(a, b + 1))
					{
						A[a][b + 1] += 1;
						b++;
					}
					else {
						b++;
						break;
					}
				}
									break;

				case Left: {
					if (Particle_Inside_The_Area(a, b - 1))
					{
						A[a][b - 1] += 1;
						b--;
					}
					else {
						b--;
						break;
					}
				}
								 break;
				}

				////ѕроверка на отпрыгивание назад
				//if (!Particle_Inside_The_Area(a, b)) {
				//	pr_jump = (double)(rand()) / RAND_MAX;
				//	if (pr_jump <= Pr_Jump(a, b)) {

				//		switch (direction) {
				//		case Down: {
				//			a--;
				//			break;
				//		}
				//		case Up: {
				//			a++;
				//			break;
				//		}
				//		case Right: {
				//			b--;
				//			break;
				//		}
				//		case Left: {
				//			b++;
				//			break;
				//		}
				//		default:
				//			break;
				//		}

				//		A[a][b] = A[a][b] + 1;
				//	}
				//}


				WEIGHT = Pr_Survive(a, b);

			}
			else
				break;

		}
	}

	double BurgersSystem::coff_D(double x, double y) 
	{
		double c;
		c = 100;
		return c;
	}

	vector<double> BurgersSystem::Exact_Solution_Burgers(const double x, const double y, const double t)
	{
		vector<double> vec(2);
		//Problem 1
		/*double tmp = (-4 * x + 4 * y - t) / (32 * coff_D(x, y));
		vec[0] = 3.0 / 4 - 1.0 / ((1 + exp(tmp)) * 4);
		vec[1] = 3.0 / 4 + 1.0 / ((1 + exp(tmp)) * 4);*/
		


		//Problem 2
		vec[0] = (x + y - 2 * x * t) / (1 - 2 * pow(t, 2));
		vec[1] = (x - y - 2 * y * t) / (1 - 2 * pow(t, 2));
		return vec;
	}

	vector<double> BurgersSystem::Initial_Condition_Burgers(const int a, const int b)
	{
		vector<double> vec(2);
		double x, y;
		Transform(a, b, x, y);
		//Problem 1
		/*double tmp = (-4 * x + 4 * y) / (32 * coff_D(x, y));
		vec[0] = 3.0 / 4 - 1.0 / ((1 + exp(tmp)) * 4);
		vec[1] = 3.0 / 4 + 1.0 / ((1 + exp(tmp)) * 4);*/
		
		//Problem 2
		vec[0] = x + y;
		vec[1] = x - y;
		return vec;
	}

	vector<double**> BurgersSystem::Massiv_Initial_Conditions()
	{
		std::vector<double**> result(2);
		result[0] = new double* [n+1];
		result[1] = new double* [n+1];
		create(result[0]);
		create(result[1]);

			for (int i = 0; i < n+1; i++) {
				for (int j = 0; j < n+1; j++) {
					result[0][i][j] = Initial_Condition_Burgers(i, j)[0];
					result[1][i][j] = Initial_Condition_Burgers(i, j)[1];
				}
			}

		return result;
	}

	vector<double**> BurgersSystem::SolveBurgersSystem()
	{

		vector<double**> result;
		//for first step
		vector<double**> field_velocity = Massiv_Initial_Conditions();
		vector<vector<double>> solution_axis_t(2, vector<double>(iterations + 1));
		solution_axis_t[0][0] = field_velocity[0][n / 2][n / 2];
		solution_axis_t[1][0] = field_velocity[1][n / 2][n / 2];

		for (int iter = 1; iter <= iterations; iter++) {
			Equation equations(field_velocity);

			if (iter == iterations) {
				result = equations.SolveEquations(iter * step_time);
				solution_axis_t[0][iterations] = field_velocity[0][n / 2][n / 2];
				solution_axis_t[1][iterations] = field_velocity[1][n / 2][n / 2];
				ExactAndNumericalSolutionsAxisT(solution_axis_t);
				return result;
			}

			field_velocity = equations.SolveEquations(iter * step_time);

			solution_axis_t[0][iter] = field_velocity[0][n / 2][n / 2];
			solution_axis_t[1][iter] = field_velocity[1][n / 2][n / 2];
			//Debug
			cout << "step time:" << iter * step_time << endl << endl;
		}

	}

}

