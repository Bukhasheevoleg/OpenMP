#pragma once

#include <cmath> 
#include <iostream> 
#include<ctime> 
#include <fstream> 
#include <math.h> 
#include <vector> 

namespace EquationData {
	using namespace std;

	extern const double side_of_the_square;
	extern const int n;
	extern const double number_of_particles;
	extern const double time_max;
	extern const double step_space;
	extern const int iterations;
	extern const double step_time;

	class Equation;
	class BurgersSystem;

	enum Direction {
		Right,
		Up,
		Left,
		Down,
	};

	enum Component {
		U_component,
		V_component,
	};

	void Transform(const int a1, const int b1, double& x, double& y);
	void Transform2(const double x, const double y, int& a, int& b);
	void create(double**& A);
	void Delete(double**& A);
	void zero(double**& A);

	class Equation
	{
	public:

		vector<double**> velocity;

		vector<double**> SolveEquations(const double time);
	/*	double** SolveEquation(const double time, Component component);*/
		Equation(const std::vector<double**> &ivelocity);

		//Чтобы можно было в main рисовать
		vector<double> Exact_Solution(double x, double y, double t);
	private:

		double coff_C(double x, double y);
		double coff_D(double x, double y);
		double coff_S(double x, double y);
		vector<double> Velocity(const int a1, const int b1);
		vector<double> Boundary_Condition(const int a1, const int b1, const double t);
		vector<double> Internal_Condition(const int a1, const int b1, const double t);


		Direction DirectionMovement(const double rand, const int a1, const int b1);
		double Pr_Survive(const int a1, const int b1);
		double Pr_Jump(const int a1, const int b1);
		void Generating_a_node(int n, int& a, int& b);
		double sgn(const double x);
		bool Particle_Inside_The_Area(const int a, const int b);
		void walk(double**& A, int& a, int& b);
	};

	class BurgersSystem
	{
	public:

		std::vector<double**> SolveBurgersSystem();
		std::vector<double> Exact_Solution_Burgers(const double x, const double y, const double t);
	private:

		double coff_D(double x, double y);
		std::vector<double> Initial_Condition_Burgers(const int a, const int b);
		std::vector<double**> Massiv_Initial_Conditions();
		void ExactAndNumericalSolutionsAxisT(const vector<vector<double>> isolution_axis_t);
	};
}
