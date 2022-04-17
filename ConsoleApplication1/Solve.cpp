
#include"Input.h"
#include"Graphics.h"
using namespace EquationData;

int main()
{
	double** u_function = new double* [n + 1];
	double** v_function = new double* [n + 1];
	create(u_function);
	create(v_function);
	std::vector<double**> result = { u_function, v_function };

	BurgersSystem burgersSystem;


	double start = clock();
	result = burgersSystem.SolveBurgersSystem();
	double end = clock();
	double time = (end - start) / CLOCKS_PER_SEC;
	cout << "time:" << time << endl;


	Compare(burgersSystem, result);
	FieldVelocity(burgersSystem, result);
	ExactAndNumericalSolutions(burgersSystem, result);
	Errors(burgersSystem, result);
}


