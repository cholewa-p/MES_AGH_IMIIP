#include <iostream>
#include "Matrix.h"
#include "Simulation.h"

using namespace std;

void simulation(double** H, double** C, double* P, double* t_start, int time, int step, int n)
{
	//cout << "Step 1\n";
	//for (int i = 0; i < n; i++)
	//	cout << t_start[i] << " ";
	//cout << endl;
	cout << "Time[s]\tMinTemp[*C]\tMaxTemp[*C]\n";
	for (int i = 0; i < time / step; i++)
	{
		double* temperature = calculate_temperature(H, C, P, t_start, step, n);
		delete[] t_start;
		/*
		
		cout << "Step " << i + 1 << "\n";
		for (int i = 0; i < n; i++)
			cout << temperature[i] << " ";
			cout << endl;
		*/
		double min = min_temperature(temperature, n);
		double max = max_temperature(temperature, n);
		cout << (time / (time / step)) * (i + 1) << "\t" << min << "\t\t" << max << "\n";
		t_start = temperature;
		temperature = nullptr;
	}
}

double* calculate_temperature(double** H, double** C, double* P, double* t_start, int step, int n)
{
	double** newC = devide_matrix_by_number(C, step, n);
	double** newMatrixHC = sum_matrix(H, newC, n);
	for (int i = 0; i < n; i++)
		delete[] newC[i];
	delete[] newC;

	//for (int i = 0; i < n; i++)
	//{
	//	for (int j = 0; j < n; j++)
	//	{
	//		cout << newMatrixHC[i][j] << " ";
	//	}
	//	cout << endl;
	//}
	//cout << endl;

	newC = devide_matrix_by_number(C, step, n);
	double* newVector = multiplication_matrix_by_vector(newC, t_start, n);
	for (int i = 0; i < n; i++)
		delete[] newC[i];
	delete[] newC;
	/*
	for (int i = 0; i < n; i++)
	{
		cout << P[i] << " ";
	}
	cout << endl;
	*/
	double* newVector2 = sum_vector(newVector, P, n);

	//for (int i = 0; i < n; i++)
	//{
	//	cout << newVector2[i] << " ";
	//}
	//cout << endl;

	double* new_temperature = Gauss_elimination(newMatrixHC, newVector2, n);

	for (int i = 0; i < n; i++)
		delete[] newMatrixHC[i];
	delete[] newMatrixHC;

	delete[] newVector;
	delete[] newVector2;

	return new_temperature;
}

double max_temperature(double* temperature, int n)
{
	double max = temperature[0];

	for (int i = 1; i < n; i++)
	{
		if (max < temperature[i])
		{
			max = temperature[i];
		}
	}
	return max;
}

double min_temperature(double* temperature, int n)
{
	double min = temperature[0];

	for (int i = 1; i < n; i++)
	{
		if (min > temperature[i])
		{
			min = temperature[i];
		}
	}

	return min;
}
