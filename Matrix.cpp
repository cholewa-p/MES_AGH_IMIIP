#include "Element.h"
#include "Grid.h"
#include "Node.h"
#include "Elem4.h"
#include "Matrix.h"
#include "Jakobian.h"
#include <iostream>
#include <cmath>
using namespace std;
#define print 0
void calculate_matrix_H(Grid& grid, Elem4 element)
{
	// Wspó³czynnik przewodzenia
	double countivity = grid.Conductivity;
	//Przechodzenie po ka¿dym elemencie
	for (int i = 0; i < grid.get_amount_elements(); i++)
	{
		double H[4][4] = { {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0} };
#if print==1
		cout << "Element " << i + 1 << endl;
#endif
		//Przechodzenie po ka¿dym wierzcho³ku
		for (int j = 0; j < element.schemat * element.schemat; j++) //Dodac wiêcej punktów
		{
#if print==1
			cout << "Point " << j + 1 << endl;
#endif

			double I[2][2] = { {0,0},{0,0} };
			double Iinv[2][2] = { {0,0},{0,0} };

			jakobian(i, j, I, Iinv, element, grid);

#if print==1
			for (int l = 0; l < 2; l++)
			{
				for (int m = 0; m < 2; m++)
					cout << I[l][m] << " ";
				cout << endl;
			}
			cout << endl;
#endif
			if (element.schemat == 2)
			{
				double dNdx[4][4] = { {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0} };
				double dNdy[4][4] = { {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0} };

				//Liczenie macierzy powierzchni po punktach ca³kowania
				for (int l = 0; l < 4; l++)
					for (int m = 0; m < 4; m++)
					{
						dNdx[l][m] = Iinv[0][0] * element.ksi[l][m] + Iinv[0][1] * element.eta[l][m];
						dNdy[l][m] = Iinv[1][0] * element.ksi[l][m] + Iinv[1][1] * element.eta[l][m];
					}

				double Hx[4][4] = { {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0} };
				double Hy[4][4] = { {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0} };
				double H1[4][4] = { {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0} };

				for (int k = 0; k < 4; k++)
				{
					for (int m = 0; m < 4; m++)
					{
						Hx[k][m] = dNdx[j][k] * dNdx[j][m];
						Hy[k][m] = dNdy[j][k] * dNdy[j][m];
						double detI = I[0][0] * I[1][1] - I[1][0] * I[0][1];

						H1[k][m] = countivity * (Hx[k][m] + Hy[k][m]) * detI;

						H[k][m] += H1[k][m];
#if print==1
						cout << H1[k][m] << " ";
#endif
					}
#if print==1
					cout << endl;
#endif

				}
#if print==1
				cout << endl;
#endif

			}
			else if (element.schemat == 3)
			{
				double dNdx[9][4] = { {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0} };
				double dNdy[9][4] = { {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0} };

				for (int l = 0; l < 9; l++)
					for (int m = 0; m < 4; m++)
					{
						dNdx[l][m] = Iinv[0][0] * element.ksi[l][m] + Iinv[0][1] * element.eta[l][m];
						dNdy[l][m] = Iinv[1][0] * element.ksi[l][m] + Iinv[1][1] * element.eta[l][m];
					}

				double Hx[4][4] = { {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0} };
				double Hy[4][4] = { {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0} };
				double H1[4][4] = { {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0} };

				for (int k = 0; k < 4; k++)
				{
					for (int m = 0; m < 4; m++)
					{
						Hx[k][m] = dNdx[j][k] * dNdx[j][m];
						Hy[k][m] = dNdy[j][k] * dNdy[j][m];
						double detI = I[0][0] * I[1][1] - I[1][0] * I[0][1];

						H1[k][m] = countivity * (Hx[k][m] + Hy[k][m]) * detI * element.weight[j / 3] * element.weight[j % 3];

						H[k][m] += H1[k][m];
#if print==1
						cout << H1[k][m] << " ";
#endif
					}
#if print==1
					cout << endl;
#endif
				}
#if print==1
				cout << endl;
#endif
			}
		}
#if print==1
		for (int l = 0; l < 4; l++)
		{
			for (int m = 0; m < 4; m++)
				cout << H[l][m] << " ";
			cout << endl;
		}
		cout << endl;
#endif
		//save matrix H to the grid
		//grid.save_H_Matrix(H, i);
		for (int j = 0; j < 4; j++)
			for (int k = 0; k < 4; k++)
				grid.elements[i].H[j][k] = H[j][k];

	}
}
#define print 0
void calculate_matrix_Hbc(Grid& grid, Elem4 element)
{
	double alfa = grid.Alfa;
	double endTemperature = grid.Tot;
	double tab[2][4] = {};
	int point;
	double points[2]= { -1 / sqrt(3), 1 / sqrt(3) };
	double points2[3]= { -sqrt(3.0 / 5.0), 0 ,sqrt(3.0 / 5.0) };
	double weights2[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
	for (int i = 0; i < grid.get_amount_elements(); i++)
	{



#ifdef print==1
		//cout << "Element " << i + 1 << endl;
#endif 
		if (grid.get_node(grid.get_element(i).get_id1() - 1).get_bc() == 1)
		{
			if (grid.get_node(grid.get_element(i).get_id2() - 1).get_bc() == 1)
			{
				double Hbc[4][4] = { {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0} };
				double P[4] = { 0,0,0,0 };

				for (int k = 0; k < 4; k++)
				{
					//double det = (grid.get_width() / (grid.get_n_width() - 1)) / 2;

					double det = calculate_determinant(
						grid.get_node(grid.get_element(i).get_id1() - 1).get_x(),
						grid.get_node(grid.get_element(i).get_id1() - 1).get_y(),
						grid.get_node(grid.get_element(i).get_id2() - 1).get_x(), 
						grid.get_node(grid.get_element(i).get_id2() - 1).get_y()) / 2;
					if (element.schemat == 2)
					{
						//calculate for which elements there is 0 in matrix
						for (int i = 0; i < 2; i++)
						{
							point = -1;
							tab[i][0] = 0.25 * (1.0 - points[i]) * (1.0 - point);
							tab[i][1] = 0.25 * (1.0 + points[i]) * (1.0 - point);
							tab[i][2] = 0.25 * (1.0 + points[i]) * (1.0 + point);
							tab[i][3] = 0.25 * (1.0 - points[i]) * (1.0 + point);
						}
						P[k] += alfa * (tab[0][k] * endTemperature + tab[1][k] * endTemperature) * det;
						for (int l = 0; l < 4; l++)
						{
							Hbc[k][l] = alfa * (tab[0][k] * tab[0][l] + tab[1][k] * tab[1][l]) * det;
						}
					}
					else if (element.schemat == 3)
					{
						for (int i = 0; i < 3; i++)
						{
							point = -1;
							tab[i][0] = 0.25 * (1.0 - points[i]) * (1.0 - point);
							tab[i][1] = 0.25 * (1.0 + points[i]) * (1.0 - point);
							tab[i][2] = 0.25 * (1.0 + points[i]) * (1.0 + point);
							tab[i][3] = 0.25 * (1.0 - points[i]) * (1.0 + point);
						}
						P[k] += alfa * (weights2[0] * tab[0][k] * endTemperature + weights2[1] * tab[1][k] * endTemperature + weights2[2] * tab[2][k] * endTemperature) * det;

						for (int l = 0; l < 4; l++)
						{
							Hbc[k][l] = alfa * (tab[0][k] * tab[0][l] * weights2[0] + tab[1][k] * tab[1][l] * weights2[1] + tab[2][k] * tab[2][l] * weights2[2]) * det;
						}
					}
				}
				
				
				//add values to Hbc matrix & P vector from the grid
				for (int k = 0; k < 4; k++)
				{
				grid.elements[i].P[k] += P[k];
					for (int l = 0; l < 4; l++)
					{
						grid.elements[i].Hbc[k][l] += Hbc[k][l];
					}
				}
				//grid.add_boundary_Condition(i, Hbc, P);
				
#ifdef print==1
				/*
				
				for (int z = 0; z < 4; z++)
				{
					for (int v = 0; v < 4; v++)
						cout << Hbc[z][v] << " ";
					cout << endl;
				}
				cout << endl;
				*/
#endif
			}

			if (grid.get_node(grid.get_element(i).get_id4() - 1).get_bc() == 1)
			{
				double Hbc[4][4] = { {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0} };
				double P[4] = { 0,0,0,0 };
				for (int k = 0; k < 4; k++)
				{
					//double det = (grid.get_height() / (grid.get_n_height() - 1)) / 2;

					double det = calculate_determinant(grid.get_node(grid.get_element(i).get_id1() - 1).get_x(), grid.get_node(grid.get_element(i).get_id1() - 1).get_y(), grid.get_node(grid.get_element(i).get_id4() - 1).get_x(), grid.get_node(grid.get_element(i).get_id4() - 1).get_y()) / 2;

					if (element.schemat == 2)
					{
						for (int i = 0; i < 2; i++)
						{
							point = -1;
							tab[i][0] = 0.25 * (1.0 - point) * (1.0 - points[i]);
							tab[i][1] = 0.25 * (1.0 + point) * (1.0 - points[i]);
							tab[i][2] = 0.25 * (1.0 + point) * (1.0 + points[i]);
							tab[i][3] = 0.25 * (1.0 - point) * (1.0 + points[i]);
						}
						P[k] += alfa * (tab[0][k] * endTemperature + tab[1][k] * endTemperature) * det;
						for (int l = 0; l < 4; l++)
						{
							Hbc[k][l] = alfa * (tab[0][k] * tab[0][l] + tab[1][k] * tab[1][l]) * det;
						}
					}
					else if (element.schemat == 3)
					{
						for (int i = 0; i < 3; i++)
						{
							point = -1;
							tab[i][0] = 0.25 * (1.0 - point) * (1.0 - points[i]);
							tab[i][1] = 0.25 * (1.0 + point) * (1.0 - points[i]);
							tab[i][2] = 0.25 * (1.0 + point) * (1.0 + points[i]);
							tab[i][3] = 0.25 * (1.0 - point) * (1.0 + points[i]);
						}
						P[k] += alfa * (weights2[0] * tab[0][k] * endTemperature + weights2[1] * tab[1][k] * endTemperature + weights2[2] * tab[2][k] * endTemperature) * det;
						for (int l = 0; l < 4; l++)
						{
							Hbc[k][l] = alfa * (tab[0][k] * tab[0][l] * weights2[0] + tab[1][k] * tab[1][l] * weights2[1] + tab[2][k] * tab[2][l] * weights2[2]) * det;
						}
					}
				}
				for (int k = 0; k < 4; k++)
				{
					grid.elements[i].P[k] += P[k];
					for (int l = 0; l < 4; l++)
					{
						grid.elements[i].Hbc[k][l] += Hbc[k][l];
					}
				}
				//grid.add_boundary_Condition(i, Hbc, P);
#ifdef print==1	
				/*
				
				for (int z = 0; z < 4; z++)
				{
					for (int v = 0; v < 4; v++)
						cout << Hbc[z][v] << " ";
					cout << endl;
				}
				cout << endl;
				*/
#endif
			}
		}

		if (grid.get_node(grid.get_element(i).get_id3() - 1).get_bc() == 1)
		{
			if (grid.get_node(grid.get_element(i).get_id2() - 1).get_bc() == 1)
			{
				double Hbc[4][4] = { {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0} };
				double P[4] = { 0,0,0,0 };
				for (int k = 0; k < 4; k++)
				{
					//double det = (grid.get_height() / (grid.get_n_height() - 1)) / 2;

					double det = calculate_determinant(grid.get_node(grid.get_element(i).get_id3() - 1).get_x(), grid.get_node(grid.get_element(i).get_id3() - 1).get_y(), grid.get_node(grid.get_element(i).get_id2() - 1).get_x(), grid.get_node(grid.get_element(i).get_id2() - 1).get_y()) / 2;
					if (element.schemat == 2)
					{
						for (int i = 0; i < 2; i++)
						{
							point = 1;
							tab[i][0] = 0.25 * (1.0 - point) * (1.0 - points[i]);
							tab[i][1] = 0.25 * (1.0 + point) * (1.0 - points[i]);
							tab[i][2] = 0.25 * (1.0 + point) * (1.0 + points[i]);
							tab[i][3] = 0.25 * (1.0 - point) * (1.0 + points[i]);
						}
						P[k] += alfa * (tab[0][k] * endTemperature + tab[1][k] * endTemperature) * det;
						for (int l = 0; l < 4; l++)
						{
							Hbc[k][l] = alfa * (tab[0][k] * tab[0][l] + tab[1][k] * tab[1][l]) * det;
						}
					}
					else if (element.schemat == 3)
					{
						for (int i = 0; i < 3; i++)
						{
							point = 1;
							tab[i][0] = 0.25 * (1.0 - point) * (1.0 - points[i]);
							tab[i][1] = 0.25 * (1.0 + point) * (1.0 - points[i]);
							tab[i][2] = 0.25 * (1.0 + point) * (1.0 + points[i]);
							tab[i][3] = 0.25 * (1.0 - point) * (1.0 + points[i]);
						}
						P[k] += alfa * (weights2[0] * tab[0][k] * endTemperature + weights2[1] * tab[1][k] * endTemperature + weights2[2] * tab[2][k] * endTemperature) * det;
						for (int l = 0; l < 4; l++)
						{
							Hbc[k][l] = alfa * (tab[0][k] * tab[0][l] * weights2[0] + tab[1][k] * tab[1][l] * weights2[1] + tab[2][k] * tab[2][l] * weights2[2]) * det;
						}
					}
				}
				for (int k = 0; k < 4; k++)
				{
					grid.elements[i].P[k] += P[k];
					for (int l = 0; l < 4; l++)
					{
						grid.elements[i].Hbc[k][l] += Hbc[k][l];
					}
				}
				//grid.add_boundary_Condition(i, Hbc, P);
				
#ifdef print==1
				/*
				
				for (int z = 0; z < 4; z++)
				{
					for (int v = 0; v < 4; v++)
						cout << Hbc[z][v] << " ";
					cout << endl;
				}
				cout << endl;
				*/
#endif	
			}

			if (grid.get_node(grid.get_element(i).get_id4() - 1).get_bc() == 1)
			{
				double Hbc[4][4] = { {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0} };
				double P[4] = { 0,0,0,0 };
				for (int k = 0; k < 4; k++)
				{
					//double det = (grid.get_width() / (grid.get_n_width() - 1)) / 2;

					double det = calculate_determinant(grid.get_node(grid.get_element(i).get_id3() - 1).get_x(), grid.get_node(grid.get_element(i).get_id3() - 1).get_y(), grid.get_node(grid.get_element(i).get_id4() - 1).get_x(), grid.get_node(grid.get_element(i).get_id4() - 1).get_y()) / 2;

					if (element.schemat == 2)
					{
						for (int i = 0; i < 2; i++)
						{
							point = 1;
							tab[i][0] = 0.25 * (1.0 - points[i]) * (1.0 - point);
							tab[i][1] = 0.25 * (1.0 + points[i]) * (1.0 - point);
							tab[i][2] = 0.25 * (1.0 + points[i]) * (1.0 + point);
							tab[i][3] = 0.25 * (1.0 - points[i]) * (1.0 + point);
						}
						P[k] += alfa * (tab[0][k] * endTemperature + tab[1][k] * endTemperature) * det;
						for (int l = 0; l < 4; l++)
						{
							Hbc[k][l] = alfa * (tab[0][k] * tab[0][l] + tab[1][k] * tab[1][l]) * det;
						}
					}
					else if (element.schemat == 3)
					{
						for (int i = 0; i < 3; i++)
						{
							point = 1;
							tab[i][0] = 0.25 * (1.0 - points[i]) * (1.0 - point);
							tab[i][1] = 0.25 * (1.0 + points[i]) * (1.0 - point);
							tab[i][2] = 0.25 * (1.0 + points[i]) * (1.0 + point);
							tab[i][3] = 0.25 * (1.0 - points[i]) * (1.0 + point);
						}
						P[k] += alfa * (weights2[0] * tab[0][k] * endTemperature + weights2[1] * tab[1][k] * endTemperature + weights2[2] * tab[2][k] * endTemperature) * det;
						for (int l = 0; l < 4; l++)
						{
							Hbc[k][l] = alfa * (tab[0][k] * tab[0][l] * weights2[0] + tab[1][k] * tab[1][l] * weights2[1] + tab[2][k] * tab[2][l] * weights2[2]) * det;
						}
					}
				}
				for (int k = 0; k < 4; k++)
				{
					grid.elements[i].P[k] += P[k];
					for (int l = 0; l < 4; l++)
					{
						grid.elements[i].Hbc[k][l] += Hbc[k][l];
					}
				}
				//grid.add_boundary_Condition(i, Hbc, P);
				
#ifdef print==1
				/*
				
				for (int z = 0; z < 4; z++)
				{
					for (int v = 0; v < 4; v++)
						cout << Hbc[z][v] << " ";
					cout << endl;
				}
				cout << endl;
				*/
#endif
			}
		}
	}
}

double calculate_determinant(double x1, double y1, double x2, double y2)
{
	double value = 0;
	double x = 0, y = 0;
	x = (x1 - x2) * (x1 - x2);
	y = (y1 - y2) * (y1 - y2);
	value = sqrt(x + y);
	return value;
}

double* aggregate_vector_P(Grid& grid)
{
	double* aggregation_matrix = new double[grid.get_amount_nodes()];
	for (int j = 0; j < grid.get_amount_nodes(); j++)
	{
		aggregation_matrix[j] = 0.0;
	}
	//Przejscie po ka¿dym elemencie
	for (int i = 0; i < grid.get_amount_elements(); i++)
	{
		for (int l = 0; l < 4; l++)
		{
			aggregation_matrix[grid.get_element(i).get_id_parameter(l) - 1] += grid.elements[i].P[l];//get P matrix at element i in row l
		}
	}
	/*

	cout << "\n\n\n";

	for (int i = 0; i < grid.get_amount_nodes(); i++)
	{
		cout << aggregation_matrix[i] << " ";
	}
	cout << "\n";
	*/
	return aggregation_matrix;

}

double** aggregate_vector_H(Grid& grid)
{
	//return nullptr;
	double** aggregation_matrix = new double* [grid.get_amount_nodes()];
	for (int i = 0; i < grid.get_amount_nodes(); i++)
	{
		aggregation_matrix[i] = new double[grid.get_amount_nodes()];
		for (int j = 0; j < grid.get_amount_nodes(); j++)
		{
			aggregation_matrix[i][j] = 0.0;
		}
	}

	//Przejscie po ka¿dym elemencie
	for (int i = 0; i < grid.get_amount_elements(); i++)
	{
		for (int l = 0; l < 4; l++)
			for (int k = 0; k < 4; k++)
			{
				aggregation_matrix[grid.get_element(i).get_id_parameter(l) - 1][grid.get_element(i).get_id_parameter(k) - 1] += grid.elements[i].H[l][k];//get element i from matrix H in row l & column k
			}
	}
	/*

	for (int i = 0; i < grid.get_amount_nodes(); i++)
	{
		for (int j = 0; j < grid.get_amount_nodes(); j++)
		{
			cout << aggregation_matrix[i][j] << " ";
		}
		cout << "\n";
	}
	*/
	return aggregation_matrix;

}

void calculate_matrix_C(Grid& grid, Elem4 element)
{
	double c = grid.SpecificHeat;
	double ro = grid.Density;
	for (int i = 0; i < grid.get_amount_elements(); i++)
	{
		double C[4][4] = { {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0} };
		for (int j = 0; j < element.schemat * element.schemat; j++) //schemat ca³kowania
		{
			double I[2][2] = { {0,0},{0,0} };
			double Iinv[2][2] = { {0,0},{0,0} };

			jakobian(i, j, I, Iinv, element, grid);
			double detI = I[0][0] * I[1][1] - I[1][0] * I[0][1];

			for (int k = 0; k < 4; k++)
			{
				for (int m = 0; m < 4; m++)
				{
					double data = 0;
					if (element.schemat == 2)
					{
						data = c * ro * (element.N[j][k] * element.N[j][m]) * detI;
					}
					else if (element.schemat == 3)
					{
						data = c * ro * (element.N[j][k] * element.N[j][m]) * detI * element.weight[j / 3] * element.weight[j % 3];
					}
					C[k][m] += data;
				}
			}

		}
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++)
			{
				grid.elements[i].C[j][k] = C[j][k];
			}
		}
	}
		//grid.save_C_matrix(C, i);
		/*
		void Grid::save_C_matrix(double H[4][4], int element)
{
	this->elements[element].set_C_matrix(H);
}
		
		void Element::set_C_matrix(double C[4][4])
{
	
		*/

	
}

double** aggregate_vector_C(Grid& grid)
{
	double** aggregation_matrix = new double* [grid.get_amount_nodes()];
	for (int i = 0; i < grid.get_amount_nodes(); i++)
	{
		aggregation_matrix[i] = new double[grid.get_amount_nodes()];
		for (int j = 0; j < grid.get_amount_nodes(); j++)
		{
			aggregation_matrix[i][j] = 0.0;
		}
	}

	for (int i = 0; i < grid.get_amount_elements(); i++)
	{
		for (int l = 0; l < 4; l++)
			for (int k = 0; k < 4; k++)
			{
				aggregation_matrix[grid.get_element(i).get_id_parameter(l) - 1][grid.get_element(i).get_id_parameter(k) - 1] += grid.elements[i].C[l][k];//get_C_matrix_at(i, l, k);
			}
	}

	//cout << "C matrix\n";

	//for (int i = 0; i < grid.get_amount_nodes(); i++)
	//{
	//	for (int j = 0; j < grid.get_amount_nodes(); j++)
	//	{
	//		cout << aggregation_matrix[i][j] << " ";
	//	}
	//	cout << "\n";
	//}

	return aggregation_matrix;
}

double** devide_matrix_by_number(double** matrix1, double number1, int size)
{
	double** temp = new double* [size];
	for (int i = 0; i < size; i++)
	{
		temp[i] = new double[size];
	}
	for (int i = 0; i < size; i++) 
	{
		for (int j = 0; j < size; j++)
		{
			temp[i][j] = matrix1[i][j] / number1;
		}
	}
	return temp;
}

double** sum_matrix(double** matrix1, double** matrix2, int size)
{
	double** temp = new double* [size];
	for (int i = 0; i < size; i++)
	{
		temp[i] = new double[size];
	}

	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
		{
			temp[i][j] = matrix1[i][j] + matrix2[i][j];
		}
	return temp;
}

double* multiplication_matrix_by_vector(double** matrix1, double* vector1, int size)
{
	double* temp = new double[size];
	for (int i = 0; i < size; i++)
		temp[i] = 0.0;


	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
		{
			temp[i] += matrix1[i][j] * vector1[j];
		}
	return temp;
}

double* sum_vector(double* vector1, double* vector2, int size)
{
	double* temp = new double[size];

	for (int i = 0; i < size; i++)
		temp[i] = vector1[i] + vector2[i];
	return temp;
}

double* Gauss_elimination(double** A, double* B, int n) {
	//double** temp = merge(A, B, N);
	//return Gauss_elimination(temp, N);

	double** temp = new double* [n];
	for (int i = 0; i < n; i++)
	{
		temp[i] = new double[n + 1];
		for (int j = 0; j < n; j++)
		{
			temp[i][j] = A[i][j];
		}
		temp[i][n] = B[i];
	}
	double* handle = Gauss_elimination(temp, n);
	for (int i = 0; i < n; i++)
	{
		delete[] temp[i];
	}
	delete[] temp;


	return handle;
}

double* Gauss_elimination(double** AB, int N) {
	const double accuracy = 1e-15;
	double* result = new double[N];
	int* vector = new int[N + 1];

	for (int i = 0; i < N + 1; i++)
		vector[i] = i;

	for (int i = 0; i < N - 1; i++) {
		bool hasChanged = false;
		int largest = i;

		for (int j = i + 1; j < N; j++)
			if (fabs(AB[i][vector[largest]]) < fabs(AB[i][vector[j]])) {
				hasChanged = true;
				largest = j;
			}

		if (hasChanged) {
			int pom = vector[i];
			vector[i] = vector[largest];
			vector[largest] = pom;
		}

		for (int j = i + 1; j < N; j++) {
			if (fabs(AB[i][vector[i]]) < accuracy)
				return NULL;

			double divisor = AB[j][vector[i]] / AB[i][vector[i]];

			for (int k = i + 1; k < N + 1; k++)
				AB[j][vector[k]] -= (AB[i][vector[k]] * divisor);
		}
	}

	for (int i = N - 1; i >= 0; i--) {
		if (fabs(AB[i][vector[i]]) < accuracy)
			return NULL;

		for (int j = N - 1; j > i; j--)
			AB[i][N] -= AB[i][vector[j]] * result[vector[j]];

		result[vector[i]] = AB[i][N] / AB[i][vector[i]];
	}

	return result;
}