#include "Element.h"
#include "Grid.h"
#include "Node.h"
#include "Elem4.h"
#include "Matrix.h"
#include "Jakobian.h"
#include <iostream>
using namespace std;
void jakobian(int i, int j, double I[2][2], double Iinv[2][2], Elem4 element, Grid& grid)
{
	double points[4][2] = { {0,0}, {0,0}, {0, 0}, {0,0} };

	points[0][0] = grid.get_node(grid.get_element(i).get_id1() - 1).get_x();
	points[0][1] = grid.get_node(grid.get_element(i).get_id1() - 1).get_y();

	points[1][0] = grid.get_node(grid.get_element(i).get_id2() - 1).get_x();
	points[1][1] = grid.get_node(grid.get_element(i).get_id2() - 1).get_y();

	points[2][0] = grid.get_node(grid.get_element(i).get_id3() - 1).get_x();
	points[2][1] = grid.get_node(grid.get_element(i).get_id3() - 1).get_y();

	points[3][0] = grid.get_node(grid.get_element(i).get_id4() - 1).get_x();
	points[3][1] = grid.get_node(grid.get_element(i).get_id4() - 1).get_y();
	for (int k = 0; k < 4; k++)
	{
		I[0][0] += element.ksi[j][k] * points[k][0];
		I[1][1] += element.eta[j][k] * points[k][1];
		I[0][1] += element.ksi[j][k] * points[k][1];
		I[1][0] += element.eta[j][k] * points[k][0];
	}

	double detI = I[0][0] * I[1][1] - I[1][0] * I[0][1];

	Iinv[1][1] = (1.0 / detI) * I[0][0];
	Iinv[0][0] = (1.0 / detI) * I[1][1];
	Iinv[0][1] = (1.0 / detI) * -I[0][1];
	Iinv[1][0] = (1.0 / detI) * -I[1][0];
}