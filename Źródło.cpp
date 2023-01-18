#include <iostream>
#include "Element.h"
#include "Grid.h"
#include "Node.h"
#include "Elem4.h"
#include "Matrix.h"
#include "Jakobian.h"
#include "Simulation.h"
using namespace std;

int main()
{
	
	Grid grid;
	
	grid.set_Grid("Test2_4_4_MixGrid.txt");
	//grid.set_Grid("Test1_4_4.txt");
	//grid.show_elements();
	//grid.show_nodes();
	Elem4 element(2);
	//element.show_Elem4();
	calculate_matrix_H(grid, element);
	//grid.show_H_Matrix();
	calculate_matrix_Hbc(grid,element);
	double* P = aggregate_vector_P(grid);
	grid.initialise_Temperature();
	grid.sum_BC();
	//grid.show_H_Matrix();
	double** H = aggregate_vector_H(grid);
	calculate_matrix_C(grid,element);
	double** C = aggregate_vector_C(grid);
	simulation(H, C, P, grid.t, grid.SimulationTime, grid.SimulationStepTime, grid.get_amount_nodes());
	return 0;
}
