#pragma once
#include "Node.h"
#include "Element.h"
#include <string>
using namespace std;
class Grid
{
public:
	int SimulationTime;
	int SimulationStepTime;
	int Conductivity;
	int Alfa;
	int Tot;
	int InitialTemp;
	int Density;
	int SpecificHeat;
	double* t;

	int nNodes;
	int nElements;

	Node* nodes;
	Element* elements;
	Grid();
	void set_Grid(string path);
	~Grid();
	void show_nodes();
	void show_elements();
	int get_amount_nodes();
	int get_amount_elements();
	Node get_node(int number);
	Element get_element(int number);
	Node* get_orginal_node(int number);
	Element* get_orginal_element(int number);
	void setnNodes(int n);
	void setnElements(int n);
	void save_H_Matrix(double H[4][4], int element);
	void show_H_Matrix();
	void show_Hbc_Matrix();
	void initialise_Temperature();
	void sum_BC();
};
