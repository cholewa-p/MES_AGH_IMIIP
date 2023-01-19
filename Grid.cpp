#include "Element.h"
#include "Grid.h"
#include "Node.h"
#include "Elem4.h"
#include "Matrix.h"
#include "Jakobian.h"
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
using namespace std;
Grid::Grid():nElements(NULL),nNodes(NULL),t(NULL)
{
	SimulationTime = NULL;
	SimulationStepTime = NULL;
	Conductivity = NULL;
	Alfa = NULL;
	Tot = NULL;
	InitialTemp = NULL;
	Density = NULL;
	SpecificHeat = NULL;
	nodes = nullptr;
	elements = nullptr;
}
Grid::~Grid()
{
	delete[] nodes;
	delete[] elements;
}
void Grid::set_Grid(string path)
{
	fstream plik;
	plik.open(path);
	if (plik.good())
	{
		string line;
		string name, value;
		int amount_of_nodes = 0;
		int amount_of_element = 0;
		int loads = 0;
		while (getline(plik, line))
		{
			if (line[0] == '*')
			{
				loads++;
				getline(plik, line);
			}
			if (loads == 0)
			{
				stringstream load(line);
				load >> name;
				load >> value;
				if (name == "SimulationTime")
					SimulationTime = stoi(value);
				else if (name == "SimulationStepTime")
					SimulationStepTime = stoi(value);
				else if (name == "Conductivity")
					Conductivity = stoi(value);
				else if (name == "Alfa")
					Alfa = stoi(value);
				else if (name == "Tot")
					Tot = stoi(value);
				else if (name == "InitialTemp")
					InitialTemp = stoi(value);
				else if (name == "Density")
					Density = stoi(value);
				else if (name == "SpecificHeat")
					SpecificHeat = stoi(value);
				else if (name == "Nodes_number")
				{
					amount_of_nodes = stoi(value);
					nNodes = amount_of_nodes;
					nodes = new Node[amount_of_nodes];
					
				}
				else if (name == "Elements_number")
				{
					amount_of_element = stoi(value);
					nElements = amount_of_element;
					elements = new Element[amount_of_element];
				}
				else
					cout << "Bledne dane" << name << " " << value << endl;
			}
			if (loads == 1)
			{
				string node;
				string x, y;
				stringstream load(line);
				load >> node;
				load >> x;
				load >> y;
				node.pop_back();
				x.pop_back();
				Node* handle = &nodes[(stoi(node) - 1)];
				handle->set_x(stof(x));
				handle->set_y(stof(y));
				handle->set_bc(0);
			}
			if (loads == 2)
			{
				string element;
				string nod1, nod2, nod3, nod4;
				stringstream load(line);
				load >> element;
				load >> nod1;
				load >> nod2;
				load >> nod3;
				load >> nod4;
				element.pop_back();
				nod1.pop_back();
				nod2.pop_back();
				nod3.pop_back();

				Element* handle = &elements[(stoi(element) - 1)];
				handle->set_id(stoi(nod1), stoi(nod2), stoi(nod3), stoi(nod4));
			}
			if (loads == 3)
			{
				string border;
				stringstream load(line);
				while (true)
				{
					load >> border;
					if (border[border.size() - 1] == ',')
					{
						border.pop_back();
						Node* handle = &nodes[(stoi(border) - 1)];
						handle->set_bc(1);
					}
					else
					{
						Node* handle = &nodes[(stoi(border) - 1)];
						handle->set_bc(1);
						break;
					}
				}
			}
		}
		t = new double[amount_of_nodes];
	}
	else
		cout << "Input ERROR\n";
}
void Grid::show_nodes()
{
	for (int i = 0; i < nNodes; i++)
	{
		cout << "Node No " << i + 1 << ":\n" << "x:" << this->nodes[i].get_x() << " y: " << this->nodes[i].get_y();
		string a = (this->nodes[i].get_bc() == 1) ? "TRUE" : "FALSE";
		cout<<"\nBC: " << a << endl;
	}
}
void Grid::show_elements()
{
	for (int i = 0; i < nElements; i++)
	{
		cout << "Element No " << i + 1 << ":\n";
		cout << this->elements[i].get_id1() << " - ";
		cout << this->elements[i].get_id2() << " - ";
		cout << this->elements[i].get_id3() << " - ";
		cout << this->elements[i].get_id4() << endl;
	}
}
int Grid::get_amount_nodes()
{
	return nNodes;
}
int Grid::get_amount_elements()
{
	return nElements;
}
Node Grid::get_node(int number)
{
	return nodes[number];
}
Element Grid::get_element(int number)
{
	return elements[number];
}
Node* Grid::get_orginal_node(int number)
{
	return &nodes[number];
}
Element* Grid::get_orginal_element(int number)
{
	return &elements[number];
}
void Grid::setnNodes(int n)
{
	nNodes = n;
	if (nodes != nullptr)
	{
		delete[] nodes;
	}
	nodes = new Node[nNodes];
}
void Grid::setnElements(int n)
{
	nElements = n;
	if (elements != nullptr)
	{
		delete[] elements;
	}
	elements = new Element[nElements];
}
void Grid::save_H_Matrix(double H[4][4], int element)
{
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			elements[element].H[i][j] = H[i][j];
}
void Grid::show_H_Matrix() {
	for (int i = 0; i < get_amount_elements(); i++)
	{
		cout << "Element " << i + 1 << endl;
		for (int l = 0; l < 4; l++)
		{
			for (int m = 0; m < 4; m++)
				cout << elements[i].H[l][m] << " ";
			cout << endl;
		}
		cout << endl;
	}
}
void Grid::show_Hbc_Matrix()
{
	for (int i = 0; i < get_amount_elements(); i++)
	{
		//cout << "Element " << i + 1 << endl;
		if (get_node(get_element(i).get_id1() - 1).get_bc() == 1)
		{
			if (get_node(get_element(i).get_id2() - 1).get_bc() == 1)
			{
				for (int z = 0; z < 4; z++)
				{
					for (int v = 0; v < 4; v++)
						cout << elements->Hbc[z][v] << " ";
					cout << endl;
				}
				cout << endl;

			}

			if (get_node(get_element(i).get_id4() - 1).get_bc() == 1)
			{
				for (int z = 0; z < 4; z++)
				{
					for (int v = 0; v < 4; v++)
						cout << elements->Hbc[z][v] << " ";
					cout << endl;
				}
				cout << endl;
				
			}
		}

		if (get_node(get_element(i).get_id3() - 1).get_bc() == 1)
		{
			if (get_node(get_element(i).get_id2() - 1).get_bc() == 1)
			{
				for (int z = 0; z < 4; z++)
				{
					for (int v = 0; v < 4; v++)
						cout << elements->Hbc[z][v] << " ";
					cout << endl;
				}
				cout << endl;
			}

			if (get_node(get_element(i).get_id4() - 1).get_bc() == 1)
			{
				for (int z = 0; z < 4; z++)
				{
					for (int v = 0; v < 4; v++)
						cout << elements->Hbc[z][v] << " ";
					cout << endl;
				}
				cout << endl;
			}
		}
	}
}
void Grid::initialise_Temperature() 
{
	for (int i = 0; i < get_amount_nodes(); i++)
	{
		t[i] = InitialTemp;
	}
}
void Grid::sum_BC() 
{
	for (int i = 0; i < nElements; i++)
	{
		for (int j = 0; j < 4; j++)
			for (int k = 0; k < 4; k++)
			{
				elements[i].H[j][k] += elements[i].Hbc[j][k];
			}
	}
}


