#include "Element.h"
#include "Grid.h"
#include "Node.h"
#include "Elem4.h"
#include "Matrix.h"
#include "Jakobian.h"
Element::Element()
{
	for (int i = 0; i < 4; i++)
	{
		ID[i] = 0;
	}

	for (int i = 0; i < 4; i++)
	{
		P[i] = 0;
		for (int j = 0; j < 4; j++)
		{
			H[i][j] = 0;
			Hbc[i][j] = 0;
			C[i][j] = 0;
		}
	}
}

Element::Element(Element& element)
{
	for (int i = 0; i < 4; i++)
	{
		this->ID[i] = element.ID[i];
	}

	for (int i = 0; i < 4; i++)
	{
		P[i] = element.P[i];
		for (int j = 0; j < 4; j++)
		{
			this->H[i][j] = element.H[i][j];
			Hbc[i][j] = element.Hbc[i][j];
			C[i][j] = element.C[i][j];
		}
	}
}

Element::Element(int* array)
{
	for (int i = 0; i < 4; i++)
	{
		ID[i] = array[i];
	}

	for (int i = 0; i < 4; i++)
	{
		P[i] = 0;
		for (int j = 0; j < 4; j++)
		{
			H[i][j] = 0;
			Hbc[i][j] = 0;
			C[i][j] = 0;
		}
	}
}

Element::Element(int id1, int id2, int id3, int id4)
{
	ID[0] = id1;
	ID[1] = id2;
	ID[2] = id3;
	ID[3] = id4;

	for (int i = 0; i < 4; i++)
	{
		P[i] = 0;
		for (int j = 0; j < 4; j++)
		{
			H[i][j] = 0;
			Hbc[i][j] = 0;
			C[i][j] = 0;
		}
	}
}

Element::Element(int* array, double H1[4][4])
{
	for (int i = 0; i < 4; i++)
	{
		ID[i] = array[i];
	}

	for (int i = 0; i < 4; i++)
	{
		P[i] = 0;
		for (int j = 0; j < 4; j++)
		{
			H[i][j] = H1[i][j];
			Hbc[i][j] = 0;
			C[i][j] = 0;
		}
	}
}

Element::Element(int id1, int id2, int id3, int id4, double H1[4][4])
{
	ID[0] = id1;
	ID[1] = id2;
	ID[2] = id3;
	ID[3] = id4;

	for (int i = 0; i < 4; i++)
	{
		P[i] = 0;
		for (int j = 0; j < 4; j++)
		{
			H[i][j] = H1[i][j];
			Hbc[i][j] = 0;
			C[i][j] = 0;
		}
	}
}


void Element::set_id(int id1, int id2, int id3, int id4)
{
	ID[0] = id1;
	ID[1] = id2;
	ID[2] = id3;
	ID[3] = id4;
}

int Element::get_id1()
{
	return ID[0];
}

int Element::get_id2()
{
	return ID[1];
}

int Element::get_id3()
{
	return ID[2];
}

int Element::get_id4()
{
	return ID[3];
}

int Element::get_id_parameter(int number)
{
	return ID[number];
}

void Element::set_H_matrix(double H1[4][4])
{
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
		{
			double temp = H1[i][j];
			this->H[i][j] = temp;
		}
}