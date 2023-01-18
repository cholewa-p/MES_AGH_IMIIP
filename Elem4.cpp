#include <iostream>
#include "Elem4.h"


Elem4::Elem4(int n):N(NULL),eta(NULL),ksi(NULL),weight(NULL)
{
	schemat = n;
	int nodes=0;
	if (schemat == 2)
	{
		nodes = 2 * 2;

		eta = new double* [nodes];
		ksi = new double* [nodes];
		N = new double* [nodes];

		double points[4][2] = {
			{-1 / sqrt(3), -1 / sqrt(3)},
			{1 / sqrt(3), -1 / sqrt(3)},
			{1 / sqrt(3), 1 / sqrt(3)},
			{-1 / sqrt(3), 1 / sqrt(3)}
		};
		weight = new double[2];

		weight[0] = 1.0;
		weight[1] = 1.0;

		for (int i = 0; i < nodes; i++)
		{
			eta[i] = new double[nodes];
			ksi[i] = new double[nodes];
			N[i] = new double[nodes];
		}

		for (int i = 0; i < nodes; i++)
		{

			ksi[i][0] = -0.25 * (1 - points[i][1]);
			ksi[i][1] = 0.25 * (1 - points[i][1]);
			ksi[i][2] = 0.25 * (1 + points[i][1]);
			ksi[i][3] = -0.25 * (1 + points[i][1]);

			eta[i][0] = -0.25 * (1 - points[i][0]);
			eta[i][1] = -0.25 * (1 + points[i][0]);
			eta[i][2] = 0.25 * (1 + points[i][0]);
			eta[i][3] = 0.25 * (1 - points[i][0]);
		}

		for (int i = 0; i < 4; i++)
		{
			N[i][0] = 0.25 * (1.0 - points[i][0]) * (1.0 - points[i][1]);
			N[i][1] = 0.25 * (1.0 + points[i][0]) * (1.0 - points[i][1]);
			N[i][2] = 0.25 * (1.0 + points[i][0]) * (1.0 + points[i][1]);
			N[i][3] = 0.25 * (1.0 - points[i][0]) * (1.0 + points[i][1]);
		}
	}

	if (schemat == 3)
	{
		nodes = 3 * 3;

		eta = new double* [nodes];
		ksi = new double* [nodes];
		N = new double* [nodes];

		weight = new double[3];

		weight[0] = 5.0 / 9.0;
		weight[1] = 8.0 / 9.0;
		weight[2] = 5.0 / 9.0;

		double points[9][2] = {
			{-sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0)},
			{-sqrt(3.0 / 5.0), 0},
			{-sqrt(3.0 / 5.0), sqrt(3.0 / 5.0)},
			{0, -sqrt(3.0 / 5.0)},
			{0,0},
			{0, sqrt(3.0 / 5.0)},
			{sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0)},
			{sqrt(3.0 / 5.0), 0},
			{sqrt(3.0 / 5.0),sqrt(3.0 / 5.0)}
		};
		for (int i = 0; i < nodes; i++)
		{
			eta[i] = new double[4];
			ksi[i] = new double[4];
			N[i] = new double[4];
		}

		for (int i = 0; i < nodes; i++)
		{
			ksi[i][0] = -0.25 * (1 - points[i][1]);
			ksi[i][1] = 0.25 * (1 - points[i][1]);
			ksi[i][2] = 0.25 * (1 + points[i][1]);
			ksi[i][3] = -0.25 * (1 + points[i][1]);

			eta[i][0] = -0.25 * (1 - points[i][0]);
			eta[i][1] = -0.25 * (1 + points[i][0]);
			eta[i][2] = 0.25 * (1 + points[i][0]);
			eta[i][3] = 0.25 * (1 - points[i][0]);
		}

		for (int i = 0; i < nodes; i++)
		{
			N[i][0] = 0.25 * (1.0 - points[i][0]) * (1.0 - points[i][1]);
			N[i][1] = 0.25 * (1.0 + points[i][0]) * (1.0 - points[i][1]);
			N[i][2] = 0.25 * (1.0 + points[i][0]) * (1.0 + points[i][1]);
			N[i][3] = 0.25 * (1.0 - points[i][0]) * (1.0 + points[i][1]);
		}
	}

}

Elem4::Elem4(Elem4& element)
{
	this->schemat = element.schemat;

	eta = new double* [schemat * schemat];
	ksi = new double* [schemat * schemat];
	N = new double* [schemat * schemat];

	weight = new double[element.schemat];

	for (int i = 0; i < schemat; i++)
	{
		weight[i] = element.weight[i];
	}

	for (int i = 0; i < schemat * schemat; i++)
	{
		eta[i] = new double[4];
		ksi[i] = new double[4];
		N[i] = new double[4];
	}

	for (int i = 0; i < schemat * schemat; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			eta[i][j] = element.eta[i][j];
			ksi[i][j] = element.ksi[i][j];
			N[i][j] = element.N[i][j];
		}
	}
}

Elem4::~Elem4()
{
	for (int i = 0; i < schemat * schemat; i++)
	{
		delete[] eta[i];
		delete[] ksi[i];
		delete[] N[i];
	}

	delete[] weight;
	delete[] eta;
	delete[] ksi;
	delete[] N;
	
}
