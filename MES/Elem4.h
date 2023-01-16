#pragma once
class Elem4
{
public:
	int schemat;

	double** eta;
	double** ksi;
	double* weight;

	double** N;

	Elem4(int nodes);

	Elem4(Elem4& element);

	~Elem4();
};