#pragma once
#include "Elem4.h"
#include "Grid.h"
#include <iostream>
using namespace std;
void calculate_matrix_H(Grid& grid, Elem4 element);
double calculate_determinant(double x1, double y1, double x2, double y2);
void calculate_matrix_Hbc(Grid& grid, Elem4 element);
double* aggregate_vector_P(Grid& grid);
double** aggregate_vector_H(Grid& grid);
void calculate_matrix_C(Grid& grid, Elem4 element);
double** aggregate_vector_C(Grid& grid);
double** divide_matrix_by_number(double** matrix1, double number1, int size);
double** sum_matrix(double** matrix1, double** matrix2, int size);
double* multiplication_matrix_by_vector(double** matrix1, double* vector1, int size);
double* sum_vector(double* vector1, double* vector2, int size);

double* Gauss_elimination(double** A, double* B, int n);
double* Gauss_elimination(double** AB, int N);