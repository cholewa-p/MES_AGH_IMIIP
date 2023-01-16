#pragma once
#include "Grid.h"
#include "Matrix.h"

void simulation(double** H, double** C, double* P, double* t_start, int time, int step, int n);
double* calculate_temperature(double** H, double** C, double* P, double* t_start, int step, int n);
double max_temperature(double* temperature, int n);
double min_temperature(double* temperature, int n);