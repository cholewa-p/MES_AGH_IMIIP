#pragma once
#include "Elem4.h"
#include "Grid.h"
#include <iostream>
using namespace std;
void jakobian(int i, int j, double I[2][2], double Iinv[2][2], Elem4 element, Grid& grid);