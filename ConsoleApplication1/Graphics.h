#pragma once

#include "Input.h"
using namespace EquationData;

void Errors(BurgersSystem& iburgersSystem, vector<double**>& iresult);
void ExactAndNumericalSolutions(BurgersSystem& iburgersSystem, vector<double**>& iresult);
void FieldVelocity(BurgersSystem& iburgersSystem, vector<double**>& iresult);
void Compare(BurgersSystem& burgersSystem, const vector<double**>& result);
