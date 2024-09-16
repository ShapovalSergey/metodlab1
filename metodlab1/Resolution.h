#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <stdio.h>
#include <vector>
#include <locale.h>

using namespace std;

class Resolution
{
private:
	vector <vector <double>> A, L, U;
	vector <double> B, Result;
	int size=0;
public: 
	Resolution();
	Resolution(vector <vector <double>> A, vector <double> B, int size);
	void FindLU();
	vector <double> FindSolution(vector <double> B);
	void FindResult();
	void ShowResult();
	void ShowMatrix(vector <vector <double>> Matrix);
	void InputMatrixA();
	void InputB();
	void InputSize();
	vector <vector <double>> GetL();
	vector <vector <double>> GetU();
	vector <vector <double>> FindInverseMatrix();
	vector <double> GetB();
	vector <double> FindDiscrepancies(); //Невязки
	double FindDeterminant();
	double FindCertainty(vector <vector <double>> InverseMatrix); //Определенность
};

