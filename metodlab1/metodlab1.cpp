#include "Resolution.h"


int main()
{
	setlocale(LC_ALL, "ru");
	Resolution res;
	res.InputSize();
	res.InputMatrixA();
	res.InputB();
	res.FindLU();
	res.FindResult();
	printf("Нижняя матрица:\n");
	res.ShowMatrix(res.GetL());
	printf("Верхняя матрица:\n");
	res.ShowMatrix(res.GetU());
	res.ShowResult();
	printf("Обратная матрица:\n");
	res.ShowMatrix(res.FindInverseMatrix());
	res.FindDeterminant();
	res.FindDiscrepancies();
	res.FindCertainty(res.FindInverseMatrix());
	return 0;
}
