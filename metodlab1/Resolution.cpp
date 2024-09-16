#include "Resolution.h"

Resolution::Resolution(){}

Resolution::Resolution(vector <vector <double>> A, vector <double> B, int size)
{
	this->A = A;
	this->B = B;
	this->size = size;
}

void Resolution::FindLU()
{
	this->U = this->A;

	for (int i = 0; i < this->size; i++)
	{
		for (int j = i; j < this->size; j++)
		{
			this->L[j][i] = this->U[j][i] / this->U[i][i];
		}
	}

	for (int k = 1; k < this->size; k++)
	{
		for (int i = k - 1; i < this->size; i++)
		{
			for (int j = i; j < this->size; j++)
			{
				this->L[j][i] = this->U[j][i] / this->U[i][i];
			}
		}

		for (int i = k; i < this->size; i++)
		{
			for (int j = k - 1; j < this->size; j++)
			{
				this->U[i][j] = this->U[i][j] - this->L[i][k - 1] * this->U[k - 1][j];
			}
		}
	}
}

vector <double> Resolution::FindSolution(vector <double> B)
{
	double sum;
	vector <double>y, x;
	y.push_back(B[0]);
	for (int i = 1; i < this->size; i++)
	{
		sum = this->L[i][0] * y[0];
		for (int j = i - 1; j > 0; j--)
		{
			sum += this->L[i][j] * y[j];
		}
		y.push_back(B[i] - sum);
	}
	x.push_back(y[this->size - 1] / this->U[this->size - 1][this->size - 1]);
	for (int i = 2; i <= this->size; i++)
	{
		sum = 0;
		for (int j = 1; j < i; j++)
		{
			sum += this->U[this->size - i][this->size - j] * x[j - 1];
		}
		x.push_back((y[this->size - i] - sum) / (this->U[this->size - i][this->size - i]));
	}
	reverse(x.begin(), x.end());
	return x;
}


void Resolution::FindResult()
{
	this->Result = FindSolution(this->B);
}

void Resolution::ShowResult()
{
	setlocale(LC_ALL, "ru");
	printf("\nРезультат:\n");
	for (int i = 0; i < this->size; i++)
	{
		printf("x%d: %lf\n",i+1,this->Result[i]);
	}
}

void Resolution::ShowMatrix(vector <vector <double>> Matrix)
{
	for (int i = 0; i < this->size; i++)
	{
		for (int j = 0; j < this->size; j++)
		{
			printf("\t%lf\t", Matrix[i][j]);
		}
		printf("\n");
	}
}

void Resolution::InputMatrixA()
{
	setlocale(LC_ALL, "ru");
	if (this->size == 0)
	{
		printf("\nНе указаны размеры матрицы!!!\n");
	}
	else
	{
		double inp;
		vector <double> str;
		printf("\nВведите матрицу\n");
		for (int i = 0; i < this->size; i++)
		{
			str.push_back(0);
		}
		for (int i = 0; i < this->size; i++)
		{
			this->A.push_back(str);
			this->L.push_back(str);
			this->U.push_back(str);
			for (int j = 0; j < this->size; j++)
			{
				inp = 0;
				scanf("%lf", &inp);
				this->A[i][j] = inp;
				this->L[i][j] = 0;
				this->U[i][j] = 0;
			}
		}
	}
}

void Resolution::InputB()
{
	setlocale(LC_ALL, "ru");
	if (this->size == 0)
	{
		printf("\nНе указаны размеры матрицы!!!\n");
	}
	else
	{
		double inp;
		vector <double> str;
		printf("\nВведите столбец b\n");
		for (int j = 0; j < this->size; j++)
		{
			scanf("%lf", &inp);
			this->B.push_back(inp);
		}
	}
}

void Resolution::InputSize()
{
	setlocale(LC_ALL, "ru");
	printf("\nВведите размерность матрицы\n");
	int inp;
	scanf("%d", &inp);
	this->size = inp;
}

vector <vector <double>> Resolution::GetL()
{
	return this->L;
}

vector <vector <double>> Resolution::GetU()
{
	return this->U;
}

double Resolution::FindDeterminant()
{
	double result = 1;
	for (int i = 0; i < this->size; i++)
	{
		result *= this->GetU()[i][i];
	}

	printf("\nОпределитель = %lf\n", result);
	return result;
}

vector <double> Resolution::GetB()
{
	return this->B;
}

vector <vector <double>> Resolution::FindInverseMatrix()
{
	vector <vector <double>> result, E;
	vector <double> column;
	for (int i = 0; i < this->size; i++)
	{
		column.push_back(0);
	}
	for (int i = 0; i < this->size; i++)
	{
		E.push_back(column);
		E[i][i] = 1;
	}
	for (int i = 0; i < this->size; i++)
	{
		result.push_back(FindSolution(E[i]));
	}

	return result;
}


vector <double> Resolution::FindDiscrepancies() 
{
	vector <double> z;
	double sum;
	for (int i = 0; i < this->size; i++)
	{
		sum = 0;
		for (int j = 0; j < this->size; j++)
		{
			sum += this->A[i][j] * this->Result[j];
		}
		z.push_back(this->B[i] - sum);
	}

	setlocale(LC_ALL, "ru");
	printf("\nНевязки:\n");
	for (int i = 0; i < this->size; i++)
	{
		printf("z%d: %.14lf\n", i + 1, z[i]);
	}
	return z;
}

double Resolution::FindCertainty(vector <vector <double>> InverseMatrix)
{
	double norm1 = 0, norm2 = 0, sum1, sum2;
	for (int i = 0; i < this->size; i++)
	{
		sum1 = 0;
		sum2 = 0;
		for (int j = 0; j < this->size; j++)
		{
			if (this->A[i][j]<0)
			{
				sum1 -= A[i][j];
			}
			else
			{
				sum1 += A[i][j];
			}

			if (InverseMatrix[i][j] < 0)
			{
				sum2 -= InverseMatrix[i][j];
			}
			else
			{
				sum2 += InverseMatrix[i][j];
			}
		}

		if (norm1<sum1)
		{
			norm1 = sum1;
		}

		if (norm2 < sum2)
		{
			norm2 = sum2;
		}
	}
	setlocale(LC_ALL, "ru");
	printf("\nЧисло обусловленности = %lf\n", norm1 * norm2);
	return norm1*norm2;
}