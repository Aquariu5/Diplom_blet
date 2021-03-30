#ifndef __MY_MATH_H
#define __MY_MATH_H

#include <iostream>
#include <vector>
#include <fstream>
#include "adept.h"
static double PI = 3.1415926535897932384626433832795;

// Печать одномерного вектора
template <class T>
inline void Print(const std::vector<T> &A, int n) {
	for (int i = 0; i < n; i++) {
		std::cout << A[i] << std::endl;
	}
}

//Печать двуменрого вектора
template <class T>
inline void Print2D(const std::vector<std::vector<T> > &A, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			std::cout.width(10);
			std::cout << A[i][j];
		}
		std::cout << std::endl;
	}
}

// Считывание с файла IRAP и запись в массив
//typedef adept::adouble T;
//template <class T>
inline void FromIrap(std::string filename, std::vector<std::vector<adept::adouble>> &arr) {
	int ch = 0;
	std::ifstream fin(filename);
	int rows, cols;
	double num;
	double tmp;

	fin >> tmp; // первое число 1ой строки
	fin >> cols; // второе число 1ой строки - и есть число столбцов
	while ((ch = fin.get()) != '\n') {
		fin >> tmp;
	} // 1я строка кончилась


	while ((ch = fin.get()) != '\n') {
		fin >> tmp;
	} // 2я строка кончилась


	fin >> rows; // первое число 3ей строки
	//cout << "Strok: " << rows << endl;
	//cout << "Stolbcov: " << cols << endl;

	while ((ch = fin.get()) != '\n') {
		fin >> tmp;
	} // 3я кончилась

	while ((ch = fin.get()) != '\n') {
		fin >> tmp;
	} // 4я кончилась

	arr.assign(rows, std::vector<adept::adouble>(cols)); // задаем размерность двумерного вектора

	std::vector<std::vector<double>>help;

	help.assign(rows, std::vector<double>(cols)); // задаем размерность двумерного вектора
	// теперь пошли числа
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			fin >> help[i][j];
		}
	}
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			arr[i][j] = help[i][j];
		}
	}
	//Print2D(arr, rows);
}
template <class T>
inline void Print2D(const std::vector<std::vector<T> > &A, int m, int n) {
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			std::cout.width(10);
			std::cout << A[i][j];
		}
		std::cout << std::endl;
	}
}
// Нахождение нормы двух векторов
template <class T>
inline T Norma(const std::vector<T> &xnew, const std::vector<T> &xold, int n) {
	T max = abs(xnew[0] - xold[0]);
	for (int i = 1; i < n; i++) {
		//max = std::max(max, abs(xnew[i] - xold[i]));
		if (abs(xnew[i] - xold[i]) > max) {
			max = abs(xnew[i] - xold[i]);
		}
	}
	return max;
}

//Нахождение определителя матрицы
double det(const std::vector<std::vector<double> > &matrix, int n) //квадратная матрица размера n*n
{
	//double **B = clone(matrix, n);
	std::vector< std::vector<double>> B(n);
	for (int i = 0; i < n; i++)
		B[i].resize(n);

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			B[i][j] = matrix[i][j];
	//приведение матрицы к верхнетреугольному виду
	for (int step = 0; step < n - 1; step++)
		for (int row = step + 1; row < n; row++)
		{
			double coeff = -B[row][step] / B[step][step]; //метод Гаусса
			for (int col = step; col < n; col++)
				B[row][col] += B[step][col] * coeff;
		}
	//Рассчитать определитель как произведение элементов главной диагонали
	double Det = 1;
	for (int i = 0; i < n; i++)
		Det *= B[i][i];
	//Очистить память
	return Det;
}

//Метод прогонки для решения СЛАУ с 3diag матрицей
 inline void Progonka(std::vector<double> &A, std::vector<double> &C, std::vector<double> &B, std::vector<double> &F, std::vector<double> &x, int n) {
	std::vector<double> alpha(n);
	std::vector<double> beta(n);

	alpha[0] = 0;
	beta[0] = 0;
	A[0] = 0;
	B[n - 1] = 0;
	
	for (int i = 0; i < n - 1; i++) {
		alpha[i + 1] = -B[i] / (A[i] * alpha[i] + C[i]);
		beta[i + 1] = (F[i] - A[i] * beta[i]) / (A[i] * alpha[i] + C[i]);
	}

	x[n - 1] = (F[n - 1] - A[n - 1] * beta[n - 1]) / (C[n - 1] + A[n - 1] * alpha[n - 1]);
	
	for (int i = n - 2; i >= 0; i--) {
		x[i] = alpha[i + 1] * x[i + 1] + beta[i + 1];
	}
}

 // Метод Гаусса для решения СЛАУ
 inline void gauss(std::vector<std::vector<double>> a, std::vector<double> y, int n, std::vector<double> &x)
{
	double  max;
	int k, index;
	const double eps = 0.00001;  // точность
								 //x = new double[n];
	k = 0;
	while (k < n)
	{
		// Поиск строки с максимальным a[i][k]
		max = abs(a[k][k]);
		index = k;
		for (int i = k + 1; i < n; i++)
		{
			if (abs(a[i][k]) > max)
			{
				max = abs(a[i][k]);
				index = i;
			}
		}
		// Перестановка строк
		if (max < eps)
		{
			// нет ненулевых диагональных элементов
			std::cout << "Решение получить невозможно из-за нулевого столбца ";
			std::cout << index << " матрицы A" << std::endl;
			return;
		}
		for (int j = 0; j < n; j++)
		{
			double temp = a[k][j];
			a[k][j] = a[index][j];
			a[index][j] = temp;
		}
		double temp = y[k];
		y[k] = y[index];
		y[index] = temp;
		// Нормализация уравнений
		for (int i = k; i < n; i++)
		{
			double temp = a[i][k];
			if (abs(temp) < eps) continue; // для нулевого коэффициента пропустить
			for (int j = 0; j < n; j++)
				a[i][j] = a[i][j] / temp;
			y[i] = y[i] / temp;
			if (i == k)  continue; // уравнение не вычитать само из себя
			for (int j = 0; j < n; j++)
				a[i][j] = a[i][j] - a[k][j];
			y[i] = y[i] - y[k];
		}
		k++;
	}
	// обратная подстановка
	for (k = n - 1; k >= 0; k--)
	{
		x[k] = y[k];
		for (int i = 0; i < k; i++)
			y[i] = y[i] - a[i][k] * x[k];
	}
	//return;
}


 // Вычисление обратной матрицы
 inline void Inversion(std::vector< std::vector<double>> A, std::vector< std::vector<double>> &E, int N)
{
	double temp;

	/*std::vector< std::vector<double>> E(N);

	for (int i = 0; i < N; i++)
	E[i].resize(N);*/

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
		{
			E[i][j] = 0.0;

			if (i == j)
				E[i][j] = 1.0;
		}

	for (int k = 0; k < N; k++)
	{
		temp = A[k][k];

		for (int j = 0; j < N; j++)
		{
			A[k][j] /= temp;
			E[k][j] /= temp;
		}

		for (int i = k + 1; i < N; i++)
		{
			temp = A[i][k];

			for (int j = 0; j < N; j++)
			{
				A[i][j] -= A[k][j] * temp;
				E[i][j] -= E[k][j] * temp;
			}
		}
	}

	for (int k = N - 1; k > 0; k--)
	{
		for (int i = k - 1; i >= 0; i--)
		{
			temp = A[i][k];

			for (int j = 0; j < N; j++)
			{
				A[i][j] -= A[k][j] * temp;
				E[i][j] -= E[k][j] * temp;
			}
		}
	}

	/*for (int i = 0; i < N; i++)
	for (int j = 0; j < N; j++)
	A[i][j] = E[i][j];*/

	/*for (int i = 0; i < N; i++)
	delete[] E[i];

	delete[] E;*/
}

//Получение одномерного вектора из двумерного
template <class T>
void
encode1d(const std::vector< std::vector <T> > &x_2d, std::vector<T> &x_1d, int n) {

	int at = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			x_1d[at++] = x_2d[i][j];
		}
	}
}

//Получение двумерного вектора из одномерного
template <class T>
void
decode2d(const std::vector<T> &x_1d, std::vector<std::vector<T>> &x_2d, int n) {

	int at = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			x_2d[i][j] = x_1d[at++];
		}
	}
}
template <class T>
void InFile2D(const std::vector< std::vector <T> > &x_2d, int n, std::ofstream &fout) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			fout.width(20);
			fout << x_2d[i][j];
		}
		fout << std::endl;
	}
}

#endif