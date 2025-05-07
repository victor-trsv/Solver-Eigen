#include <iostream>
#include <vector>
#include <chrono>
#include <stdexcept>
#include <algorithm>
#include <iomanip>
#include <cmath>
#pragma once
using namespace std;
using namespace chrono;
const double EPS = 1e-10; // Определите константу для сравнения с нулем
typedef std::vector<double> DoubleVector;
typedef std::vector<DoubleVector> DoubleVector2D;
void generateMatrix(int n, DoubleVector2D& A) // создание матрицы А
{
	for (int i = 1; i <= n; ++i)
	{
		DoubleVector Arow;
		for (int j = 1; j <= n; ++j)
		{
			double aVal;
			aVal = 1.0 / (2 + i + 3 * j);
			Arow.push_back(aVal);
		}

		A.push_back(Arow);
		Arow.clear();
	}
}

void generateRightHandSideVector(int n, DoubleVector2D& A, DoubleVector& b, DoubleVector& f) // создание вектора правой части f
{
	for (int i = 0; i < n; ++i) {
		double temp = 0.0;
		for (int j = 0; j < n; ++j) {
			temp += A[i][j] * b[i];
		}
		f.push_back(temp);
	}
}

class Matrix {
public:
	int M, N; // Размеры матрицы
	std::vector<std::vector<double>> Elem;

	Matrix(int m, int n) : M(m), N(n), Elem(m, std::vector<double>(n, 0)) {}

	void Copy(const Matrix& A) {
		for (int i = 0; i < M; ++i) {
			for (int j = 0; j < N; ++j) {
				Elem[i][j] = A.Elem[i][j];
			}
		}
	}
};

class Vector {
public:
	int size;
	std::vector<double> Elem;

	Vector(int s) : size(s), Elem(s, 0) {}
};

void Householder_Orthogonalization(Matrix& A, Matrix& Q, Matrix& R) {
	Vector p(A.M);

	double s, beta, mu;


	for (int i = 0; i < R.N - 1; i++) {
		// Находим квадрат нормы столбца для обнуления
		s = 0;
		for (int I = i; I < R.M; ++I) {
			s += R.Elem[I][i] * R.Elem[I][i];
		}

		if (std::sqrt(std::abs(s - R.Elem[i][i] * R.Elem[i][i])) > EPS) {
			// Выбор знака
			if (R.Elem[i][i] < 0) {
				beta = std::sqrt(s);
			}
			else {
				beta = -std::sqrt(s);
			}

			mu = 1.0 / (beta * (beta - R.Elem[i][i]));

			// Формирование вектора p
			for (int I = 0; I < R.M; ++I) {
				p.Elem[I] = (I >= i) ? R.Elem[I][i] : 0.0;
			}

			p.Elem[i] -= beta;

			// Вычисление новых компонент матрицы R
			for (int m = i; m < R.N; ++m) {
				// Произведение S = At * p
				s = 0;
				for (int n = i; n < R.M; ++n) {
					s += R.Elem[n][m] * p.Elem[n];
				}
				s *= mu;

				for (int n = i; n < R.M; ++n) {
					R.Elem[n][m] -= s * p.Elem[n];
				}
			}

			// Вычисление новых компонент матрицы Q
			for (int m = 0; m < R.M; ++m) {
				// Произведение Q * p
				s = 0;
				for (int n = i; n < R.M; ++n) {
					s += Q.Elem[m][n] * p.Elem[n];
				}
				s *= mu;

				for (int n = i; n < R.M; ++n) {
					Q.Elem[m][n] -= s * p.Elem[n];
				}
			}
		}
	}
}
typedef std::vector<double> DoubleVector;
typedef std::vector<DoubleVector> DoubleVector2D;


void generateRightHandSideVector(int n, Matrix& A, DoubleVector& b, DoubleVector& f) // создание вектора правой части f
{
	for (int i = 0; i < n; ++i) {
		double temp = 0.0;
		for (int j = 0; j < n; ++j) {
			temp += A.Elem[i][j] * b[i];
		}
		f.push_back(temp);
	}
}
Matrix GenerateMatrix(int N) {
	Matrix A(N, N);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			A.Elem[i][j] = 1.0 / (2 + (i + 1) + 3 * (j + 1));
		}
	}
	return A;
}

void Back_Row_Substitution(const Matrix& R, Vector& X, const Vector& B) {
	for (int i = R.M - 1; i >= 0; i--) {
		X.Elem[i] = B.Elem[i];
		for (int j = i + 1; j < R.N; j++)
			X.Elem[i] -= R.Elem[i][j] * X.Elem[j];
		X.Elem[i] /= R.Elem[i][i];
	}
}
int main() {
	setlocale(LC_ALL, "Ru");
	vector<int> sizes = { 5,10,20 };
	cout << setw(50);
	cout << "Метод Гаусса (част.выбор) || " << "QR-Хаусхолдер || " << " SVD-разложение" << endl;
	for (int n : sizes) {
		int h = 0;
		int k = 0;
		DoubleVector2D A;
		generateMatrix(n, A);

		DoubleVector b;
		for (int i = 0; i < n; ++i)
		{
			double bVal = 1.0;
			b.push_back(bVal);
		}
		DoubleVector f;
		generateRightHandSideVector(n, A, b, f);
		auto start = steady_clock::now();
		while ((h < n) && (k < n))
		{
			// нахождение пивотирующего элемента
			int iMax = h;
			double AcolMax = std::fabs(A[iMax][k]);
			for (int i = (h + 1); i < n; ++i)
			{
				if (std::fabs(A[i][k]) > AcolMax)
				{
					AcolMax = std::fabs(A[i][k]);
					iMax = i;
				}
			}

			if (A[iMax][k] == 0.0)
			{
				// если нет пивотирующего элемента, переходим к следующему столбцу
				k++;
			}
			else
			{
				// перестановка строк
				for (int j = 0; j < n; ++j)
				{
					double temp = A[h][j];
					A[h][j] = A[iMax][j];
					A[iMax][j] = temp;
				}
				double temp = f[h];
				f[h] = f[iMax];
				f[iMax] = temp;


				for (int i = (h + 1); i < n; ++i)
				{
					double m = A[i][k] / A[h][k];

					// явное обнуление обнуляемых компонент
					A[i][k] = 0.0;

					for (int j = (k + 1); j < n; ++j)
					{
						A[i][j] -= (A[h][j] * m);
					}
					f[i] -= (f[h] * m);
				}
				h++;
				k++;
			}
		}

		// обратная замена
		DoubleVector x_gauss;
		for (int j = 0; j < n; ++j)
		{
			x_gauss.push_back(0.0);
		}
		for (int i = (n - 1); i >= 0; --i)
		{
			double s = 0.0;
			for (int j = (i + 1); j < n; ++j)
			{
				s += A[i][j] * x_gauss[j];
			}
			x_gauss[i] = (f[i] - s) / A[i][i];
		}
		auto end = steady_clock::now();
		auto total_s = duration_cast<microseconds>(end - start);

		Matrix A_hous = GenerateMatrix(n);
		DoubleVector b_Hous;
		Matrix Q(n, n);
		Matrix R(n, n);
		Vector b_(n);
		DoubleVector f_hous;
		Vector result_hous = (n);
		for (int i = 0; i < n; i++) {
			result_hous.Elem[i] = 0;
		}
		generateRightHandSideVector(n, A_hous, b, f_hous);
		for (int i = 0; i < n; i++) {
			b_.Elem[i] = f_hous[i];
		}
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				Q.Elem[i][j] = (i == j) ? 1.0 : 0.0;
			}
		}
		R.Copy(A_hous);
		auto start1 = steady_clock::now();
		Householder_Orthogonalization(A_hous, Q, R);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				result_hous.Elem[i] += Q.Elem[j][i] * b_.Elem[j];
			}
		}
		Back_Row_Substitution(R, result_hous, result_hous);
		auto end1 = steady_clock::now();
		auto total_s1 = duration_cast<microseconds>(end1 - start1);
		Vector x_hous = result_hous;
		double residual = 0.0;
		double x_norm = 0.0;
		for (int i = 0; i < n; i++) {
			residual += (x_hous.Elem[i] - b[i]) * (x_hous.Elem[i] - b[i]);
			x_norm += b[i] * b[i];
		}
		double upper_gauss = 0.0;
		double lower_gauss = 0.0;
		for (int i = 0; i < n; i++) {
			upper_gauss += (x_gauss[i] - b[i]) * (x_gauss[i] - b[i]);
			lower_gauss += b[i] * b[i];
		}
		cout << "N: " << n << endl;
		cout << "Погрешность решения " << setw(20) << sqrt(upper_gauss) / sqrt(lower_gauss) << setw(22) << sqrt(residual) / sqrt(x_norm);
		cout << endl;
		cout << "Время решения (в секундах) " << setw(13) << total_s.count() * 1e-6 << setw(22) << total_s1.count() * 1e-6;
		cout << endl;
		cout << "---------------" << endl;
	}
}
