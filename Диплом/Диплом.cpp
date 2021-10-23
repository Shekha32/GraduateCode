// Диплом.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//


//Программа реализует метод решения дифференциального уравнения математической физики (уравнение переноса тепла или теплопроводности) с помощью численного
//метода сеток. Уравнение является двумерным (рассматриваемая область - пластина) и динамическим (нестационарным, зависящим от времени).
//В качестве входных данных вводятся параметры для описания граничных условий, геометрические размеры пластины, шаги по сеткам.
//Данная программа позволяет заменить дорогостоящий физический эксперимент вычислительным. 
//Таким образом, изменив исходные данные, программа дает возможность моделировать всевозможные тепловые процессы в пластине.
//В программе реализованы два класса (Vector и Matrix) и связь между ними - агрегация.


//#include "pch.h"
#include <iostream>
#include <stdlib.h>
#include <cstdlib>
#include <conio.h>
#include <stdio.h>
#include <iomanip>
#include <math.h>
#include <string.h>
#include <locale.h>
#include <time.h>
#include <vector>
#include <algorithm>
#include <functional>
#include <fstream>
using namespace std;

typedef unsigned int ui;
using namespace std;
const double PI = 3.141592653589793;
const double eps = 0.0000001;

//поскольку схема неявная, то она абсолютно устойчива = > можно брать любые значения шагов по x, y, t
const double hx = 0.1;            //шаг по x
const double hy = 0.05;           //шаг по y
const double tau = 1.0;           //шаг по t

const double lx = 3.0;            //длина пластины (x)   [м]
const double ly = 1.0;            //ширина пластины (y)  [м]

const int K1 = (ly / hy) + 1.0;   //количество узлов по x
const int K2 = ly / hy;           //количество отрезков hx по x

const int M1 = (lx / hx) + 1.0;   //количество узлов по y
const int M2 = lx / hx;           //количество отрезков hy по y

int q = 1500;                     //воздействие теплового потока (градиент температуры)  [Вт/(м^2)]
int alpha = 10;                   //коэффициент теплоотдачи от стенки к среде  [Вт/(град*м^2)]                  
int Tc = 15;                      //температура окружающей среды  [град]


class Vector
{
private:
	ui size;
	long double *mas;

public:

	Vector()  //конструктор по умолчанию
	{
		size = 0;
		mas = NULL;
	}


	Vector(const Vector &obj)  //копирующий конструктор
	{
		size = obj.size;
		mas = new long double[size];
		for (int i(0); i < size; i++)
			mas[i] = obj[i];
	}


	Vector(ui n, const int k = 0)  //конструктор, принимающий размерность массива и заполняющий его нулями
	{
		size = n;
		mas = new long double[size];

		for (int i(0); i < size; i++)
			mas[i] = k;
	}


	~Vector() { if (mas != NULL) delete[] mas; }  //деструктор


	void Vvod1(ui n1, const int k = 0)  //метод для ввода (создания) вектора
	{
		if (mas != NULL)
		{
			delete[] mas;
			size = 0;
		}

		size = n1;

		mas = new long double[size];
		for (int i(0); i < size; i++)
			mas[i] = k;
	}


	void Vvod2(ui n1, long double *Y)  //метод для заполнения вектора в СЛАУ
	{
		if (mas != NULL)
		{
			delete[] mas;
			size = 0;
		}

		size = n1;

		mas = new long double[size];
		for (int i(0); i < size; i++)
			mas[i] = Y[i];
	}


	Vector &operator=(const Vector &obj)  //оператор присваивания = (приравнивает один вектор к другому)
	{
		if (mas != NULL)
		{
			delete[] mas;
			size = 0;
		}

		size = obj.size;
		if (size > 0) mas = new long double[size];
		else mas = NULL;

		for (int i(0); i < size; i++)
			mas[i] = obj[i];

		return *this;
	}


	long double &operator[](ui ind) const  //оператор индексирования [] (возвращает элемент с введенным (переданным) индексом)
	{
		if (ind < 0) cout << "Введенный индекс меньше нуля";
		if (ind > size) cout <<"Введенный индекс больше размерности вектора";
		return mas[ind];
	}


	long double &operator()(ui ind) const  //оператор индексирования () (возвращает элемент с введенным (переданным) индексом)
	{
		if (ind < 0) cout << "Введенный индекс меньше нуля";
		if (ind > size) cout << "Введенный индекс больше размерности вектора";
		return mas[ind];
	}


	Vector operator-()  //унарный минус (умножает каждый элемент вектора на -1) 
	{
		for (int i(0); i < size; i++)
			mas[i] = -mas[i];

		return *this;
	}


	Vector operator+(Vector &obj)  //бинарная операция сложения векторов + (складывает два вектора любых размеров)
	{
		ui n;

		if (size <= obj.getSize()) n = obj.size;
		else n = size;

		Vector temp(n);

		for (int i(0); i < size; i++)
			temp[i] = mas[i];

		for (int i(0); i < obj.size; i++)
			temp[i] += obj[i];

		return temp;
	}


	Vector operator-(Vector &obj)  //бинарная операция разность векторов - (поэлементная разность векторов)
	{
		ui n;

		if (size <= obj.getSize()) n = obj.size;
		else n = size;

		Vector temp(n);

		for (int i(0); i < size; i++)
			temp[i] = mas[i];

		for (int i(0); i < obj.size; i++)
			temp[i] -= obj[i];

		return temp;
	}


	long double operator*(Vector &b)  //бинарная операция умножения * (скалярное произведение векторов)
	{
		ui n;
		long double sum = 0;
		
		if (size <= b.size) n = size;
		else n = b.size;

		for (int i(0); i < n; i++)
			sum += mas[i] * b.mas[i];

		return sum;
	}


	Vector operator*(long double b)  //бинарная операция умножения * (умножает вектор на введенное (переданное) число)
	{
		Vector temp(size);

		for (int i(0); i < size; i++)
			temp[i] = mas[i] * b;

		return temp;
	}

	//friend Vector operator*(float a, Vector &b);
	friend Vector operator*(long double a, Vector &b);
	friend ostream &operator<<(ostream &os, Vector obj);
	ui getSize() const { return size; }  //метод, возвращающий текущий размер массива (геттер)
};


Vector operator*(long double a, Vector &b)  //бинарная операция умножения * (умножает введенное (переданное) число на вектор)
{
	Vector temp(b.getSize());

	for (int i(0); i < b.getSize(); i++)
		temp[i] = b[i] * a;

	return temp;
}


ostream &operator<<(ostream &os, Vector obj)  //оператор вывода в поток << (печать всего вектора)
{
	for (int i(0); i < obj.getSize(); i++)
		printf("%1.3lf\n", obj[i]);

	return os;
};



class Matrix
{
protected:
	ui n;         //количество строк
	ui m;         //количество столбцов
	Vector *mas;  //агрегация

public:
	Matrix()  //конструктор по умолчанию
	{
		n = m = 0;
		mas = NULL;
	}


	Matrix(const Matrix &obj)  //копирующий конструктор
	{
		n = obj.getN();
		m = obj.getM();
		mas = new Vector[n];

		for (int i(0); i < n; i++)
		{
			mas[i].Vvod1(m);

			for (int j = 0; j < m; j++)
				mas[i][j] = obj[i][j];
		}
	}


	Matrix(ui n1, ui m1, const int k = 0)  //конструктор, принимающий размерности матрицы
	{
		n = n1;
		m = m1;
		mas = new Vector[n];

		for (int i(0); i < n; i++)
			mas[i].Vvod1(m);
	}


	~Matrix() { if (mas != NULL) delete[] mas; }  //деструктор


	Matrix &operator=(const Matrix &obj)  //оператор присваивания = (приравнивает одну матрицу к другой)
	{
		if (mas != NULL)
		{
			n = m = 0;
			delete[] mas;
		}

		n = obj.getN();
		m = obj.getM();
		mas = new Vector[n];

		for (int i(0); i < n; i++)
		{
			mas[i].Vvod1(m);

			for (int j(0); j < m; j++)
				mas[i][j] = obj[i][j];
		}

		return *this;
	}


	Vector &operator[](ui ind) const  //оператор индексирования [] (возвращает элемент с введенным (переданным) индексом)
	{
		return mas[ind];
	}


	Vector &operator()(ui ind)  //оператор индексирования () (возвращает элемент с введенным (переданным) индексом)
	{
		return mas[ind];
	}


	Matrix operator+(Matrix &obj)  //бинарная операция сложения матриц + (складывает две матрицы любых размеров согласно правилу)
	{
		if (m == obj.getM())
		{
			ui n1, m1;
			if (n >= obj.getN()) n1 = n;
			else n1 = obj.getN();
			if (m >= obj.getM()) m1 = m;
			else m1 = obj.getM();

			Matrix temp(n1, m1);

			for (int i(0); i < n; i++)
				temp[i] = mas[i];

			for (int i(0); i < obj.getN(); i++)
				temp[i] = temp[i] + obj[i];

			return temp;
		}
		else cout << "\nОперация невозможна! Размеры матриц не совпадают\n";
	}


	Matrix operator-()  //унарный минус (умножает каждый элемент матрицы на -1)
	{
		for (int i(0); i < n; i++)
			for (int j(0); j < m; j++)
				mas[i][j] = -mas[i][j];

		return *this;
	}


	Matrix operator-(Matrix &obj)  //бинарная операция разности матриц - (поэлементная разность матриц согласно правилу)
	{
		if (m == obj.getN())
		{
			ui n1, m1;
			if (n >= obj.getN()) n1 = n;
			else n1 = obj.getN();
			if (m >= obj.getM()) m1 = m;
			else m1 = obj.getM();

			Matrix temp(n1, m1);

			for (int i(0); i < n; i++)
				temp[i] = mas[i];

			for (int i(0); i < obj.getN(); i++)
				temp[i] = temp[i] - obj[i];

			return temp;
		}
		else cout << "\nОперация невозможна! Размеры матриц не совпадают\n";
	}


	Matrix operator*(Matrix &obj)  //бинарная операция умножения * (перемножает две матрицы согласно правилу)
	{
		if (m == obj.getN())
		{
			ui n1, m1;
			n1 = n;
			m1 = obj.getM();
			Matrix temp(n1, m1);

			for (int i(0); i < n1; i++)
				for (int j(0); j < m1; j++)
				{
					temp[i][j] = 0;
					for (int k(0); k < m1; k++)
						temp[i][j] = temp[i][j] + mas[i][k] * obj[k][j];
				}

			return temp;
		}

		else cout << "\nОперация невозможна! Размеры матриц не совпадают\n";
	}

	
	Vector operator*(Vector &obj)  //бинарная операция умножения * (умножение матрицы на вектор)
	{
		ui n1, m1;
		n1 = n;
		m1 = obj.getSize();
		Vector temp(n1);

		for (int i(0); i < n1; i++)
		{
			temp[i] = 0;

			for (int j(0); j < m1; j++)
				temp[i] = temp[i] + mas[i][j] * obj[j];
		}

		return temp;
	}


	Matrix operator*(long double b)  //бинарная операция умножения * (умножает матрицу (поэлементно) на введенное (переданное) число)
	{
		Matrix temp(n, m);

		for (int i(0); i < n; i++)
			for (int j(0); j < m; j++)
				temp[i][j] = mas[i][j] * b;

		return temp;
	}


	Matrix Transposition()  //метод транспонирования матрицы
	{
		Matrix temp(m, n);

		for (int i(0); i < m; i++)
			for (int j(0); j < n; j++)
				temp[i][j] = mas[j][i];
		
		return temp;
	}


	long double Determinant()  //вычисление детерминанта. матрица должна быть квадратной, тогда берем количество строк
	{
		Matrix A(*this);

		if (n == 1) return mas[0][0];
		if (n == 2) return mas[0][0] * mas[1][1] - mas[0][1] * mas[1][0];

		else
		{
			int k = 0;
			while (mas[0][k] == 0) k++;
			if (k == n) return 0;
			else
				return pow(-1, k)*A.GausStep(0, k).mas[0][k] * (A.GausStep(0, k).Minor(0, k)).Determinant();
		}
	}


	Matrix Minor(const int p, const int q)  //поиск минора матрицы
	{
		Matrix minor(n - 1, m - 1);

		for (int i = 0; i < p; i++)
		{
			for (int j = 0; j < q; j++)
				minor[i][j] = mas[i][j];

			for (int j = q; j < n - 1; j++)
				minor[i][j] = mas[i][j + 1];
		}

		for (int i = p; i < n - 1; i++)
		{
			for (int j = q; j < n - 1; j++)
				minor[i][j] = mas[i + 1][j + 1];

			for (int j = 0; j < q; j++)
				minor[i][j] = mas[i + 1][j];
		}

		Matrix result(minor);
		return result;
	}


	Matrix GausStep(int p, int q)  //функция для СЛАУ
	{
		Matrix result(n, m);

		for (int i(0); i < n; i++)
			for (int j(0); j < m; j++)
				result[i][j] = mas[i][j];

		for (int i(0); i < n; i++)
			for (int j(0); j < m; j++)
			{
				if (i != p)
					result[i][j] = mas[i][j] - (mas[p][j] * mas[i][q] / mas[p][q]);
			}

		Matrix A(result);
		return A;
	}


	Matrix ReverseMatrix()  //поиск обратной матрицы
	{
		Matrix A(n, m);

		if (Determinant() != 0)
		{
			long double det = Determinant();

			for (int i(0); i < n; i++)
				for (int j(0); j < m; j++)
					A.mas[i][j] = Minor(i, j).Determinant()*pow(-1, i + j) / det;
		}
		
		return A.Transposition();
	}


	void Vvod3(float **A, ui n1) //для заполнения матрицы в СЛАУ
	{
		n = n1;
		m = n1;
		mas = new Vector[n];

		for (int i(0); i < n; i++)
		{
			mas[i].Vvod1(m);

			for (int j(0); j < m; j++)
				mas[i][j] = A[i][j];
		}
	}


	Matrix Distribution(double t, Vector lambda_x, Vector lambda_y, Vector a_x, Vector a_y, double heat, double density, int &num, char f_)  //функция, реализующая метод переменных направлений при введенном времени t
	{
		//t - введенное пользователем время
		//lambda - коэффициент теплопроводности вещества
		//a - коэффициент температуропроводности
		//num - счетчик итераций по времени

		Matrix T(K1, M1);             //основная матрица распределения температуры для введенного времени
		Matrix T_(K1, M1);            //вспомогательная матрица для промежуточного слоя n+1/2
		Matrix Lx(K1 - 2, M1 - 1);    //матрица коэффициентов для неявного шага по x
		Matrix Mx(K1 - 2, M1 - 1);    //матрица коэффициентов для неявного шага по x
		Matrix Ly(K1 - 1, M1 - 2);    //матрица коэффициентов для неявного шага по y
		Matrix My(K1 - 1, M1 - 2);    //матрица коэффициентов для неявного шага по y
		Matrix f(K1, M1);             //функция плотности внешних потоков

		//введем константы для функции f (плотность тепловых источников):
		double A = 0.4 * heat * density;
		double psi = 1;
		
		long double alpha_x;          //замена
		long double beta_x;           //замена
		long double alpha_y;          //замена
		long double beta_y;           //замена

		for (int i(0); i < K1; i++)
		{
			for (int j(0); j < M1; j++)
			{
				T[i][j] = 10;         //начальная температура пластины в нулевой момент времени (начальное условие)
				T_[i][j] = 10;        //начальная температура пластины в нулевой момент времени (начальное условие)
			}
		}


		for(int k(1); tau * k <= t; k++)  //k - счетчик итераций по времени        
		{
			//1 этап. Вычисление температуры в пластине на n+1/2 слое
			alpha_x = -tau / (2 * pow(hx, 2));    //замена
			beta_x = tau / (pow(hx, 2));          //замена
			alpha_y = tau / (2 * pow(hy, 2));     //замена

			if (f_ == '1')  //если есть f, то вычисляем функцию f
			{
				for (double i(ly), i1(0); i >= -eps; i -= hy, i1++)  //копится погрешность, поэтому добавляем eps
				{
					for (double j(0), j1(0); j <= lx + eps; j += hx, j1++)  //копится погрешность, поэтому добавляем eps
					{
						//берем [tau*k - tau/2] из-за n+1/2 слоя
						f[i1][j1] = (A * ((tau * k) - (tau / 2)) * exp((-psi * ((tau * k) - (tau / 2))) / 10) * sin(PI * j * ((tau * k) - (tau / 2)) / 4) * cos(PI * i * ((tau * k) - (tau / 2)) / 2)) / (heat * density);  //вычисление функции f на n+1/2 слое
					}
				}
			}

			for (int i(0); i <= K1 - 3; i++)
			{
				Lx[i][0] = 1;
				Mx[i][0] = ((q + 25 * sqrt((tau * k) - (tau / 2))) * hx) / lambda_x[1];
			}

			//прямой ход прогонки - вычисление коэффициентов Lx, Mx
			for (int i(0); i <= K1 - 3; i++)
			{
				for (int j(1); j <= M1 - 2; j++)
				{
					long double Fx = T[i + 1][j + 1] + a_y[i + 1] * alpha_y * (T[i][j + 1] - 2 * T[i + 1][j + 1] + T[i + 2][j + 1]);  //вычисление функции с предыдущего слоя

					Lx[i][j] = (-a_x[j + 1] * alpha_x) / ((1 + a_x[j + 1] * beta_x) + a_x[j + 1] * alpha_x * Lx[i][j - 1]);
					Mx[i][j] = (Fx - a_x[j + 1] * alpha_x * Mx[i][j - 1]) / ((1 + a_x[j + 1] * beta_x) + a_x[j + 1] * alpha_x * Lx[i][j - 1]);
				}
			}

			//обратный ход прогонки - вычисление матрицы T_ на n+1/2 слое
			for (int i(1); i <= K1 - 2; i++)
			{
				for (int j(M1 - 1); j >= 0; j--)
				{
					if (j == M1 - 1)
						T_[i][j] = (alpha * hx * Tc + lambda_x[j] * Mx[i - 1][j - 1]) / (lambda_x[j] + alpha * hx - lambda_x[j] * Lx[i - 1][j - 1]);
				
					else
						T_[i][j] = Lx[i - 1][j] * T_[i][j + 1] + Mx[i - 1][j];
				}
			}

			//из граничных условий по y вычисляем 0 и K1-1 строки
			for (int j(0), i1(0), i2(K1 - 1); j <= M1 - 1; j++)
			{
				T_[i1][j] = T_[i1 + 1][j];
				T_[i2][j] = T_[i2 - 1][j];
			}

					   
			//2 этап. Вычисление температуры в пластине на n+1 слое
			alpha_y = -tau / (2 * pow(hy, 2));
			beta_y = tau / (pow(hy, 2));
			alpha_x = tau / (2 * pow(hx, 2));

			for (int j(0); j <= M1 - 3; j++)
			{
				Ly[K1 - 2][j] = 1;
				My[K1 - 2][j] = 0; //можно не писать (0 по умолчанию)
			}

			//прямой ход прогонки - вычисление коэффициентов Ly, My
			for (int j(0); j <= M1 - 3; j++)
			{
				for (int i(K1 - 3); i >= 0; i--)
				{
					long double Fy = (T_[i][j + 1]) + a_x[j] * alpha_x * (T_[i][j + 2] - 2 * T_[i][j + 1] + T_[i][j]);  //вычисление функции с предыдущего слоя

					Ly[i][j] = (-a_y[i] * alpha_y) / ((1 + a_y[i] * beta_y) + a_y[i] * alpha_y * Ly[i + 1][j]);
					My[i][j] = (Fy - a_y[i] * alpha_y * My[i + 1][j]) / ((1 + a_y[i] * beta_y) + a_y[i] * alpha_y * Ly[i + 1][j]);
				}
			}

			//обратный ход прогонки - вычисление матрицы T на n+1 слое
			for (int j(1); j <= M1 - 2; j++)
			{
				for (int i(0); i <= K1 - 1; i++)
				{
					if (i == 0)
						T[i][j] = My[i][j - 1] / (1 - Ly[i][j - 1]);

					else
						T[i][j] = Ly[i - 1][j - 1] * T[i - 1][j] + My[i - 1][j - 1];
				}
			}


			//из граничных условий по x вычисляем 0 и M1-1 столбцы 
			for (int i(0), j1(0), j2(M1 - 1); i <= K1 - 1; i++)
			{ 
				T[i][j1] = T[i][j1 + 1] + ((q + 25 * sqrt(tau * k)) * hx / lambda_x[j1]);
				T[i][j2] = (lambda_x[j2] / (lambda_x[j2] + alpha * hx)) * T[i][j2 - 1] + ((alpha * hx * Tc) / (lambda_x[j2] + alpha * hx));
			}

			T = T + f;    //учитываем функцию f на полуцелом слое

			T_ = T;       //приравниваем матрицы для следующей итерации
			num = k;      //счетчик итераций
		}

		return T;
	}


	void Write_in(ofstream &f, Matrix T)  //функция ввода данных в файл
	{
		for (int i(0); i < K1; i++)
		{
			for (int j(0); j < M1; j++)
				f << T[i][j] << "  ";
			
			f << endl;
		}
	}


	friend Matrix operator*(long double a, Matrix &b);
	friend Matrix operator*(Vector &v, Matrix &m);
	friend ostream &operator<<(ostream &os, Matrix obj);
	ui getN() const { return n; } //метод, возвращающий количество строк
	ui getM() const { return m; } //метод, возвращающий количество столбцов
};


Matrix operator*(long double a, Matrix &b) //бинарная операция умножения * (умножает (поэлементно) введенное (переданное) число на матрицу)
{
	Matrix temp(b.getN(), b.getM());

	for (int i(0); i < b.getN(); i++)
		for (int j(0); j < b.getM(); j++)
			temp[i][j] = b[i][j] * a;

	return temp;
}


Matrix operator*(Vector &V, Matrix &M) //бинарная операция умножения * (умножение вектора на матрицу)
{
	ui n1, m1;
	n1 = M.getN();
	m1 = V.getSize();
	Matrix temp(n1, 1);

	for (int i(0); i < n1; i++)
	{
		temp[i][0] = 0;

		for (int j(0); j < m1; j++)
			temp[i][0] = temp[i][0] + M.mas[i][j] * V[j];
	}

	return temp;
}


ostream &operator<<(ostream &os, Matrix obj) //оператор вывода в поток << (печать всей матрицы)
{
	for (int i(0); i < obj.getN(); i++)
	{
		for (int j(0); j < obj.getM(); j++)
			printf("%1.3lf   ", obj[i][j]);

		os << endl;
	}

	return os;
}


long double Recovery1(double t, long double Ti, Vector lambda_x, Vector lambda_y, Vector a_x, Vector a_y, Vector &lam_x, Vector &lam_y, Vector &a_xx, Vector &a_yy, double heat, double density, char f, char m, int &k)  //метод восстановления для констант
{
	Matrix T(K1, M1);
	const double eps = 0.00000001;    //погрешность для минимизации разности
	const double delta1 = 0.01;       //шаг для коэффициентов теплопроводности

	int n = 0;                        //нумерация итерация для функции распределения Distribution
	int a1 = 70, b1 = 200;            //допустимые пределы для коэффициентов теплопроводности

	for (int j(0); j < M1; j++)
	{
		lam_x[j] = a1;
		a_xx[j] = 10 * lam_x[j] / (heat * density);
	}

	for (int i(0); i < K1; i++)
	{
		lam_y[i] = a1;
		a_yy[i] = 10 * lam_y[i] / (heat * density);
	}

	if (lambda_x[0] > lambda_y[0])
	{
		while (pow(Ti - T.Distribution(t, lam_x, lam_y, a_xx, a_yy, heat, density, n, f)[12][2], 2) > eps && abs(lambda_x[0] - lam_x[0]) > delta1)
		{
			k++;

			if (abs(lambda_x[0] - lam_x[0]) > delta1 && lam_x[0] <= b1)
			{
				for (int j(0); j < M1; j++)
				{
					lam_x[j] += delta1;
					a_xx[j] = 10 * lam_x[j] / (heat * density);
				}
			}

			if (abs(lambda_y[0] - lam_y[0]) > delta1 && lam_y[0] <= b1)
			{
				for (int i(0); i < K1; i++)
				{
					lam_y[i] += delta1;
					a_yy[i] = 10 * lam_y[i] / (heat * density);
				}
			}
		}
	}

	else
	{		
		while (pow(Ti - T.Distribution(t, lam_x, lam_y, a_xx, a_yy, heat, density, n, f)[12][2], 2) > eps && abs(lambda_y[0] - lam_y[0]) > delta1)
		{
			k++;

			if (abs(lambda_x[0] - lam_x[0]) > delta1 && lam_x[0] <= b1)
			{
				for (int j(0); j < M1; j++)
				{
					lam_x[j] += delta1;
					a_xx[j] = 10 * lam_x[j] / (heat * density);
				}
			}

			if (abs(lambda_y[0] - lam_y[0]) > delta1 && lam_y[0] <= b1)
			{
				for (int i(0); i < K1; i++)
				{
					lam_y[i] += delta1;
					a_yy[i] = 10 * lam_y[i] / (heat * density);
				}
			}
		}
	}

	return T.Distribution(t, lam_x, lam_y, a_xx, a_yy, heat, density, n, f)[12][2]; //возвращает восстановленное значение функции в заданной точке
}


void Recovery2(double t, Vector lambda_x, Vector lambda_y, Vector a_x, Vector a_y, Vector &lam_x, Vector &lam_y, Vector &a_xx, Vector &a_yy, double heat, double density, char f, char m, int &k)  //метод восстановления для функций
{
	Matrix T(K1, M1);
	const double eps = 0.000000000001;  //погрешность для минимизации разности
	const double delta1 = 0.06;         //шаг для коэффициентов теплопроводности
	const double delta2 = 0.04;         //погрешность для коэффициентов теплопроводности

	int n = 0;                          //нумерация итерация для функции распределения Distribution
	int a1 = 80, b1 = 200;              //допустимые пределы для коэффициентов теплопроводности

	for (int j(0); j < M1; j++)
	{
		lam_x[j] = a1;
		a_xx[j] = 10 * lam_x[j] / (heat * density);
	}

	for (int i(0); i < K1; i++)
	{
		lam_y[i] = a1;
		a_yy[i] = 10 * lam_y[i] / (heat * density);
	}

	//отсюда находим только lam_y[0]
	for (int j(0), i(0); j < 1; j++)    
	{
		if (lambda_x[j] > lambda_y[i])
		{
			while (pow(T.Distribution(t, lambda_x, lambda_y, a_x, a_y, heat, density, n, f)[i][j] - T.Distribution(t, lam_x, lam_y, a_xx, a_yy, heat, density, n, f)[i][j], 2) > eps && abs(lambda_x[j] - lam_x[j]) > delta2)  //специальное условие для констант
			{
				k++;

				if (abs(lambda_x[j] - lam_x[j]) > delta2 && lam_x[j] <= b1)
				{
					lam_x[j] += delta1;
					a_xx[j] = 10 * lam_x[j] / (heat * density);
				}

				if (abs(lambda_y[i] - lam_y[i]) > delta2 && lam_y[i] <= b1)
				{
					lam_y[i] += delta1;
					a_yy[i] = 10 * lam_y[i] / (heat * density);
				}
			}
		}

		else
		{
			while (pow(T.Distribution(t, lambda_x, lambda_y, a_x, a_y, heat, density, n, f)[i][j] - T.Distribution(t, lam_x, lam_y, a_xx, a_yy, heat, density, n, f)[i][j], 2) > eps && abs(lambda_y[i] - lam_y[i]) > delta2)  //специальное условие для констант
			{
				k++;

				if (abs(lambda_x[j] - lam_x[j]) > delta2 && lam_x[j] <= b1)
				{
					lam_x[j] += delta1;
					a_xx[j] = 10 * lam_x[j] / (heat * density);
				}

				if (abs(lambda_y[i] - lam_y[i]) > delta2 && lam_y[i] <= b1)
				{
					lam_y[i] += delta1;
					a_yy[i] = 10 * lam_y[i] / (heat * density);
				}
			}
		}
	}

	//восстановление вектора lam_x
	for (int j(0), i(0); j < M1; j++)
	{
		while (pow(T.Distribution(t, lambda_x, lambda_y, a_x, a_y, heat, density, n, f)[i][j] - T.Distribution(t, lam_x, lam_y, a_xx, a_yy, heat, density, n, f)[i][j], 2) > eps && abs(lambda_x[j] - lam_x[j]) > delta2)  //специальное условие для констант
		{
			k++;

			if (abs(lambda_x[j] - lam_x[j]) > delta2 && lam_x[j] <= b1)
			{
				lam_x[j] += delta1;
				a_xx[j] = 10 * lam_x[j] / (heat * density);
			}
		}
	}
	
	//восстановление вектора lam_y
	for (int i(0), j(0); i < K1; i++)
	{
		while (pow(T.Distribution(t, lambda_x, lambda_y, a_x, a_y, heat, density, n, f)[i][j] - T.Distribution(t, lam_x, lam_y, a_xx, a_yy, heat, density, n, f)[i][j], 2) > eps && abs(lambda_y[i] - lam_y[i]) > delta2)  //специальное условие для констант
		{
			k++;

			if (abs(lambda_y[i] - lam_y[i]) > delta2 && lam_y[i] <= b1)
			{
				lam_y[i] += delta1;
				a_yy[i] = 10 * lam_y[i] / (heat * density);
			}
		}
	}
}



int main()
{
	//setlocale(0, "Russian");
	srand((unsigned)time(NULL));
	system("color 4B");            //цвет и фон консоли

	double t;                      //пользователь вводит время, для которого нужно произвести расчет (неотрицательная целочисленная переменная)
	int num = 0, k = 0;            //количество итераций по времени t и для обратной задачи
	char m, f;                     //переменные для меню и функции f

	Matrix T(K1, M1);              //матрица распределения температуры
	Vector lambda_x(M1), lam_x(M1);
	Vector lambda_y(K1), lam_y(K1);
	Vector a_x(M1), a_xx(M1);
	Vector a_y(K1), a_yy(K1);

	ofstream file("Matrix_Distribution.txt");  //создание нового файла в текущем каталоге

	//для обратной задачи задаем узловую точку и время
	int t_ = 3;
	int x_ = 2.6;
	int y_ = 1.2;

	cout << "Введите время t: ";   //пользователь вводит время, для которого необходимо произвести расчет
	cin >> t;
	cout << endl;

	if (t < 0)                     //проверка на правильность введения данных
	{
		cout << "Вы ввели отрицательное значение!\n\n";
		system("pause");
		return 0;
	}

	cout << "Добавить функцию внутренних источников теплоты f(x, y, t)?  [1 - да, иначе - нет]  ";
	cin >> f;
	cout << endl;

	cout << "1. Тест №1 - константы" << endl;
	cout << "2. Тест №2 - линейные функции" << endl;
	cout << "3. Тест №3 - тригонометрические функции" << endl << endl;
	cout << "Выберите пункт меню: ";
	cin >> m;

	switch (m)          //открытие меню
	{
	case '1':
	{
		system("cls");  //очистка экрана

		double conductivity_x = 80.4;                          //теплопроводность  [Вт/(м*град)]
		double conductivity_y = 95.7;                          //теплопроводность  [Вт/(м*град)]
		double heat = 0.460;                                   //удельная теплоемкость  [кДж/(кг*град)]
		double density = 7870;                                 //плотность  [кг/(м^3)]

		for (int j(0); j < M1; j++)
		{
			lambda_x[j] = conductivity_x;
			a_x[j] = 10 * lambda_x[j] / (heat * density);
		}

		for (int i(0); i < K1; i++)
		{
			lambda_y[i] = conductivity_y;
			a_y[i] = 10 * lambda_y[i] / (heat * density);
		}

		//прямая задача
		printf("\tТест №1 - константы\t\t\t\t\tВведенное время t = %.1f\n", t);
		printf("\tКоэффициенты теплопроводности lambda = %.1lf; %.1lf\tДлина пластины lx = %.0lf\n", conductivity_x, conductivity_y, lx);
		printf("\tКоэффициенты температуропроводности a = %.2lf; %.2lf\tШирина пластины ly = %.0lf\n", a_x[0], a_y[0], ly);
		printf("\tВоздействие теплового потока q = %d + 25*t\t\tКоэффициент теплоотдачи от стенки к среде alpha = %d\n", q, alpha);
		printf("\tТемпература окружающей среды Tc = %d\t\t\tШаг по времени t, tau = %.3lf\n", Tc, tau);
		printf("\tШаг по оси x, hx = %.3lf\t\t\t\tШаг по оси y, hy = %.3lf\n\n", hx, hy);

		if(f == '1')
			printf("\tФункция внутренних источников теплоты f(x, y, t) = A*t*EXP(-psi*t/10)*SIN(PI*x*t/4)*COS(PI*y*t/2)\n\n");
		else
			printf("\tФункция внутренних источников теплоты f(x, y, t) = 0\n\n");

		printf("\tРаспределение температуры в пластине (матрица T[%dx%d]):\n", K1, M1);
		cout << endl << T.Distribution(t, lambda_x, lambda_y, a_x, a_y, heat, density, num, f);
		cout << endl << endl << "\tКоличество итераций по времени t: " << num;		

		//записываем распределение в файл
		T.Write_in(file, T.Distribution(t, lambda_x, lambda_y, a_x, a_y, heat, density, num, f));
		printf("\n\n\n\n\n\n\n");

		//обратная задача
		//выбираем точку (0,4; 0,1) (для матрицы T это [12][2]), и для заданного времени t=3 находим значение функции T
		//значение функции в этой точке
		long double T0 = T.Distribution(t_, lambda_x, lambda_y, a_x, a_y, heat, density, num, f)[12][2];
		//вызов функции восстановления
		long double Tr = Recovery1(t_, T0, lambda_x, lambda_y, a_x, a_y, lam_x, lam_y, a_xx, a_yy, heat, density, f, m, k);

		cout << "\tРешение обратной задачи - восстановление теплофизических характеристик:\n\n";

		cout << "Заданное значение T(0,4; 0,1; 3) = " << T0 << endl;
		cout << "Восстановленное значение T_(0,4; 0,1; 3) = " << Tr << endl;
		cout << "Минимизация: " << abs(T0 - Tr) << "\nКоличество итераций: " << k << endl << endl;

		cout << "\tКоэффициент теплопроводности по x:\n";
		cout << "Заданное значение lambda_x = " << conductivity_x << endl;
		cout << "Восстановленное значение lam_x = " << lam_x[0] << endl;
		cout << "Разность по модулю: " << abs(conductivity_x - lam_x[0]) << endl << endl;

		cout << "\tКоэффициент теплопроводности по y:\n";
		cout << "Заданное значение lambda_y = " << conductivity_y << endl;
		cout << "Восстановленное значение lam_y = " << lam_y[0] << endl;
		cout << "Разность по модулю: " << abs(conductivity_y - lam_y[0]) << endl << endl;

		cout << "\tКоэффициент температуропроводности по x:\n";
		cout << "Заданное значение a_x = " << a_x[0] << endl;
		cout << "Восстановленное значение a_xx = " << a_xx[0] << endl;
		cout << "Разность по модулю: " << abs(a_x[0] - a_xx[0]) << endl << endl;

		cout << "\tКоэффициент температуропроводности по y:\n";
		cout << "Заданное значение a_y = " << a_y[0] << endl;
		cout << "Восстановленное значение a_yy = " << a_yy[0] << endl;
		cout << "Разность по модулю: " << abs(a_y[0] - a_yy[0]) << endl;

		cout << "\n\n\tИзмеренное распределение температуры при t=3:\n\n";
		cout << T.Distribution(t_, lambda_x, lambda_y, a_x, a_y, heat, density, num, f) << endl << endl;
		cout << "\tВосстановленное распределение температуры при t=3:\n\n";
		cout << T.Distribution(t_, lam_x, lam_y, a_xx, a_yy, heat, density, num, f);
		printf("\n\n\n\n\n\n\n");
		break;
	}

	case '2':
	{
		system("cls");  //очистка экрана

		double conductivity_x = 100;                          //теплопроводность  [Вт/(м*град)]
		double conductivity_y = 120;                          //теплопроводность  [Вт/(м*град)]
		double heat = 0.460;                                  //удельная теплоемкость  [кДж/(кг*град)]
		double density = 7870;                                //плотность  [кг/(м^3)]

		for (double j(0), j1(0); j < M1, j1 <= lx + eps; j++, j1 += hx)
		{
			lambda_x[j] = conductivity_x - 5 * j1;
			a_x[j] = 10 * lambda_x[j] / (heat * density);
		}

		for (double i(0), i1(0); i < K1, i1 <= ly + eps; i++, i1 += hy)
		{
			lambda_y[i] = conductivity_y + 15 * i1;
			a_y[i] = 10 * lambda_y[i] / (heat * density);
		}

		//прямая задача
		printf("\tТест №2 - линейные функции\t\t\t\t\t\t\t\t\t\tВведенное время t = %.1f\n", t);
		printf("\tКоэффициенты теплопроводности lambda_x = %.0lf-5x; lambda_y = %.0lf+15y\t\t\t\t\tДлина пластины lx = %.0lf\n", conductivity_x, conductivity_y, lx);
		printf("\tКоэффициенты температуропроводности a_x = lambda_x/(heat*density); a_y = lambda_y/(heat*density)\tШирина пластины ly = %.0lf\n", ly);
		printf("\tВоздействие теплового потока q = %d + 25*t\t\t\t\t\t\t\t\tКоэффициент теплоотдачи от стенки к среде alpha = %d\n", q, alpha);
		printf("\tТемпература окружающей среды Tc = %d\t\t\t\t\t\t\t\t\tШаг по времени t, tau = %.3lf\n", Tc, tau);
		printf("\tШаг по оси x, hx = %.3lf\t\t\t\t\t\t\t\t\t\tШаг по оси y, hy = %.3lf\n\n", hx, hy);

		if (f == '1')
			printf("\tФункция внутренних источников теплоты f(x, y, t) = A*t*EXP(-psi*t/10)*SIN(PI*x*t/4)*COS(PI*y*t/2)\n\n");
		else
			printf("\tФункция внутренних источников теплоты f(x, y, t) = 0\n\n");

		printf("\tРаспределение температуры в пластине (матрица T[%dx%d]):\n", K1, M1);
		cout << endl << T.Distribution(t, lambda_x, lambda_y, a_x, a_y, heat, density, num, f);
		cout << endl << endl << "\tКоличество итераций по времени t: " << num;

		//записываем распределение в файл
		T.Write_in(file, T.Distribution(t, lambda_x, lambda_y, a_x, a_y, heat, density, num, f));
		printf("\n\n\n\n\n\n\n");

		//обратная задача
		//для заданного времени t=3 находим значение функции T
		//вызов функции восстановления
		Recovery2(t_, lambda_x, lambda_y, a_x, a_y, lam_x, lam_y, a_xx, a_yy, heat, density, f, m, k);

		cout << "\tРешение обратной задачи - восстановление теплофизических характеристик:\n\n";

		cout << "\tКоэффициент теплопроводности по x:\n";
		cout << "Заданное значение lambda_x:\n" << lambda_x << endl;
		cout << "Восстановленное значение lam_x:\n" << lam_x << endl;

		cout << "\tКоэффициент теплопроводности по y:\n";
		cout << "Заданное значение lambda_y\n" << lambda_y << endl;
		cout << "Восстановленное значение lam_y\n" << lam_y << endl;

		cout << "\tКоэффициент температуропроводности по x:\n";
		cout << "Заданное значение a_x\n" << a_x << endl;
		cout << "Восстановленное значение a_xx\n" << a_xx << endl;

		cout << "\tКоэффициент температуропроводности по y:\n";
		cout << "Заданное значение a_y\n" << a_y << endl;
		cout << "Восстановленное значение a_yy\n" << a_yy << endl;

		cout << "\n\n\tИзмеренное распределение температуры при t=3:\n\n";
		cout << T.Distribution(t_, lambda_x, lambda_y, a_x, a_y, heat, density, num, f) << endl << endl;
		cout << "\tВосстановленное распределение температуры при t=3:\n\n";
		cout << T.Distribution(t_, lam_x, lam_y, a_xx, a_yy, heat, density, num, f);
		printf("\n\n\n\n\n\n\n");
		break;
	}	

	case '3':
	{
		system("cls");  //очистка экрана

		double conductivity_x = 130;                          //теплопроводность  [Вт/(м*град)]
		double conductivity_y = 180;                          //теплопроводность  [Вт/(м*град)]
		double heat = 0.460;                                  //удельная теплоемкость  [кДж/(кг*град)]
		double density = 7870;                                //плотность  [кг/(м^3)]

		for (double j(0), j1(0); j < M1, j1 <= lx + eps; j++, j1 += hx)
		{
			lambda_x[j] = conductivity_x + 20 * sin(4 * j1);
			a_x[j] = 10 * lambda_x[j] / (heat * density);
		}

		for (double i(0), i1(0); i < K1, i1 <= ly + eps; i++, i1 += hy)
		{
			lambda_y[i] = conductivity_y - 20 * sqrt(i1) + 10 * sin(i1);
			a_y[i] = 10 * lambda_y[i] / (heat * density);
		}

		//прямая задача
		printf("\tТест №3 - тригонометрические функции\t\t\t\t\t\t\t\t\tВведенное время t = %.1f\n", t);
		printf("\tКоэффициенты теплопроводности lambda_x = %.0lf+20sin(4x); lambda_y = %.0lf-20sqrt(y)+10sin(y)\t\tДлина пластины lx = %.0lf\n", conductivity_x, conductivity_y, lx);
		printf("\tКоэффициенты температуропроводности a_x = lambda_x/(heat*density); a_y = lambda_y/(heat*density)\tШирина пластины ly = %.0lf\n", ly);
		printf("\tВоздействие теплового потока q = %d + 25*t\t\t\t\t\t\t\t\tКоэффициент теплоотдачи от стенки к среде alpha = %d\n", q, alpha);
		printf("\tТемпература окружающей среды Tc = %d\t\t\t\t\t\t\t\t\tШаг по времени t, tau = %.3lf\n", Tc, tau);
		printf("\tШаг по оси x, hx = %.3lf\t\t\t\t\t\t\t\t\t\tШаг по оси y, hy = %.3lf\n\n", hx, hy);

		if (f == '1')
			printf("\tФункция внутренних источников теплоты f(x, y, t) = A*t*EXP(-psi*t/10)*SIN(PI*x*t/4)*COS(PI*y*t/2)\n\n");
		else
			printf("\tФункция внутренних источников теплоты f(x, y, t) = 0\n\n");

		printf("\tРаспределение температуры в пластине (матрица T[%dx%d]):\n", K1, M1);
		cout << endl << T.Distribution(t, lambda_x, lambda_y, a_x, a_y, heat, density, num, f);
		cout << endl << endl << "\tКоличество итераций по времени t: " << num;

		//записываем распределение в файл
		T.Write_in(file, T.Distribution(t, lambda_x, lambda_y, a_x, a_y, heat, density, num, f));
		printf("\n\n\n\n\n\n\n");


		//обратная задача
		//для заданного времени t=3 находим значение функции T
		//вызов функции восстановления
		Recovery2(t_, lambda_x, lambda_y, a_x, a_y, lam_x, lam_y, a_xx, a_yy, heat, density, f, m, k);

		cout << "\tРешение обратной задачи - восстановление теплофизических характеристик:\n\n";

		cout << "\tКоэффициент теплопроводности по x:\n";
		cout << "Заданное значение lambda_x:\n" << lambda_x << endl;
		cout << "Восстановленное значение lam_x:\n" << lam_x << endl;

		cout << "\tКоэффициент теплопроводности по y:\n";
		cout << "Заданное значение lambda_y\n" << lambda_y << endl;
		cout << "Восстановленное значение lam_y\n" << lam_y << endl;

		cout << "\tКоэффициент температуропроводности по x:\n";
		cout << "Заданное значение a_x\n" << a_x << endl;
		cout << "Восстановленное значение a_xx\n" << a_xx << endl;

		cout << "\tКоэффициент температуропроводности по y:\n";
		cout << "Заданное значение a_y\n" << a_y << endl;
		cout << "Восстановленное значение a_yy\n" << a_yy << endl;

		cout << "\n\n\tИзмеренное распределение температуры при t=3:\n\n";
		cout << T.Distribution(t_, lambda_x, lambda_y, a_x, a_y, heat, density, num, f) << endl << endl;
		cout << "\tВосстановленное распределение температуры при t=3:\n\n";
		cout << T.Distribution(t_, lam_x, lam_y, a_xx, a_yy, heat, density, num, f);
		printf("\n\n\n\n\n\n\n");
		break;
	}

	default:
	{
		cout << "\nНеверный пункт меню!\n\n\n\n\n\n\n";
	}
	}

	file.close();  //закрытие файла на запись

	system("pause");
	return 0;
}