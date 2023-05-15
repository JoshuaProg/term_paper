#include <iostream>
#include <vector>
#include <math.h>
using namespace std;

const float e = 0.0001;


float holschmidt_method(float, float);
float geron(float);
vector<float> square_root_method(vector<vector<float>>, vector<float>);
bool check_matrix(vector<vector<float>>);
void print_matrix(vector<vector<float>>);
void print_matrix(vector<float>);
float determinant(vector<vector<float>>);

int main() {
	setlocale(LC_ALL, "ru");

	bool exit = false;
	while (!exit) {
		int number = 0;
		cout << "Выберите алгоритм решения СЛАУ методом квадратного корня: \n1. Метод Герона \n2 Алгоритм Гольдшмидта \n3. Метод квадратного корня \n4. Выход \nВвод: ";
		cin >> number;
		if (number == 1) { // алгоритм Гольдшмидта
			float N, D;
			cout << "Введите делимое: ";
			cin >> N;
			cout << "Введите делитель: ";
			cin >> D;
			cout << "Результат алгоритма: \n";
			cout << holschmidt_method(N, D);
			cout << "\n\n";
		}
		else if (number == 2) { // Метод Герона
			float b;
			cout << "Введите число: ";
			cin >> b;
			cout << "Результат алгоритма: \n";
			if (geron(b) == -1)
				cout << "ЧИСЛО ВВЕДЕНО НЕВЕРНО";
			else cout << geron(b);
			cout << "\n\n";
		}

		else if (number == 3) { // Метод квадратного корня
			int rows, columns;
			cout << "Введите число строк матрицы A: ";
			cin >> rows;
			cout << "Введите число столбцов матрицы A: ";
			cin >> columns;
			cout << "Введите матрицу: \n";
			vector<vector<float>> matrix_A(rows, vector<float>(columns));
			for (int i = 0; i < rows; i++)
				for (int j = 0; j < columns; j++)
					cin >> matrix_A[i][j];
			cout << "Введите матрицу B: \n";
			vector<float> matrix_B(rows);
			for (int i = 0; i < rows; i++)
				cin >> matrix_B[i];
			cout << "Результат алгоритма: \n";
			vector<float> result = square_root_method(matrix_A, matrix_B);
			if (result.size() != 0)
				print_matrix(result);
			else cout << "МАТРИЦА ВВЕДЕНА НЕВЕРНО";
			cout << "\n\n";
		}
		else if (number == 0)
			exit = true;
		else
			cout << "ТАКОГО НОМЕРА НЕТ\n\n";
	}
}

float holschmidt_method(float N, float D) {
	float result;
	N /= D * 1.5;
	D /= D * 1.5; 
	float old_N;

	do {
		float F = 2 - D;
		old_N = N;
		N *= F;
		D *= F;

	} while (abs(old_N - N) >= e);

	return N;
}

bool check_matrix(vector<vector<float>> matrix) {
	if (matrix.size() != matrix.back().size())
		return false;

	for (int i = 0; i < matrix.size(); i++)
		for (int j = i + 1; j < matrix.size(); j++)
			if (matrix[i][j] != matrix[j][i])
				return false;

	if (determinant(matrix) == 0)
		return false;

	return true;
}

float determinant(vector<vector<float>> matrix) {
	float temp = 0;

	for (int k = 0; k < matrix.size() - 1; k++)
		for (int i = k + 1; i < matrix.size(); i++) {
			temp = -matrix[i][k] / matrix[k][k];
			for (int j = 0; j < matrix.size(); j++)
				matrix[i][j] += matrix[k][j] * temp;
		}

	float d = 1;
	for (int i = 0; i < matrix.size(); i++)
		d *= matrix[i][i];

	return d;
}


void print_matrix(vector<vector<float>> matrix) {
	for (int i = 0; i < matrix.size(); i++) {
		for (int k = 0; k < matrix[i].size(); k++)
			cout << matrix[i][k] << " ";
		cout << endl;
	}
}

void print_matrix(vector<float> matrix) {
	for (int i = 0; i < matrix.size(); i++)
		cout << matrix[i] << " ";
}
vector<float> square_root_method(vector<vector<float>> A, vector<float> B) {
	if (check_matrix(A)) {
		vector<vector<float>> new_T(A.size(), vector<float>(A.size()));
		vector<float> y(A.size()), x(A.size());

		for (int i = 0; i < new_T.size(); i++)
			for (int j = 0; j < new_T[i].size(); j++)
				new_T[i][j] = 0;
		for (int i = 0; i < A.size(); i++) {
			float temp = 0;
			for (int k = 0; k < i; k++)
				temp += new_T[k][i] * new_T[k][i];
			new_T[i][i] = sqrt(A[i][i] - temp);
			for (int j = i; j < A.size(); j++) {
				temp = 0;
				for (int k = 0; k < i; k++)
					temp += new_T[k][i] * new_T[k][j];
				new_T[i][j] = (A[i][j] - temp) / new_T[i][i];
			}
		}

		for (int i = 0; i < A.size(); i++)
		{
			float temp = 0;
			for (int k = 0; k < i; k++)
				temp = temp + new_T[k][i] * y[k];
			y[i] = (B[i] - temp) / new_T[i][i];
		}
		for (int i = A.size() - 1; i >= 0; i--)
		{
			float temp = 0;
			for (int k = i + 1; k < A.size(); k++)
				temp = temp + new_T[i][k] * x[k];
			x[i] = (y[i] - temp) / new_T[i][i];
		}
		return x;
	}
	else return vector<float>(0);
}
float geron(float number) {
	if (number < 0) return -1;
	float y = 1, z = 0;
	while (abs(z - y) >= e) {
		z = y;
		y = (y + number / y) / 2;
	}
	return y;
}
