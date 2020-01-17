#include"matrix_lib.h"

void eigenvalues::Power_method::run(matr matrix) {
	std::vector<double> y_prev;
	std::vector<double> y_current(matrix.size(), 1);
	double lambda_prev;
	double lambda_current = 0;
	double esp = std::pow(10, -6);
	do
	{
		this->operation_count++;
		lambda_prev = lambda_current;
		y_prev = y_current;
		y_current = slay::matrix::multiplicationMatrix(matrix, y_current);
		double temp = 0;
		for (int i = 0; i < matrix.size(); i++) {
			temp += y_current[i] / y_prev[i];
		}
		lambda_current = temp / matrix.size();
		y_current = slay::matrix::normalization_vector(y_current);//нормализуем вектор
	} while (std::abs(lambda_current - lambda_prev) > esp);
	this->lambda = lambda_current;
	this->e_vector = y_current;
}

double eigenvalues::Power_method::get_lambda() {
	return this->lambda;
}

std::vector<double> eigenvalues::Power_method::get_eigen_vector() {
	return this->e_vector;
}