#include"matrix_lib.h"

std::vector<double>eigenvalues::common::get_residual_vector(std::vector<std::vector<double>> matrix, std::vector<double> eigen_vector, double lambda) {
	for (int i = 0; i < matrix.size(); i++) {
		matrix[i][i] -= lambda;
	}
	return slay::matrix::multiplicationMatrix(matrix, eigen_vector);
}