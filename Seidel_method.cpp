#include"matrix_lib.h"

std::vector<double> slay::Seidel_method::get_solution(matr matrix, std::vector<double> value_vector) {
	std::vector<double> solution = value_vector;
	slay::matrix mat;

	while (mat.GetVectorNorm(mat.GetResidualVector(matrix, value_vector, solution)) > std::pow(0.1, 5)) {
		this->operation_count++;
		std::vector<double> tempSol(matrix.size(), 0);
		for (int i = 0; i < matrix.size(); i++) {
			tempSol[i] = value_vector[i];
			for (int j = 0; j < matrix.size(); j++) {
				if (i != j) {
					tempSol[i] -= (matrix[i][j] * solution[j]);
				}
			}
			solution[i] = tempSol[i] / matrix[i][i];
		}
	}
	return solution;
}
