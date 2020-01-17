#include"matrix_lib.h"

std::vector<double> slay::Relaxation_method::get_solution(matr matrix, std::vector<double> value_vector) {

	std::vector<double> solution = value_vector;
	double w = 1.4;
	double esp = std::powf(0.1, 5);
	auto resid = slay::matrix::GetResidualVector(matrix, value_vector, solution);
	double norm = slay::matrix::GetVectorNorm(resid);

	while (norm > esp) {
		this->operation_count++;
		std::vector<double> temp_solution(matrix.size(), 0);
		for (int i = 0; i < matrix.size(); i++) {

			double first_sum = 0;
			for (int j = 0; j < i; j++) {
				first_sum += ((matrix[i][j] / matrix[i][i]) * temp_solution[j]);
			}

			double second_sum = 0;
			for (int j = i + 1; j < matrix.size(); j++) {
				second_sum += ((matrix[i][j] / matrix[i][i]) * solution[j]);
			}

			double third_sum = 0;
			third_sum = value_vector[i] / matrix[i][i];

			temp_solution[i] = (1 - w) * solution[i] - w * first_sum - w * second_sum + w * third_sum;
		}
		solution.assign(temp_solution.begin(), temp_solution.end());
		resid = slay::matrix::GetResidualVector(matrix, value_vector, solution);
		norm = slay::matrix::GetVectorNorm(resid);
	}
	return solution;
}