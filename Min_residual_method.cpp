#include"matrix_lib.h"

std::vector<double> slay::Min_residual_method::get_solution(matr matrix, std::vector<double> valueVector) {

	valueVector = matrix::multiplicationMatrix(matrix::GetTranspose(matrix), valueVector);
	matrix = matrix::multiplicationMatrix(matrix::GetTranspose(matrix), matrix);

	std::vector<double> solution = valueVector;
	slay::matrix mat;
	std::vector<double> resid = mat.GetResidualVector(matrix, valueVector, solution);
	double esp = powf(0.1, 5);
	while (mat.GetVectorNorm(resid) > esp) {
		this->operation_count++;
		std::vector<double> mr = mat.multiplicationMatrix(matrix, resid);
		double t = (mat.scalar_prodact(mr, resid)) / (mat.scalar_prodact(mr, mr));

		for (int i = 0; i < matrix.size(); i++) {
			solution[i] -= t * resid[i];
		}

		resid = mat.GetResidualVector(matrix, valueVector, solution);
	}
	return solution;
}