#include"matrix_lib.h"

std::vector<double> eigenvalues::Danilevski_method::get_polinom_coef(std::vector<std::vector<double>> matrix) {

	std::vector<double> coef;
	this->last_matrix = matrix;

	for (int i = matrix.size() - 2; i >= 0; i--) {
		
		matr copy_matrix = matrix;
		std::vector<double> m_obr = matrix[i + 1];
		double devider = matrix[i + 1][i];

		matr to_eigen_vectors;
		{
			std::vector<double> init(matrix.size(), 0);
			for (int k = 0; k < matrix.size(); k++) {
				init[k] = 1;
				to_eigen_vectors.push_back(init);
				init[k] = 0;
			}

			for (int k = 0; k < matrix.size(); k++) {
				if (k == i) {
					to_eigen_vectors[i][k] = 1 / devider;
				}
				else {
					to_eigen_vectors[i][k] = -1 * matrix[i + 1][k] / devider;
				}
			}
		}


		this->transform.push_back(to_eigen_vectors);
		
		for (int k = 0; k <= (i + 1); k++) {
			matrix[k][i] /= devider;
		}
		
		for (int k = 0; k < matrix.size(); k++) {
			if (k != i) {
				double div = matrix[i + 1][k];
				for (int j = 0; j <= (i + 1); j++) {
					matrix[j][k] -= (div * matrix[j][i]);
				}
			}
		}

		std::vector<double> initialaizer(matrix.size(), 0);
		std::vector<std::vector<double>> obr;
		for (int k = 0; k < matrix.size(); k++) {
			if (k == i) {
				obr.push_back(m_obr);
				continue;
			}
			initialaizer[k] = 1;
			obr.push_back(initialaizer);
			initialaizer[k] = 0;
		}
		matrix = slay::matrix::multiplicationMatrix(obr, matrix);
	}
	for (int i = 0; i < matrix.size(); i++) {
		coef.push_back(-matrix[0][i]);
	}
	std::reverse(transform.begin(), transform.end());
	return coef;
}


std::vector<double> eigenvalues::Danilevski_method::get_eigen_vector(double lambda) {
	std::vector<double> eigen_vector;
	double item = 1;
	for (int i = 0; i < last_matrix.size(); i++) {
		eigen_vector.push_back(item);
		item *= lambda;
	}
	std::reverse(eigen_vector.begin(), eigen_vector.end());
	std::vector<double> copy = eigen_vector;


	for (auto it : transform) {
		eigen_vector = slay::matrix::multiplicationMatrix(it, eigen_vector);
	}

	return eigen_vector;
}