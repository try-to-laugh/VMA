#include"matrix_lib.h"


std::vector<double> eigenvalues::Krilov_method::get_polinom_coef(matr matrix, std::vector<double> random_vector)
{
	if (random_vector.empty()) {
		random_vector.assign(matrix.size(), 0);
		random_vector[0] = 1;
	}

	std::vector<std::vector<double>> cf;
	cf.push_back(random_vector);
	for (int i = 0; i < matrix.size() - 1; i++) {
		random_vector = slay::matrix::multiplicationMatrix(matrix, random_vector);
		cf.push_back(random_vector);
	}

	random_vector = slay::matrix::multiplicationMatrix(matrix, random_vector);
	std::reverse(cf.begin(), cf.end());
	cf = slay::matrix::GetTranspose(cf);
	this->solve_vectors = cf;


	slay::Gaus temp;
	auto coef = temp.getSolution(cf, random_vector);
	for (int i = 0; i < coef.size(); i++) {
		coef[i] *= -1;
	}
	this->cof = coef;
	return coef;
}


std::vector<double> eigenvalues::Krilov_method::get_eigen_vector(double lambda) {
	std::vector<double> for_beta;
	double beta = 1;
	for_beta.push_back(beta);
	for (int i = 0; i < this->cof.size(); i++) {
		beta = beta * lambda + cof[i];
		for_beta.push_back(beta);
	}

	std::vector<double> eigen_vector;
	for (int i = 0; i < this->solve_vectors.size(); i++) {
		double temp = 0;
		for (int j = 0; j < this->solve_vectors.size(); j++) {
			temp += this->solve_vectors[i][j] * for_beta[j];
		}
		eigen_vector.push_back(temp);
	}
	return eigen_vector;
}
