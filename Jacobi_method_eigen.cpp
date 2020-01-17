#include"matrix_lib.h"


std::vector<double> eigenvalues::Jacobi_method::get_eigen_values(std::vector<std::vector<double>> matrix) {

	matr eigen_vectors;
	{
		std::vector<double> init(matrix.size(), 0);
		for (int i = 0; i < matrix.size(); i++) {
			init[i] = 1;
			eigen_vectors.push_back(init);
			init[i] = 0;
		}
	}

	std::vector<double> lambdas;
	double two_tan;
	double sin;
	double cos;
	double esp = std::pow(10, -6);
	auto max_el = this->find_appropriate(matrix);
	double max = max_el.first;
	std::pair<int, int> index = max_el.second;


	while (max > esp) {
		this->operation_count++;
		if (matrix[index.first][index.first] == matrix[index.second][index.second]) {
			two_tan = std::tan(M_PI / 2);
			sin = std::sin(M_PI / 4);
			cos = std::cos(M_PI / 4);
		}
		else {
			two_tan = 2 * matrix[index.first][index.second] /
				(matrix[index.first][index.first] - matrix[index.second][index.second]);


			double temp = 1 / std::sqrt(1 + std::pow(two_tan, 2));
			int sign = two_tan / std::abs(two_tan);
			cos = std::sqrt(0.5 * (1 + temp));
			sin = sign * std::sqrt(0.5 * (1 - temp));
		}

		std::vector<double> init(matrix.size(), 0);
		matr obr(matrix.size(), init);
		matr norm(matrix.size(), init);
		norm[index.first][index.first] = norm[index.second][index.second] = cos;
		norm[index.first][index.second] = -1 * sin;
		norm[index.second][index.first] = sin;
		for (int i = 0; i < matrix.size(); i++) {
			if (norm[i][i] == 0) {
				norm[i][i] = 1;
			}
		}

		obr = slay::matrix::GetTranspose(norm);

		matrix = slay::matrix::multiplicationMatrix(obr, matrix);
		matrix = slay::matrix::multiplicationMatrix(matrix, norm);

		eigen_vectors = slay::matrix::multiplicationMatrix(eigen_vectors, norm);

		max_el = this->find_appropriate(matrix);
		max = max_el.first;
		index = max_el.second;
	}

	this->e_vectors = eigen_vectors;
	for (int i = 0; i < matrix.size(); i++) {
		lambdas.push_back(matrix[i][i]);
	}

	return lambdas;
}

std::pair<double, std::pair<int, int>> eigenvalues::Jacobi_method::find_appropriate(matr matrix) {

	if (matrix.empty() || matrix.size() == 1) {
		std::cout << "eigenvalues::Jacobi_method - \"invalid size\"" << std::endl;
		system("pause");
	}

	std::pair<int, int> index = { 0,1 };
	double max = std::abs(matrix[0][1]);
	for (int i = 0; i < matrix.size(); i++) {
		for (int j = 0; j < matrix.size(); j++) {
			if ((i != j) && (std::abs(matrix[i][j]) > max)) {
				max = std::abs(matrix[i][j]);
				index.first = i;
				index.second = j;
			}
		}
	}
	return std::make_pair(max, index);
}


matr eigenvalues::Jacobi_method::get_all_eigen_vecor() {
	return this->e_vectors;
}













//попытался решить через формулки, ничего не получилось
		/*std::vector<double> k_str;
		std::vector<double> l_str;
		std::vector<double> k_column;
		std::vector<double> l_column;


		std::cout << "def: " << std::endl;
		custom_in_out::ConsolePrintMatrix<double>(matrix);

		for (int k = 0; k < matrix.size(); k++) {
			k_str.push_back(matrix[index.first][k] * cos + matrix[index.second][k] * sin);
			l_str.push_back(matrix[index.first][k] * -1 * sin + matrix[index.second][k] * cos);
		}



		matrix[index.first] = k_str;
		matrix[index.second] = l_str;

		std::cout << "swap lines: " << std::endl;
		custom_in_out::ConsolePrintMatrix<double>(matrix);

		for (int k = 0; k < matrix.size(); k++) {
			k_column.push_back(matrix[k][index.first] * cos + matrix[k][index.second] * sin);
			l_column.push_back(matrix[k][index.first] * -1 * sin + matrix[k][index.second] * cos);
		}


		for (int k = 0; k < matrix.size(); k++) {
			matrix[k][index.first] = k_column[k];
			matrix[k][index.second] = l_column[k];
		}

		std::cout << "swap col: " << std::endl;
		custom_in_out::ConsolePrintMatrix<double>(matrix);*/