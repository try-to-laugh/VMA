#include"matrix_lib.h"

std::vector<double> slay::matrix::normalization_vector(std::vector<double> vect) {
	double coff = std::sqrt(matrix::scalar_prodact(vect, vect));
	for (int i = 0; i < vect.size(); i++) {
		vect[i] /= coff;
	}
	return vect;
}

double slay::matrix::scalar_prodact(std::vector<double> first, std::vector<double> second) {
	if (first.size() != second.size()) {
		std::cout << "Exeption! slay::matrix::multiplicationMatrix - \"unacceptable matrix size\"" << std::endl;
		exit(3);
	}

	double result = 0;
	for (int i = 0; i < first.size(); i++) {
		result += first[i] * second[i];
	}
	return result;
}

void slay::matrix::SwapColumnes(int IndexFirst, int IndexSecond, matr* matrix) {
	for (int i = 0; i < (*matrix).size(); i++) {
		std::swap((*matrix)[i][IndexFirst], (*matrix)[i][IndexSecond]);
	}
}

void slay::matrix::SwapLines(int indexFirst, int indexSecond, matr* matrix) {
	(*matrix)[indexFirst].swap((*matrix)[indexSecond]);
}

std::pair<int, int> slay::matrix::FindIndexMaxInMainMinor(int minorNumber, matr matrix)
{
	double max = std::abs(matrix[minorNumber][minorNumber]);
	std::pair<int, int> MaxIndex = std::make_pair(minorNumber, minorNumber);
	for (int i = minorNumber; i < matrix.size(); i++) {
		for (int j = 0; j < matrix.size(); j++) {
			if (std::abs(matrix[i][j]) > max) {
				max = std::abs(matrix[i][j]);
				MaxIndex = std::make_pair(i, j);
			}
		}
	}
	return MaxIndex;
}

std::vector<double> slay::matrix::GetResidualVector(matr matrix, std::vector<double> value_vector, std::vector<double> solution)
{
	std::vector<double> error;
	for (size_t i = 0; i < matrix.size(); i++) {
		double temp = 0;
		for (int g = 0; g < matrix[i].size(); g++) {
			temp += (matrix[i][g] * solution[g]);
		}
		error.push_back(temp - value_vector[i]);
	}
	return error;
}

std::vector<std::vector<double>> slay::matrix::getResidualMatrix(matr reverseMatrix)
{
	for (int i = 0; i < reverseMatrix.size(); i++) {
		reverseMatrix[i][i] -= 1;
	}
	return reverseMatrix;
}


//make this only for martix
std::vector<std::vector<double>> slay::matrix::multiplicationMatrix(matr first, matr second) {

	if (first.empty() || second.empty()) {
		std::cout << "Exeption! slay::matrix::multiplicationMatrix - \"empty matrix\"" << std::endl;
	}

	bool correct = true;
	if (first.size() != second.size()) {
		correct = false;
	}
	int size = first.size();
	for (int i = 0; i < size; i++) {
		if ((first[i].size() != size) || (second[i].size() != size)) {
			correct = false;
			break;
		}
	}

	if (!correct) {
		std::cout << "Exeption! slay::matrix::multiplicationMatrix - \"unacceptable matrix size\"" << std::endl;
		exit(3);
	}


	std::vector<double> initializer;
	matr result(size, initializer);

	for (int i = 0; i < size; i++) {
		for (int d = 0; d < size; d++) {
			double term = 0;
			for (int k = 0; k < size; k++) {
				term += (first[d][k] * second[k][i]);
			}
			result[d].push_back(term);
		}
	}
	return result;
}


std::vector<double> slay::matrix::multiplicationMatrix(matr first, std::vector<double> second) {
	std::vector<double> result;
	if (first.size() != second.size()) {
		std::cout << "Exeption! slay::matrix::multiplicationMatrix - \"unacceptable matrix size\"" << std::endl;
		exit(3);
	}

	for (int d = 0; d < second.size(); d++) {
		double term = 0;
		for (int k = 0; k < second.size(); k++) {
			term += (first[d][k] * second[k]);
		}
		result.push_back(term);
	}

	return result;
}


std::vector<std::vector<double>>slay::matrix::multiplicationMatrix(std::vector<std::vector<double>> matrix, double value)
{
	for (int i = 0; i < matrix.size(); i++) {
		for (int j = 0; j < matrix.size(); j++) {
			matrix[i][j] *= value;
		}
	}
	return matrix;
}

double slay::matrix::GetVectorNorm(std::vector<double> solution)
{
	double max = std::abs(solution[0]);
	for (int i = 1; i < solution.size(); i++) {
		if (std::abs(solution[i]) > max) {
			max = std::abs(solution[i]);
		}
	}
	return max;
}

double slay::matrix::GetMatrixNorm(std::vector<std::vector<double>> residualMatrix)
{
	double max = 0;
	for (int i = 0; i < residualMatrix.size(); i++) {
		max += std::abs(residualMatrix[0][i]);
	}

	for (int i = 1; i < residualMatrix.size(); i++) {
		double temp = 0;
		for (int j = 0; j < residualMatrix.size(); j++) {
			temp += std::abs(residualMatrix[i][j]);
		}
		if (temp > max) {
			max = temp;
		}
	}
	return max;
}

bool slay::matrix::equels(matr first, matr second) {
	if ((first.size() != second.size()) || (first[0].size() != second.size())) {
		return false;
	}
	for (int i = 0; i < first.size(); i++) {
		for (int j = 0; j < second.size(); j++) {
			if (first[i][j] != second[i][j]) {
				return false;
			}
		}
	}
	return true;
}

std::vector<double> slay::matrix::getSolutionFromUpperTriangular(matr matrix, std::vector<double> valueVector) {
	std::vector<double> solution(matrix.size(), 0);
	for (int i = matrix.size() - 1; i >= 0; i--) {
		double temp = valueVector[i];
		for (int j = matrix.size() - 1; j > i; j--) {
			temp -= (solution[j] * matrix[i][j]);
		}
		solution[i] = temp;
	}
	return solution;
}

matr slay::matrix::GetTranspose(matr matrix) {
	std::vector<double> initialize;
	matr transpose(matrix.size(), initialize);
	for (int i = 0; i < matrix.size(); i++) {
		for (int j = 0; j < matrix.size(); j++) {
			transpose[j].push_back(matrix[i][j]);
		}
	}
	return transpose;
}

