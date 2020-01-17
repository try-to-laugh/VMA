#include"matrix_lib.h"

slay::Gaus::Gaus() {
	this->determinant = 1;
}

std::vector<double> slay::Gaus::getSolution(matr matrix, std::vector<double> valueVector) {

	if ((matrix.size() != matrix[0].size()) || matrix.size() != valueVector.size()) {
		std::cout << "slay::Gaus::getSolution - \"invalid matrix/vector size\"" << std::endl;
		exit(3);
	}

	std::vector<std::pair<int, int>> permutation;
	this->determinant = 1;
	this->lastSolve = matrix;
	slay::matrix def;

	for (int i = 0; i < matrix.size(); i++) {

		std::pair<int, int> IndexMax = def.FindIndexMaxInMainMinor(i, matrix);
		def.SwapLines(i, IndexMax.first, &matrix);
		def.SwapColumnes(i, IndexMax.second, &matrix);
		std::swap(valueVector[i], valueVector[IndexMax.first]);
		permutation.push_back(std::make_pair(i, IndexMax.second));

		double LeadElement = matrix[i][i];
		if (LeadElement == 0) {
			std::cout << "slay::Gaus::getSolution - \"degenerate matrix\" " << std::endl;
			exit(3);
		}

		this->determinant *= LeadElement;

		for (int j = i; j < matrix.size(); j++) {
			matrix[i][j] /= LeadElement;
		}
		valueVector[i] /= LeadElement;

		for (int f = (i + 1); f < matrix.size(); f++) {
			double factor = matrix[f][i];
			for (int k = i; k < matrix.size(); k++) {
				matrix[f][k] -= (matrix[i][k] * factor);
			}
			valueVector[f] -= (valueVector[i] * factor);
		}
	}


	this->upperTriangular = matrix;
	std::vector<double> solution(matrix.size(), 0);
	for (int i = matrix.size() - 1; i >= 0; i--) {
		double temp = valueVector[i];
		for (int j = matrix.size() - 1; j > i; j--) {
			temp -= (solution[j] * matrix[i][j]);
		}
		solution[i] = temp;
	}


	std::reverse(permutation.begin(), permutation.end());
	for (int i = 0; i < permutation.size(); i++) {
		std::swap(solution[permutation[i].first], solution[permutation[i].second]);
	}
	return solution;
	std::reverse(permutation.begin(), permutation.end());
	for (auto psw : permutation) {
		std::swap(solution[psw.first], solution[psw.second]);
	}
}

std::vector<std::vector<double>> slay::Gaus::GetReverseMatrix(matr matrix)
{
	matr reverseMatrix;
	for (int i = 0; i < matrix.size(); i++) {
		std::vector<double> temp(matrix.size(), 0);
		reverseMatrix.push_back(temp);
	}

	for (int columnCount = 0; columnCount < matrix.size(); columnCount++) {

		std::vector<double> tempValueVector(matrix.size(), 0);
		tempValueVector[columnCount] = 1;
		std::vector<double> tempSolution = this->getSolution(matrix, tempValueVector);
		for (int i = 0; i < matrix.size(); i++) {
			reverseMatrix[i][columnCount] = tempSolution[i];
		}
	}
	return reverseMatrix;
}

double slay::Gaus::GetDet(matr matrix) {
	if (slay::matrix::equels(matrix, lastSolve)) {
		return this->determinant;
	}

	Gaus temp;
	std::vector<double> valV(matrix.size(), 0);
	valV[0] = 1;
	temp.getSolution(matrix, valV);


	return temp.GetDet(matrix);
}

matr slay::Gaus::GetLastSolveTriangular() {
	return lastSolve;
}