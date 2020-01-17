#pragma once
#include<iostream>
#include<vector>
#include<fstream>
#include<iomanip>

namespace custom_in_out {
	template<class T>
	void ConsolePrintVector(std::vector<T> vect) {
		for (int i = 0; i < vect.size(); i++) {
			std::cout << std::setw(15) << std::right << vect[i] << " ";
		}
		std::cout << std::endl;
	}

	template<class T>
	void ConsolePrintMatrix(std::vector<std::vector<T>> mat) {
		for (int i = 0; i < mat.size(); i++) {
			for (int j = 0; j < mat[i].size(); j++) {
				std::cout << std::setw(15) << std::right << mat[i][j] << " ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}

	template<class T>
	std::vector<std::vector<T>> inputMatrix(std::istream& fin) {

		std::vector<std::vector<T>> matrix;
		int size;
		fin >> size;
		for (int i = 0; i < size; i++) {
			std::vector<T> temp;
			for (int j = 0; j < size; j++) {
				T var;
				fin >> var;
				temp.push_back(var);
			}
			matrix.push_back(temp);
		}

		return matrix;
	}

	template<class T>
	std::vector<T> inputVector(std::istream& fin) {
		std::vector<T> vect;
		int size;
		fin >> size;
		for (int i = 0; i < size; i++) {
			T temp;
			fin >> temp;
			vect.push_back(temp);
		}
		return vect;
	}
}

	


