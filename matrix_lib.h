#pragma once
#include<iostream>
#include<vector>
#include<algorithm>
#include<cmath>

#define matr std::vector<std::vector<double>>
#define M_PI 3.141592


namespace slay {

	class matrix
	{
	public:
		static std::vector<double> normalization_vector(std::vector<double> vect);
		static double scalar_prodact(std::vector<double> first, std::vector<double> second);
		static void SwapColumnes(int IndexFirst, int IndexSecond, matr* matrix);
		static void SwapLines(int indexFirst, int indexSecond, matr* matrix);
		static std::pair<int, int> FindIndexMaxInMainMinor(int minorNumber, matr matrix);
		static std::vector<double> GetResidualVector(matr matrix, std::vector<double> value_vector, std::vector<double> solution);
		static std::vector<std::vector<double>> getResidualMatrix(matr reverseMatrix);
		//make this only for martix
		static std::vector<std::vector<double>> multiplicationMatrix(matr first, matr second);
		//matrix x vector
		static std::vector<double> multiplicationMatrix(matr first, std::vector<double> second);
		static std::vector<std::vector<double>> multiplicationMatrix(std::vector<std::vector<double>> matrix, double value);
		static double GetVectorNorm(std::vector<double> solution);
		static double GetMatrixNorm(std::vector<std::vector<double>> residualMatrix);
		static bool equels(matr first, matr second);
		static std::vector<double> getSolutionFromUpperTriangular(matr matrix, std::vector<double> valueVector);
		static matr GetTranspose(matr matrix);
	};

	class Gaus {
	public:
		Gaus();
		std::vector<double> getSolution(matr matrix, std::vector<double> valueVector);
		std::vector<std::vector<double>> GetReverseMatrix(matr matrix);
		double GetDet(matr matrix);
		matr GetLastSolveTriangular();

	private:
		double determinant;
		matr upperTriangular;
		matr lastSolve;
	};


	class Jacobi_method {
	public:
		std::vector<double> GetSolution(matr matrix, std::vector<double> value_vector);
		double get_count() {
			return this->operation_count;
		}
		int get_oreor_am_iter(matr matrix,std::vector<double> value_vector, double esp = 0) {
			if (esp == 0) {
				esp = std::pow(0.1, 5);
			}
			double k = std::log((esp * (1 - slay::matrix::GetMatrixNorm(matrix)) / slay::matrix::GetVectorNorm(value_vector)))
				/ std::log(slay::matrix::GetMatrixNorm(matrix));

			k = std::abs(k);
			return k;
		}
	private:
		int operation_count = 0;
	};


	class Seidel_method {
	public:
		std::vector<double> get_solution(matr matrix, std::vector<double> value_vector);
		double get_count() {
			return this->operation_count;
		}
	private:
		int operation_count = 0;
	};


	class Min_residual_method {
	public:
		std::vector<double> get_solution(matr matrix, std::vector<double> valueVector);
		double get_count() {
			return this->operation_count;
		}
	private:
		int operation_count = 0;
	};

	//для симметрических положительно определенных матриц
	class Relaxation_method {
	public:
		std::vector<double> get_solution(matr matrix, std::vector<double> value_vector);
		double get_count() {
			return this->operation_count;
		}
	private:
		int operation_count = 0;
	};
}


namespace eigenvalues {

	class common {
	public:
		static std::vector<double> get_residual_vector(std::vector<std::vector<double>> matrix, std::vector<double> eigen_vector, double lambda);
		static double get_resid_from_polinom(double lambda, std::vector<double> coef) {
			std::vector<double> new_coef;
			new_coef.push_back(1);
			for (auto it : coef) {
				new_coef.push_back(it);
			}
			double temp = 0;
			for (int i = new_coef.size() -1, d = 0; i >= 0; i--,d++) {
				temp += std::pow(lambda, i +1) * new_coef[d];
			}
			return temp;
		}
	};



	class Krilov_method {
	public:
		std::vector<double> get_polinom_coef(matr matrix, std::vector<double> random_vector = {});
		std::vector<double> get_eigen_vector(double lambda);

	private:
		std::vector<std::vector<double>> solve_vectors;
		std::vector<double> cof;
	};


	class Danilevski_method {
	public:
		std::vector<double> get_polinom_coef(std::vector<std::vector<double>> matrix);
		std::vector<double> get_eigen_vector(double lambda);
		

	private:
		std::vector<matr> transform;
		std::vector<std::vector<double>> last_matrix;
	};

	class Jacobi_method {
	public:

		std::vector<double> get_eigen_values(std::vector<std::vector<double>> matrix);
		std::pair<double, std::pair<int, int>> find_appropriate(matr matrix);
		matr get_all_eigen_vecor();
		double get_count() {
			return this->operation_count;
		}
	private:
		int operation_count = 0;
		matr e_vectors;
	};

	class Power_method {
	public:

		void run(matr matrix);
		double get_lambda();
		std::vector<double> get_eigen_vector();
		double get_count() {
			return this->operation_count;
		}

	private:
		int operation_count = 0;
		double lambda;
		std::vector<double> e_vector;
	};
}