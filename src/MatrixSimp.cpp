#include "MatrixSimp.h"




double determinant(const Matrix& m) {

	int rows = m.getRows();
	int cols = m.getCols();
	double retVal = 0.f;
	double coef;
	double f;

	if (rows == 2) {
		retVal = m(0, 0)*m(1, 1) - m(0, 1)*m(1, 0);
	}
	else {
		for (int c = 0;c < cols;c++) {
			if ((c + 2) % 2 == 0) { coef = 1.f; }
			else { coef = -1.f; }
			if (m(0, c) != 0){
				f = coef * m(0, c) * determinant(m.copyAndCut(0, c));
			}
			else {
				f = 0.f;
			}
			retVal += f;
		}
	}



	return retVal;

}


const Matrix adjugate(const Matrix& m) {

	int rows = m.getRows();
	int cols = m.getCols();

	Matrix adj_m(rows, cols);

	for (int r = 0;r < rows;r++) {
		for (int c = 0;c < cols;c++) {

			adj_m(c, r) = pow(-1.f, (r + c)) * determinant(m.copyAndCut(r, c));

		}
	}

	return adj_m;
}



const Matrix inverse(const Matrix& m) {
	double det = determinant(m);
	if (det == 0) {
		std::cout << m;
		assert("determinant equal to" && 0);
	}
	return adjugate(m)*(1 / det);
}


Matrix inverse(Matrix& m) {
	double det = determinant(m);
	if (det == 0) {
		std::cout << m;
		assert("determinant equal to" && 0);
	}
	return adjugate(m)*(1 / det);
}






std::ostream& operator<<(std::ostream& out, const Matrix& m)
{
	for (int r = 0;r < m.getRows();r++) {
		for (int c = 0;c < m.getCols();c++) {
			if (c != 0) {
				out << " , ";
			}
			out << m(r, c);
		}
		out << "\n";

	}
	// write obj to stream
	return out;
}


Matrix operator*(Matrix& m, double f) {
	Matrix new_m(m.getRows(), m.getCols());
	for (int r = 0;r < m.getRows();r++) {
		for (int c = 0;c < m.getCols();c++) {
			new_m(r, c) = m(r, c) * f;
		}
	}
	return new_m;
}


Matrix operator*(const Matrix& m, double f) {
	Matrix new_m(m.getRows(), m.getCols());
	for (int r = 0;r < m.getRows();r++) {
		for (int c = 0;c < m.getCols();c++) {
			new_m(r, c) = m(r, c) * f;
		}
	}
	return new_m;
}



std::vector<double> operator*(Matrix& m, std::vector<double> vec) {

	std::vector<double> new_vec;
	new_vec.resize(m.getRows(), 0.f);

	for (int r = 0;r < m.getRows();r++) {
		for (int c = 0;c < m.getCols();c++) {
			new_vec[r] += vec[c] * m(r, c);
		}
	}
	return new_vec;
}


std::vector<double> operator*(const Matrix& m, std::vector<double> vec) {

	std::vector<double> new_vec;
	new_vec.resize(m.getRows(), 0.f);

	for (int r = 0;r < m.getRows();r++) {
		for (int c = 0;c < m.getCols();c++) {
			new_vec[r] += vec[c] * m(r, c);
		}
	}
	return new_vec;
}



