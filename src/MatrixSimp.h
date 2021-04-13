#pragma once

#include <iostream>
#include <vector>
#include <fstream>
#include <cassert>


class Matrix {
private:
	int m_rows;
	int m_cols;
	std::vector<std::vector<double>> data;
public:
	Matrix(int r, int c) : m_rows(r), m_cols(c) {

		data.resize(r);
		for (int r_c = 0;r_c < r;r_c++) {
			for (int c_c = 0;c_c < c;c_c++) {
				data[r_c].push_back(0);
			}
		}
	}
	Matrix() : m_rows(0), m_cols(0) {}

	Matrix& operator=(Matrix m) {
		if (&m == this) {
			return *this;
		}
		resize(m.getRows(), m.getCols());
		for (int r = 0;r < m.getRows();r++) {
			for (int c = 0;c < m.getCols();c++) {
				data[r][c] = m(r, c);
			}

		}
	}




	void clear() {
		if (data.size() != 0) {
			for (auto i : data) {
				if (i.size() != 0) {
					i.resize(0);
				}
			}
			data.resize(0);
		}
	}

	void resize(int r, int c) {
		if(m_rows!=0){
			clear();
		}
		m_rows = r;
		m_cols = c;
		
		data.resize(r);
		for (auto &i : data) {
			i.resize(c);
		}
	}

	double& operator()(int r, int c) {
		assert(c >= 0 && c <= m_cols);
		assert(r >= 0 && r <= m_rows);
		return data[r][c];
	}

	const double& operator()(int r, int c) const {
		assert(c >= 0 && c <= m_cols);
		assert(r >= 0 && r <= m_rows);

		return data[r][c];
	}

	const int getRows() const {	return m_rows; }
	const int getCols() const { return m_cols; }

	void setRows(int r) { 
		m_rows = r;

	}
	void setCols(int c) { 
		m_cols = c;
	}
	void lateInitilization(int r, int c) {
		setRows(r); setCols(c);
		data.resize(m_rows);
		for (int r_c = 0;r_c < m_rows;r_c++) {
			for (int c_c = 0;c_c < m_cols;c_c++) {
				data[r_c].push_back(0);
			}
		}
	}



	const Matrix copyAndCut(int r_tc, int c_tc) const {
		int new_rows = m_rows - 1;
		int new_cols = m_cols - 1;

		Matrix new_m(new_rows, new_cols);


		new_rows = 0;
		new_cols = 0;
		for (int r = 0;r < m_rows;r++) {
			if (r != r_tc) {
				new_cols = 0;
				for (int c = 0;c < m_cols;c++) {
					if (c != c_tc) {
						new_m(new_rows, new_cols) = (*this)(r,c);
					}
					else {
						new_cols--;
					}
					new_cols++;
				}
			}
			else {
				new_rows--;
			}
			new_rows++;
		}
		return new_m;

	}
	



	bool filedWithZeros() {
		for (int r = 0;r < m_rows;r++) {
			for (int c = 0;c < m_cols;c++) {
				if (data[r][c] != 0) {
					return false;
				}
			}
		}

		return true;
	}



	friend std::ostream& operator<<(std::ostream& out, const Matrix& m);

	friend Matrix operator*(Matrix& m, double f);
	friend Matrix operator*(const Matrix& m, double f);

	friend std::vector<double> operator*(Matrix& m, std::vector<double> vec);
	friend std::vector<double> operator*(const Matrix& m, std::vector<double> vec);



	Matrix& operator+=(Matrix& m) {
		assert(m.getCols() >= 0 && m.getCols() <= m_cols);
		assert(m.getRows() >= 0 && m.getRows() <= m_rows);

		for (int r = 0;r < m_rows;r++) {
			for (int c = 0;c < m_cols;c++) {
				(this->data)[r][c] += m(r, c);
			}
		}
		return *this;
	}




	friend double determinant(const Matrix& m);
	friend const Matrix inverse(const Matrix& m);
	friend Matrix inverse( Matrix& m);



};





