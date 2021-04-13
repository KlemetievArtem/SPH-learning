#pragma once



class FiniteParticleMethod {

//      f                                Sum_W_dx                                                                     Sum_f_W
//                                                                                                            -1
//   | fi   |      | sum(Wij*Vj)     sum((xj-xi)*Wij*Vj)     sum((yj-yi)*Wij*Vj)      sum((zj-zi)*Wij*Vj)    |   | sum(fj*Wij*Vj)      |
//   | fi,x |   =  | sum(Wij,x*Vj)   sum((xj-xi)*Wij,x*Vj)   sum((yj-yi)*Wij,x*Vj)    sum((zj-zi)*Wij,x*Vj)  | X | sum(fj*Wij,x*Vj)    |
//   | fi,y |      | sum(Wij,y*Vj)   sum((xj-xi)*Wij,y*Vj)   sum((yj-yi)*Wij,y*Vj)    sum((zj-zi)*Wij,y*Vj)  |   | sum(fj*Wij,y*Vj)    |
//   | fi,z |      | sum(Wij,z*Vj)   sum((xj-xi)*Wij,z*Vj)   sum((yj-yi)*Wij,z*Vj)    sum((zj-zi)*Wij,z*Vj)  |   | sum(fj*Wij,z*Vj)    |
//
//


private:
	int m_dimension;

	std::vector<part_prec>  fi;
	//                          =    -1
	Matrix        Sum_W_dx;
	//                          * 
	std::vector<part_prec>   Sum_fj_W;


	std::vector<std::vector<DIMENSIONS>> dx;








public:
	FiniteParticleMethod(DIMENSIONS dim) : m_dimension(dim){
		int m_size = 0 ;
		for (int i = 1;i <= dim + 1;i++) {
			m_size += i;
		}

		Sum_fj_W.resize(m_size);
		Sum_W_dx.resize(m_size, m_size);

		fi.resize(m_size);

		dx.resize(m_size);
		for (int i = 0;i < m_size;i++) {
			switch (m_dimension)
			{
			case(1):
				switch (i) {
				case(0): dx[i].push_back(D0); break;
				case(1): dx[i].push_back(D1); break;
				case(2): dx[i].push_back(D1); dx[i].push_back(D1); break;
				default: break;
				}
				break;
			case(2):
				switch (i) {
				case(0): dx[i].push_back(D0); break;
				case(1): dx[i].push_back(D1); break;
				case(2): dx[i].push_back(D2); break;
				case(3): dx[i].push_back(D1); dx[i].push_back(D1); break;
				case(4): dx[i].push_back(D2); dx[i].push_back(D2); break;
				case(5): dx[i].push_back(D1); dx[i].push_back(D2); break;
				default: break;
				}
				break;
			case(3):
				switch (i) {
				case(0): dx[i].push_back(D0); break;
				case(1): dx[i].push_back(D1); break;
				case(2): dx[i].push_back(D2); break;
				case(3): dx[i].push_back(D3); break;
				case(4): dx[i].push_back(D1); dx[i].push_back(D1); break;
				case(5): dx[i].push_back(D2); dx[i].push_back(D2); break;
				case(6): dx[i].push_back(D3); dx[i].push_back(D3); break;
				case(7): dx[i].push_back(D1); dx[i].push_back(D2); break;
				case(8): dx[i].push_back(D1); dx[i].push_back(D3); break;
				case(9): dx[i].push_back(D2); dx[i].push_back(D3); break;
				default: break;
				}
				break;
			default:
				break;
			}
		}
	}
	
	//FiniteParticleMethod(): fi(std::vector<float>()), Sum_W_dx(Matrix()), Sum_fj_W(std::vector<float>()),dx(std::vector<std::vector<DIMENSIONS>>()) {}


	~FiniteParticleMethod() {
	//	delete Sum_fj_W;
	}





	void addDimensions(DIMENSIONS dim) {
		m_dimension = dim;
		int m_size = 0;
		for (int i = 1;i <= dim + 1;i++) {
			m_size += i;
		}

		Sum_fj_W.resize(m_size);
		Sum_W_dx.resize(m_size, m_size);
		fi.resize(m_size);

		dx.resize(m_size);
		for (int i = 0;i < m_size;i++) {
			switch (m_dimension)
			{
			case(1):
				switch (i) {
				case(0): dx[i].push_back(D0); break;
				case(1): dx[i].push_back(D1); break;
				case(2): dx[i].push_back(D1); dx[i].push_back(D1); break;
				default: break;
				}
				break;
			case(2):
				switch (i) {
				case(0): dx[i].push_back(D0); break;
				case(1): dx[i].push_back(D1); break;
				case(2): dx[i].push_back(D2); break;
				case(3): dx[i].push_back(D1); dx[i].push_back(D1); break;
				case(4): dx[i].push_back(D2); dx[i].push_back(D2); break;
				case(5): dx[i].push_back(D1); dx[i].push_back(D2); break;
				default: break;
				}
				break;
			case(3):
				switch (i) {
				case(0): dx[i].push_back(D0); break;
				case(1): dx[i].push_back(D1); break;
				case(2): dx[i].push_back(D2); break;
				case(3): dx[i].push_back(D3); break;
				case(4): dx[i].push_back(D1); dx[i].push_back(D1); break;
				case(5): dx[i].push_back(D2); dx[i].push_back(D2); break;
				case(6): dx[i].push_back(D3); dx[i].push_back(D3); break;
				case(7): dx[i].push_back(D1); dx[i].push_back(D2); break;
				case(8): dx[i].push_back(D1); dx[i].push_back(D3); break;
				case(9): dx[i].push_back(D2); dx[i].push_back(D3); break;
				default: break;
				}
				break;
			default:
				break;
			}
		}

	}


	void addMatrix(Matrix mat) {
		//Sum_W_dx.lateInitilization(mat.getRows(), mat.getCols());
		//std::cout << mat << "\n";
		Sum_W_dx = mat;
	}
	void addVec(std::vector<part_prec> vec) {
		//Sum_fj_W.assign(vec.begin(), vec.begin() + vec.size());
		//Sum_fj_W.resize(vec.size());
		//std::cout << Sum_fj_W.size() << "\n";
		//std::cout << "addVec " << Sum_fj_W << "  " << vec << "\n";
		Sum_fj_W = vec; 
	}

	void addMatVec(Matrix mat, std::vector<part_prec> vec) { addVec(vec); addMatrix(mat); }

	void Calculate() {

		if (Sum_W_dx.filedWithZeros()) {
			for (auto i : fi) {
				i = 0.f;
			}
		}
		else {
			fi = inverse(Sum_W_dx)*(Sum_fj_W);
		}

		/*
		Matrix inv = inverse(Sum_W_dx);
		
		for (int r = 0; r < fi.size(); r++)	{
			std::cout << "| " << fi[r] << " | " << " | ";
			for (int c = 0; c < fi.size(); c++) {
				std::cout << inv(r, c) << "  ";
			}
			std::cout << "| " << " | " << Sum_fj_W[r] << " |\n";
		}
		*/
	}


	std::vector<part_prec>* Results() {
		//fi.resize(Sum_fj_W.size());
		fi = inverse(Sum_W_dx)*(Sum_fj_W);
		return &fi;
	}

	part_prec getd(DIMENSIONS d1) {
		assert(m_dimension >= d1);
		switch (d1) {
		case(D0): return fi[0]; break;
		case(D1): return fi[1]; break;
		case(D2): return fi[2]; break;
		case(D3): return fi[3]; break;
		default:
			break;
		}
		return 0;
	}

	part_prec getd(DIMENSIONS d1, DIMENSIONS d2) {
		assert(m_dimension >= d1);
		assert(m_dimension >= d2);
		switch (m_dimension) {
		case(D1):
			if ((d1 == D1) and (d2 == D1)) { return fi[2]; }
			break;
		case(D2):
			if ((d1 == D1) and (d2 == D1)) { return fi[3]; }
			if ((d1 == D2) and (d2 == D2)) { return fi[4]; }
			if ((d1 == D1) and (d2 == D2)) { return fi[5]; }
			break;
		case(D3):
			if ((d1 == D1) and (d2 == D1)) { return fi[4]; }
			if ((d1 == D2) and (d2 == D2)) { return fi[5]; }
			if ((d1 == D3) and (d2 == D3)) { return fi[6]; }
			if ((d1 == D1) and (d2 == D2)) { return fi[7]; }
			if ((d1 == D1) and (d2 == D3)) { return fi[8]; }
			if ((d1 == D2) and (d2 == D3)) { return fi[9]; }
			break;
		default:
			break;
		}
		return 0;
	}

};
