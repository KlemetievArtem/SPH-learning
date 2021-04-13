#pragma once
#include "libs.h"


enum DIMENSIONS {
	D1 = 1,
	D2,
	D3
};

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
	float fi;
	float fix;
	float fiy;
	float fiz;

	glm::vec2 f_1D;/* = */glm::mat2x2 Sum_W_dx_1D;/* X */glm::vec2 Sum_f_W_1D;

	glm::vec3 f_2D;/* = */glm::mat3x3 Sum_W_dx_2D;/* X */glm::vec3 Sum_f_W_2D;

	glm::vec4 f_3D;/* = */glm::mat4x4 Sum_W_dx_3D;/* X */glm::vec4 Sum_f_W_3D;


public:
	FiniteParticleMethod(DIMENSIONS dim) : m_dimension(dim) {
		switch (dim) {
		case(1):
			CreatingMatrix_1D();
			CreatingVectors_1D();
			break;
		case(2):
			CreatingMatrix_2D();
			CreatingVectors_2D();
			break;
		case(3):
			CreatingMatrix_3D();
			CreatingVectors_3D();
			break;
		default:
			break;
		}
	};

/*
	void Adding(glm::mat4x4 mat4, glm::vec4 vec4) {
		switch (m_dimension) {
		case(1):
			glm::mat2x2 mat2;
			mat2[0, 0] = mat4[0, 0];	mat2[0, 1] = mat4[0, 1];
			mat2[1, 0] = mat4[1, 0];	mat2[1, 1] = mat4[1, 1];
			glm::vec2 vec2;
			vec2[0] = vec4[0]; vec2[1] = vec4[1];
			AddingToMatrix_1D(mat2);
			AddingToVector_1D(vec2);
			break;
		case(2):
			glm::mat3x3 mat3;
			mat3[0, 0] = mat4[0, 0];	mat3[0, 1] = mat4[0, 1];	mat3[0, 2] = mat4[0, 2];
			mat3[1, 0] = mat4[1, 0];	mat3[1, 1] = mat4[1, 1];	mat3[1, 2] = mat4[1, 2];
			mat3[2, 0] = mat4[2, 0];	mat3[2, 1] = mat4[2, 1];	mat3[2, 2] = mat4[2, 2];
			glm::vec3 vec3;
			vec3[0] = vec4[0]; vec3[1] = vec4[1]; vec3[2] = vec4[2];
			AddingToMatrix_2D(mat3);
			AddingToVector_2D(vec3);
			break;
		case(3):
			AddingToMatrix_3D(mat4);
			AddingToVector_3D(vec4);
			break;
		default:
			break;
		}

	}
*/


	void Refresh() {
		switch (m_dimension) {
		case(1):
			CreatingMatrix_1D();
			CreatingVectors_1D();
			break;
		case(2):
			CreatingMatrix_2D();
			CreatingVectors_2D();
			break;
		case(3):
			CreatingMatrix_3D();
			CreatingVectors_3D();
			break;
		default:
			break;
		}
	}

	void CreatingMatrix_1D() { Sum_W_dx_1D = glm::mat2x2(0.f); }
	void CreatingMatrix_2D() { Sum_W_dx_2D = glm::mat3x3(0.f); }
	void CreatingMatrix_3D() { Sum_W_dx_3D = glm::mat4x4(0.f); }

	void CreatingVectors_1D() { f_1D = glm::vec2(0.f); Sum_f_W_1D = glm::vec2(0.f); }
	void CreatingVectors_2D() { f_2D = glm::vec3(0.f); Sum_f_W_2D = glm::vec3(0.f); }
	void CreatingVectors_3D() { f_3D = glm::vec4(0.f); Sum_f_W_3D = glm::vec4(0.f); }

//	void AddingToMatrix_1D(glm::mat2x2 add_1D) { Sum_W_dx_1D += add_1D; }
//	void AddingToMatrix_2D(glm::mat3x3 add_2D) { Sum_W_dx_2D += add_2D; }
//	void AddingToMatrix_3D(glm::mat4x4 add_3D) { Sum_W_dx_3D += add_3D; }

//	void AddingToVector_1D(glm::vec2 add_1D) { Sum_f_W_1D += add_1D; }
//	void AddingToVector_2D(glm::vec3 add_2D) { Sum_f_W_2D += add_2D; }
//	void AddingToVector_3D(glm::vec4 add_3D) { Sum_f_W_3D += add_3D; }


	glm::vec4 Results(glm::mat4x4 mat4, glm::vec4 vec4) {
		switch (m_dimension) {
		case(1):{
			glm::mat2x2 mat2;
			mat2[0][0] = mat4[0][0];	mat2[0][1] = mat4[0][1];
			mat2[1][0] = mat4[1][0];	mat2[1][1] = mat4[1][1];
			Sum_W_dx_1D = mat2;
			glm::vec2 vec2;
			vec2[0] = vec4[0]; vec2[1] = vec4[1];
			Sum_f_W_1D = vec2;

			glm::vec2 retVal = FPM_1D_result();
			return glm::vec4(retVal.x, retVal.y, 0.f, 0.f);
			break;
		}
		case(2):{
			glm::mat3x3 mat3;
			mat3[0][0] = mat4[0][0];	mat3[0][1] = mat4[0][1];	mat3[0][2] = mat4[0][2];
			mat3[1][0] = mat4[1][0];	mat3[1][1] = mat4[1][1];	mat3[1][2] = mat4[1][2];
			mat3[2][0] = mat4[2][0];	mat3[2][1] = mat4[2][1];	mat3[2][2] = mat4[2][2];
			Sum_W_dx_2D = mat3;
			glm::vec3 vec3;
			vec3[0] = vec4[0]; vec3[1] = vec4[1]; vec3[2] = vec4[2];
			Sum_f_W_2D = vec3;

			glm::vec3 retVal = FPM_2D_result();
			return glm::vec4(retVal[0], retVal[1], retVal[2], 0.f);
			break;
		}
		case(3):{
			Sum_W_dx_3D = mat4;
			Sum_f_W_3D = vec4;
			return FPM_3D_result();
			break;
		}
		default:
			break;
		}
	}


	glm::vec2 FPM_1D_result() { return glm::inverse(Sum_W_dx_1D) * Sum_f_W_1D; }
	glm::vec3 FPM_2D_result() {
		//std::cout << Sum_W_dx_2D[0][0] << "  " << Sum_W_dx_2D[0][1] << "  " << Sum_W_dx_2D[0][2] <<"-1   " << Sum_f_W_3D[0] << "\n";
		//std::cout << Sum_W_dx_2D[1][0] << "  " << Sum_W_dx_2D[1][1] << "  " << Sum_W_dx_2D[1][2] <<"   x " << Sum_f_W_3D[1] << "\n";
		//std::cout << Sum_W_dx_2D[2][0] << "  " << Sum_W_dx_2D[2][1] << "  " << Sum_W_dx_2D[2][2] <<"     " << Sum_f_W_3D[2] << "\n";
		//std::cout << "\n";
		return glm::inverse(Sum_W_dx_2D) * Sum_f_W_2D; 
	}
	glm::vec4 FPM_3D_result() { return glm::inverse(Sum_W_dx_3D) * Sum_f_W_3D; }


};
