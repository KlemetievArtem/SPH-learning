#pragma once

#include "virtualComputationalDomain.h"
#include "Particle.h"
#include "FiniteParticleMethodFULL.h"

enum NBSearchAlgorithm {
	DIRECT
};
enum DISTRIBUTION {
	RANDOM,
	UNIFORM
};
enum TIME_INTEGRATION_SCHEME {
	EXPLICIT,
	IMPLICIT,
	SEMI_IMPICIT,
	SECOND_ORDER_SCEME,
	NONE
};
enum VAL_CHANGE_USING{
	VAL,
	DVAL_DT
};

enum GLOBAL_METHOD {
	STANDART_SPH,
	CORRECTIVE_SPH,
	FINITE_PARTICLE_METHOD,
	THE_ONLY_RIGHT
};

enum BOUNDARIES_TYPES {
	FIRSTORDER_FIRSTORDER,
	PERIODIC_PERIODIC,
	PERIODIC_FIRSTORDER,
	FIRSTORDER_PERIODIC
};

enum BOUNDARY_HANDLING {
	DUMMY_PARTICLES,
	MIRROR_PARTICLES,
	REPULSIVE_FORCE,
	RENORMALIZATION
};



struct ParticleRendererBuffer {
	part_prec_3 position;
	glm::vec3 color;
	part_prec size;
};


// Передавать одни и теже опции в SPH и в пары
struct SPH_OPTIONS {
	unsigned int nrOfParticles[ALLTYPES];
	NBSearchAlgorithm NBSAlg = DIRECT;
	DIMENSIONS KernelDimension = D1;

	part_prec smoothingKernelLengthCoefficient = 3.0;

	GLOBAL_METHOD SPH_algorithm = THE_ONLY_RIGHT;
	GLOBAL_METHOD SPH_iteration_algorithm = FINITE_PARTICLE_METHOD;


	VAL_CHANGE_USING densityChangeUsing = DVAL_DT; //VAL;
	int densityChangeAlgorithm = 0;
	int velocityChangeAlgorithm_pressurePart = 1;
	int velocityChangeAlgorithm_viscosityPart = 1;

	DISTRIBUTION distributionREAL = UNIFORM;//RANDOM;
	DISTRIBUTION distributionBOUNDARY = UNIFORM;
	TIME_INTEGRATION_SCHEME timeIntegrationScheme = SECOND_ORDER_SCEME;


	BOUNDARY_HANDLING boundary_handling = RENORMALIZATION;


	bool firstCycle = true;

	bool cornerVP = true;

	bool IsStab = false;

};

#include "ParticlePair.h"
#include "BoudaryMentor.h"

struct SPH_ESENTIALS {
	std::vector<Particle*> Particles;
	std::vector<ParticlePair*> ParticlePairs;
	~SPH_ESENTIALS() {
		for (auto*&i : Particles)
			delete i;
		for (auto*&i : ParticlePairs)
			delete i;
	}
};

class SPH_CD : public CompDomain {
private:
	SPH_ESENTIALS SPH;
	SPH_OPTIONS m_options;

	bool FirstSycle = true;

	float BH_D;
	float BH_r;
	BoundaryMentor BM;
	//std::unique_ptr<BoundaryMentor> BM;
	std::vector<ParticleRendererBuffer> PRB;


public:






	SPH_CD(SPH_OPTIONS options, DIMENSIONS dim) {
		setTypeTo(MICROCOSME::MC_SPH);
		nrOfDim = dim;
		m_options = options;
	}
	void Initilization();
	//void InitialRendering(std::vector<Mesh*>* meshes);

	void UpdateRendering(std::vector<Model*>* models);
	void AfterRendering(std::vector<Model*>* models);



	void singlestep();

	void timeStep(cd_prec dt);

	void neighbourSearch();
	void creatingVirtualParticles();
	void addingVirtualTOPairs();
	void deletingParticlePairs();

	void AfterRendering();


	void DensityVariation();
	void EvaluateEps();
	void InternalForces();

	void ExternalForces();

	void SaveMaxVelocity();

	void Boundary_Handling();

	//void CorrectiveKernelCalculation();

	void TAU(std::vector<glm::mat3x3>& tau, unsigned int m, unsigned int n);

	void Dissipation(std::vector<float>& dissip);



	void SPH_Iterations();
	//
	void TimeDerivative();
	void PPC_Density(ParticlePair* pp, std::vector<part_prec>* drhodt);
	void PPC_InternalForces(ParticlePair* pp, std::vector<part_prec_3>* dvdt);
	
	// В идеале
	//  function approximation                1-st derivative approximation           2-nd derivative approximation
	//
	//  Sum(Wij*Vj) = 1                       Sum(Wij,x*Vj) = 0                       Sum(Wij,xx*Vj) = 0
	//  Sum((xj-xi)*Wij*Vj) = 0				  Sum((xj-xi)*Wij,x*Vj) = 1				  Sum((xj-xi)*Wij,xx*Vj) = 0
	//  Sum((xj-xi)^2*Wij*Vj) = 0			  Sum((xj-xi)^2*Wij,x*Vj) = 0			  Sum((xj-xi)^2*Wij,xx*Vj) = 2
	//         ...							        ...								        ...
	//         ...							        ...								        ...
	//  Sum((xj-xi)^n*Wij*Vj) = 0			  Sum((xj-xi)^n*Wij,x*Vj) = 0             Sum((xj-xi)^n*Wij,xx*Vj) = 0
	//
	//FINITE PARTICLE METHOD
	//STANDART                                                                                                          в идеале матрица должна быть следующего вида:
	//  | Sum(fj*Wij*Vj)   |     | fi   |   | Sum(Wij*Vj)      Sum((xj-xi)*Wij*Vj)      Sum((yj-yi)*Wij*Vj)    |                      | 1    0   0  |
	//  | Sum(fj*Wij,x*Vj) |  =  | fi,x | * | Sum(Wij,x*Vj)    Sum((xj-xi)*Wij,x*Vj)    Sum((yj-yi)*Wij,x*Vj)  |					  | 0    1   ?  |
	//  | Sum(fj*Wij,y*Vj) |     | fi,y |   | Sum(Wij,y*Vj)    Sum((xj-xi)*Wij,y*Vj)    Sum((yj-yi)*Wij,y*Vj)  |					  | 0    ?   1  |
	void PPC_FPMStandart_General(ParticlePair* pp, std::vector<Matrix>* FPM_matrix);
	void PPC_FPMStandart_vec(ParticlePair* pp, std::vector<std::vector<part_prec>>* FPM_matrix, std::string param);

	//Difference
	//  | Sum((fj-fi)*Wij*Vj)   |     |   0  |   | 1    Sum((xj-xi)*Wij*Vj)      Sum((yj-yi)*Wij*Vj)    |
	//  | Sum((fj-fi)*Wij,x*Vj) |  =  | fi,x | * | 1    Sum((xj-xi)*Wij,x*Vj)    Sum((yj-yi)*Wij,x*Vj)  |
	//  | Sum((fj-fi)*Wij,y*Vj) |     | fi,y |   | 1    Sum((xj-xi)*Wij,y*Vj)    Sum((yj-yi)*Wij,y*Vj)  |
	void PPC_FPMDifference_General(ParticlePair* pp, std::vector<Matrix>* FPM_matrix);
	void PPC_FPMDifference_vec(ParticlePair* pp, std::vector<std::vector<part_prec>>* FPM_matrix, std::string param);

	//Difference + Trim
	//  | Sum((fj-fi)*Wij,x*Vj) |  =  | fi,x | * | 1    Sum((xj-xi)*Wij,x*Vj)    Sum((yj-yi)*Wij,x*Vj)  |
	//  | Sum((fj-fi)*Wij,y*Vj) |     | fi,y |   | 1    Sum((xj-xi)*Wij,y*Vj)    Sum((yj-yi)*Wij,y*Vj)  |
	void PPC_FPMDifTrim_General(ParticlePair* pp, std::vector<Matrix>* FPM_matrix);
	void PPC_FPMDifTrim_vec(ParticlePair* pp, std::vector<std::vector<part_prec>>* FPM_matrix, std::string param);


	//FOR CORRECTIVE SPH
	//  | Sum((fj-fi)*Wij,x*Vj) |  =  | fi,x | * | Sum((xj-xi)*Wij,x*Vj)  Sum((yj-yi)*Wij,x*Vj)  Sum((zj-zi)*Wij,x*Vj)  |
	//  | Sum((fj-fi)*Wij,y*Vj) |  =  | fi,y | * | Sum((xj-xi)*Wij,y*Vj)  Sum((yj-yi)*Wij,y*Vj)  Sum((zj-zi)*Wij,y*Vj)  |
	//  | Sum((fj-fi)*Wij,z*Vj) |  =  | fi,z | * | Sum((xj-xi)*Wij,z*Vj)  Sum((yj-yi)*Wij,z*Vj)  Sum((zj-zi)*Wij,z*Vj)  |
	void PPC_CorSPH_General(ParticlePair* pp, std::vector<Matrix>* FPM_matrix);
	void PPC_CorSPH_vec(ParticlePair* pp, std::vector<std::vector<part_prec>>* FPM_matrix, std::string param);





	void DensityCalculation(std::vector<part_prec>* drhodt, ParticlePair* pp);
	void PressurePartCalculation(std::vector<part_prec_3>* dvdt, ParticlePair* pp);
	void ViscosityPartCalculation(std::vector<part_prec_3>* dvdt, ParticlePair* pp);

	void DensityRecalculation();
	void VelocityRecalculation();
	void DensityAndVelocityRecalculation();
	void RecalcPrep();


	void RenormFactorCalc(std::vector<part_prec>* gamma, ParticlePair* pp);
	void RenormFactorDerivCalc(std::vector<part_prec_3>* gamma_deriv, ParticlePair* pp);
	void RenormPressurePartCalculation(std::vector<part_prec_3>* dvdt, ParticlePair* pp, std::vector<part_prec>* gamma);
	

	



	~SPH_CD();

	//COLORING
	void Coloring();
	void ColoringBtType();

	void punctualColorChange(int number, glm::vec3 color);

	void PRB_refresh();


	int getNrOfParticles(PARTICLETYPE type);
	float getGlobalStats(int number);
	std::string getLocalStats(int number);
};
