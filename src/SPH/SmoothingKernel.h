#pragma once

#include <vector>
#include "Particle.h"


class SmoothingKernel {
private:
	Particle* NeighbourParticle;

	float m_W;
	glm::vec3 m_dWdp;
	glm::vec3 m_d2Wdp2;
	float m_PairSmR = 0; //(Mpart_ptr->getSmR() + NBpart_ptr->getSmR()) / 2;
	float m_q = 0; //(Mpart_ptr->Getdr(*NBpart_ptr)) / m_PairSmR;
	float m_factor;



public:
	SmoothingKernel();




};