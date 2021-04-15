#pragma once
#include "Particle.h"
//#include "Kernel.h"

enum DiffAxis {
	X,
	Y,
	Z,
	R
};



class Kernel {
public:
	virtual part_prec_3 getFactor(part_prec SmR) = 0;
	virtual part_prec W(part_prec Rds) = 0;
	virtual part_prec dW(part_prec Rds, part_prec SmR) = 0;
	virtual part_prec d2W(part_prec Rds, part_prec SmR) = 0;
	virtual part_prec d2Wdr2(part_prec Rds, part_prec SmR) = 0;
	virtual part_prec d2W_add(part_prec Rds, part_prec SmR) = 0;
};


class QubicSpline : public Kernel {
public:
	QubicSpline() {}
	part_prec_3 getFactor(part_prec SmR) override {
		return { (1.0 / SmR), (15.0 / (7.0 * M_PI*pow(SmR, 2))), (3.0 / (2.0 * M_PI*pow(SmR, 3))) };
	}
	part_prec W(part_prec Rds) override {
		part_prec retVal;
		if (Rds >= 0 and Rds < 1) retVal = ((2.0 / 3.0 - pow(Rds, 2) + pow(Rds, 3) / 2.0));
		if (Rds >= 1 and Rds < 2) retVal = (1.0 / 6.0 * pow((2.0 - Rds), 3));
		if (Rds >= 2)             retVal = 0.0;
		return retVal;
	}
	part_prec dW(part_prec Rds, part_prec SmR) override {
		part_prec retVal;
		if (Rds >= 0 and Rds < 1) retVal = (-2.0*Rds / SmR + 3.0 / 2.0*pow(Rds, 2) / SmR);
		if (Rds >= 1 and Rds < 2) retVal = -1.0 / 2.0 * pow((2.0 - Rds), 2) / SmR;
		if (Rds >= 2)             retVal = 0.0;
		return retVal;
	}
	part_prec d2W(part_prec Rds, part_prec SmR) override {
		part_prec retVal;
		if (Rds >= 0 and Rds < 1) retVal = (3.0 / 2.0*Rds / pow(SmR, 2));
		if (Rds >= 1 and Rds < 2) retVal = (1.f / 2.f*(2.f*(2.f - Rds) / pow(SmR, 2) + pow((2.f - Rds), 2) / (SmR*SmR*Rds)));
		if (Rds >= 2)             retVal = 0.0;
		return retVal;
	}
	part_prec d2Wdr2(part_prec Rds, part_prec SmR) override {
		part_prec retVal;
		if (Rds >= 0 and Rds < 1) retVal = (-2.0 / pow(SmR, 2) + 3.0 * Rds / pow(SmR, 2));
		if (Rds >= 1 and Rds < 2) retVal = (2.f - Rds) / pow(SmR, 2);
		if (Rds >= 2)             retVal = 0.0;
		return retVal;
	}
	part_prec d2W_add(part_prec Rds, part_prec SmR) override {
		part_prec retVal;
		if (Rds >= 0 and Rds < 1) retVal = (-2.0*Rds / (SmR*SmR*Rds) + 3.0 / 2.0*pow(Rds, 2) / SmR * SmR*Rds);
		if (Rds >= 1 and Rds < 2) retVal = pow(2.0 - Rds, 2) / (SmR * SmR*Rds);
		if (Rds >= 2)             retVal = 0.0;
		return retVal;
	}
};
class QuinticSpline : public Kernel {
public:
	QuinticSpline() {}
	part_prec_3 getFactor(part_prec SmR) override {
		return { (120.0 / SmR), (7.0 / (478.0 * M_PI*pow(SmR, 2))), (3.0 / (359.0 * M_PI*pow(SmR, 3))) };
	}
	part_prec W(part_prec Rds) override {
		part_prec retVal;
		if (Rds >= 0 and Rds < 1) retVal = (pow(3 - Rds, 5) - 6.0*pow(2 - Rds, 5) + 15.0*pow(1 - Rds, 5));
		if (Rds >= 1 and Rds < 2) retVal = (pow(3 - Rds, 5) - 6.0*pow(2 - Rds, 5));
		if (Rds >= 2 and Rds < 3) retVal = (pow(3 - Rds, 5));
		if (Rds >= 3)             retVal = 0.0;
		return retVal;
	}
	part_prec dW(part_prec Rds, part_prec SmR) override {
		part_prec retVal;
		if (Rds >= 0 and Rds < 1) retVal = (-5.0/SmR*pow(3 - Rds, 4) + 30.0 / SmR * pow(2 - Rds, 4) - 75.0 / SmR * pow(1 - Rds, 4));
		if (Rds >= 1 and Rds < 2) retVal = (-5.0 / SmR * pow(3 - Rds, 4) + 30.0 / SmR * pow(2 - Rds, 4));
		if (Rds >= 2 and Rds < 3) retVal = (-5.0 / SmR * pow(3 - Rds, 4));
		if (Rds >= 3)             retVal = 0.0;
		return retVal;
	}
	part_prec d2W(part_prec Rds, part_prec SmR) override {
		part_prec retVal;
		if (Rds >= 0 and Rds < 1) retVal = (5.0 / Rds/ pow(SmR, 2) * pow(3 - Rds, 4) + 20.0 / pow(SmR, 2) * pow(3 - Rds, 3) - 30.0 / Rds / pow(SmR, 2) * pow(2 - Rds, 4) - 120.0 / pow(SmR, 2) * pow(2 - Rds, 3) + 75.0 / Rds / pow(SmR, 2) * pow(1 - Rds, 4) + 300.0 / pow(SmR, 2) * pow(1 - Rds, 3));
		if (Rds >= 1 and Rds < 2) retVal = (5.0 / Rds / pow(SmR, 2) * pow(3 - Rds, 4) + 20.0 / pow(SmR, 2) * pow(3 - Rds, 3) - 30.0 / Rds / pow(SmR, 2) * pow(2 - Rds, 4) - 120.0 / pow(SmR, 2) * pow(2 - Rds, 3));
		if (Rds >= 2 and Rds < 3) retVal = (5.0 / Rds / pow(SmR, 2) * pow(3 - Rds, 4) + 20.0 / pow(SmR, 2) * pow(3 - Rds, 3));
		if (Rds >= 3)             retVal = 0.0;
		return retVal;
	}
	part_prec d2Wdr2(part_prec Rds, part_prec SmR) override {
		part_prec retVal;
		if (Rds >= 0 and Rds < 1) retVal = (20.0 / pow(SmR,2) * pow(3 - Rds, 3) - 120.0 / pow(SmR, 2) * pow(2 - Rds, 3) + 300.0 / pow(SmR, 2) * pow(1 - Rds, 3));
		if (Rds >= 1 and Rds < 2) retVal = (20.0 / pow(SmR, 2) * pow(3 - Rds, 3) - 120.0 / pow(SmR, 2) * pow(2 - Rds, 3));
		if (Rds >= 2 and Rds < 3) retVal = (20.0 / pow(SmR, 2) * pow(3 - Rds, 3));
		if (Rds >= 3)             retVal = 0.0;
		return retVal;
	}
	part_prec d2W_add(part_prec Rds, part_prec SmR) override {
		part_prec retVal;
		if (Rds >= 0 and Rds < 1) retVal = (-5.0 / Rds / pow(SmR, 2) * pow(3 - Rds, 4) + 30.0 / Rds / pow(SmR, 2) * pow(2 - Rds, 4) - 75.0 / Rds / pow(SmR, 2) * pow(1 - Rds, 4));
		if (Rds >= 1 and Rds < 2) retVal = (-5.0 / Rds / pow(SmR, 2) * pow(3 - Rds, 4) + 30.0 / Rds / pow(SmR, 2) * pow(2 - Rds, 4));
		if (Rds >= 2 and Rds < 3) retVal = (-5.0 / Rds / pow(SmR, 2) * pow(3 - Rds, 4));
		if (Rds >= 3)             retVal = 0.0;
		return retVal;
	}
};
class Gaussian : public Kernel {
public:
	Gaussian() {}
	part_prec_3 getFactor(part_prec SmR) override {
		return { (1.0 /(pow(M_PI,0.5)* SmR)), (1.0 / (M_PI*pow(SmR, 2))), (1.0 / (pow(M_PI,1.5)*pow(SmR, 3))) };
	}
	part_prec W(part_prec Rds) override {
		part_prec retVal = exp(-pow(Rds,2));
		return retVal;
	}
	part_prec dW(part_prec Rds, part_prec SmR) override {
		part_prec retVal = exp(-pow(Rds, 2))*(-2.0*Rds/SmR);
		return retVal;
	}
	part_prec d2W(part_prec Rds, part_prec SmR) override {
		part_prec retVal = exp(-pow(Rds, 2))*(4.0*pow(Rds/ SmR,2));
		return retVal;
	}
	part_prec d2Wdr2(part_prec Rds, part_prec SmR) override {
		part_prec retVal = exp(-pow(Rds, 2))*(4.0*pow(Rds / SmR, 2) - 2.0 / pow(SmR, 2));
		return retVal;
	}
	part_prec d2W_add(part_prec Rds, part_prec SmR) override {
		part_prec retVal = exp(-pow(Rds, 2))*( -2.0/pow(SmR,2));
		return retVal;
	}
};



class ParticlePair {
private:
	part_prec m_aveSmR;
	part_prec m_Rds;
	part_prec m_distance;
	//SPH_OPTIONS m_options;
	DIMENSIONS m_dim;
	Kernel* mode;
public:
	int m_PairId; // Номер пары 
	Particle* Mpart_ptr;
	Particle* NBpart_ptr;
	void ChangeKernel(Kernel* k) { mode = k; }
	void saveAveSmR() { m_aveSmR = (Mpart_ptr->m_SmR + NBpart_ptr->m_SmR) / 2.0; }
	void saveDistance() { m_distance = Mpart_ptr->distance(NBpart_ptr); }
	void saveDimensionlessR() { m_Rds = m_distance / m_aveSmR; }

	part_prec calcFactor(DIMENSIONS SpecDim = D0) {
		part_prec_3 factor_3 = mode->getFactor(m_aveSmR);
		if (SpecDim == D0) {
			switch (m_dim)
			{
			case 1:
				return factor_3[0];
				break;
			case 2:
				return factor_3[1];
				break;
			case 3:
				return factor_3[2];
				break;
			default:
				assert("ParticlePair::CalcFactor" && 0);
				break;
			}
		}
		else {
			return factor_3[SpecDim-1];
		}
	}

	part_prec W(Particle* part_ptr) {
		part_prec factor = calcFactor();
		part_prec retVal;
		retVal = mode->W(m_Rds);;
		return retVal * factor;
	}
	part_prec dWd(DiffAxis axis, Particle* part_ptr) {
		part_prec factor = calcFactor();
		part_prec dWd;
		dWd = mode->dW(m_Rds, m_aveSmR);		
		switch (axis){
		case X:
			if (part_ptr == Mpart_ptr)   dWd *= part_ptr->dx(NBpart_ptr) / m_distance * factor;
			else if (part_ptr == NBpart_ptr)   dWd *= part_ptr->dx(Mpart_ptr) / m_distance * factor;
			break;
		case Y:
			if (part_ptr == Mpart_ptr)   dWd *= part_ptr->dy(NBpart_ptr) / m_distance * factor;
			else if (part_ptr == NBpart_ptr)   dWd *= part_ptr->dy(Mpart_ptr) / m_distance * factor;
			break;
		case Z:
			if (part_ptr == Mpart_ptr)   dWd *= part_ptr->dz(NBpart_ptr) / m_distance * factor;
			else if (part_ptr == NBpart_ptr)   dWd *= part_ptr->dz(Mpart_ptr) / m_distance * factor;
			break;
		case R:
			if (part_ptr == Mpart_ptr) dWd= dWd* calcFactor(D2);
			else if (part_ptr == NBpart_ptr) dWd= dWd* calcFactor(D2);
			
			break;
		default:
			assert("ParticlePair::dWd" && 0);
			break;
		}
		return dWd;
	}
	part_prec d2Wd(DiffAxis axis1, DiffAxis axis2, Particle* part_ptr) {
		part_prec factor = calcFactor();
		part_prec d2Wd, d2Wd_add;
		d2Wd = mode->d2W(m_Rds, m_aveSmR);
		d2Wd_add = mode->d2W_add(m_Rds, m_aveSmR);
		if ((axis1 == axis2)) {
			switch (axis1) {
			case X:
				if (part_ptr == Mpart_ptr) d2Wd *= part_ptr->dx(NBpart_ptr) / m_distance * part_ptr->dx(NBpart_ptr) / m_distance;
				else if (part_ptr == NBpart_ptr) d2Wd *= part_ptr->dx(Mpart_ptr) / m_distance * part_ptr->dx(Mpart_ptr) / m_distance;
				d2Wd += d2Wd_add;
				d2Wd *= factor;
				break;
			case Y:
				if (part_ptr == Mpart_ptr)  d2Wd *= part_ptr->dy(NBpart_ptr) / m_distance * part_ptr->dy(NBpart_ptr) / m_distance;
				else if (part_ptr == NBpart_ptr)  d2Wd *= part_ptr->dy(Mpart_ptr) / m_distance * part_ptr->dy(Mpart_ptr) / m_distance;
				d2Wd += d2Wd_add;
				d2Wd *= factor;
				break;
			case Z:
				if (part_ptr == Mpart_ptr)  d2Wd *= part_ptr->dz(NBpart_ptr) / m_distance * part_ptr->dz(NBpart_ptr) / m_distance;
				else if (part_ptr == NBpart_ptr)  d2Wd *= part_ptr->dz(Mpart_ptr) / m_distance * part_ptr->dz(Mpart_ptr) / m_distance;
				d2Wd += d2Wd_add;
				d2Wd *= factor;
				break;
			case R:
				if (part_ptr == Mpart_ptr)  d2Wd = mode->d2Wdr2(m_Rds, m_aveSmR)* calcFactor();
				else if (part_ptr == NBpart_ptr)  d2Wd = mode->d2Wdr2(m_Rds, m_aveSmR)* calcFactor();
				
				break;
			default:
				assert("ParticlePair::d2Wd" && 0);
				break;
			}
		}
		if (((axis1 == X) and (axis2 == Y)) or ((axis1 == Y) and (axis2 == X))) {
			if (part_ptr == Mpart_ptr)  d2Wd *= part_ptr->dx(NBpart_ptr) / m_distance * part_ptr->dy(NBpart_ptr) / m_distance * factor;
			else if (part_ptr == NBpart_ptr)  d2Wd *= part_ptr->dx(Mpart_ptr) / m_distance * part_ptr->dy(Mpart_ptr) / m_distance * factor;
		}
		if (((axis1 == X) and (axis2 == Z)) or ((axis1 == Z) and (axis2 == X))) {
			if (part_ptr == Mpart_ptr)  d2Wd *= part_ptr->dx(NBpart_ptr) / m_distance * part_ptr->dz(NBpart_ptr) / m_distance * factor;
			else if (part_ptr == NBpart_ptr)  d2Wd *= part_ptr->dx(Mpart_ptr) / m_distance * part_ptr->dz(Mpart_ptr) / m_distance * factor;
		}
		if (((axis1 == Y) and (axis2 == Z)) or ((axis1 == Y) and (axis2 == X))) {
			if (part_ptr == Mpart_ptr)  d2Wd *= part_ptr->dy(NBpart_ptr) / m_distance * part_ptr->dz(NBpart_ptr) / m_distance * factor;
			else if (part_ptr == NBpart_ptr)  d2Wd *= part_ptr->dy(Mpart_ptr) / m_distance * part_ptr->dz(Mpart_ptr) / m_distance * factor;
		}
		return d2Wd;
	}


	inline part_prec_3 vec_e(Particle* part_ptr) {
		if (part_ptr == Mpart_ptr) return part_prec_3{ part_ptr->dx(NBpart_ptr) / (m_distance*(1 + m_aveSmR*0.01)),part_ptr->dy(NBpart_ptr) / (m_distance*(1 + m_aveSmR * 0.01)) ,part_ptr->dz(NBpart_ptr) / (m_distance*(1 + m_aveSmR * 0.01)) };
		else if (part_ptr == NBpart_ptr) return part_prec_3{ part_ptr->dx(Mpart_ptr) / (m_distance*(1 + m_aveSmR * 0.01)),part_ptr->dy(Mpart_ptr) / (m_distance*(1 + m_aveSmR * 0.01)) ,part_ptr->dz(Mpart_ptr) / (m_distance*(1 + m_aveSmR * 0.01)) };
	}

	inline part_prec getDist() {
		return m_distance;
	}


public:

	typedef part_prec(Particle::*dsmthFunc)(Particle*);
	typedef part_prec(Particle::*dsmthFuncvec3pos)(glm::vec3);

	ParticlePair(Particle* mpart, Particle* nbpart, DIMENSIONS dim)
		: Mpart_ptr(mpart), NBpart_ptr(nbpart), mode(new QuinticSpline()){
		//std::cout << "PariclePair constructor\n";
		if(!((mpart->m_type == BOUNDARY) or (nbpart->m_type == BOUNDARY))){
			mpart->addNeighbour();
			nbpart->addNeighbour();
		}
		m_dim = dim;
		saveAveSmR();
		saveDistance();
		saveDimensionlessR();
		/*
		float smr = 0.025f;
		float r = 0.f;
		float dr = 0.01f * smr;
		float factor = (15.f / (7.f * M_PI*pow(smr, 2)));
		while (r < 2.f * smr) {
			float W;
			float dW;
			float d2W;
			if ((r / smr >= 0) and (r / smr <= 1)){
				W = factor * ((2.0 / 3.0 - pow(r / smr, 2) + pow(r / smr, 3) / 2.0));
				dW = factor * (-2.0*r / smr / smr + 3.0 / 2.0*pow(r / smr, 2) / smr);
				d2W = factor* (-2.0 / pow(smr, 2) + 3.0 * r / smr / pow(smr, 2));
			}
			if ((r / smr > 1) and (r / smr <= 2)){
				W = factor * (1.0 / 6.0 * pow((2.0 - r / smr), 3));
				dW = -factor / 2.0 * pow((2.0 - r / smr), 2) / smr;
				d2W = factor * (2.f - r / smr) / pow(smr, 2);
			}
			if (r / smr > 2){
				W = 0.f;
				dW = 0.f;
				d2W = 0;
			}
			std::cout << r << "  " << W << "  " << dW << "  " << d2W << "\n";
			r += dr;
		}
		*/
	}


	float getPairDistance() { return m_distance; }


	~ParticlePair() {
		//std::cout << "PariclePair destructor\n";
		//std::cout << "~ParticlePair()\n";
		//std::cout << "  "<< Mpart_ptr->m_id << " " << NBpart_ptr->m_id << "\n";
		//delete Mpart_ptr;
		//delete NBpart_ptr;
		delete mode;
	}


	part_prec dvfunc(dsmthFunc pfcn) {
		part_prec retVal = (Mpart_ptr->*pfcn)(NBpart_ptr);
		return retVal;
	}

};
