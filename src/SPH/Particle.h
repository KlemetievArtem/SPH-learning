#pragma once


typedef double part_prec;
typedef glm::highp_dvec3 part_prec_3;



enum PARTICLETYPE
{
	REAL,
	BOUNDARY,
	VIRTUAL,
	ALLTYPES
};

template <class T>
struct Ch_ {
public:
	Ch_(T v=static_cast<T>(0)) : val(v) {}
	T val;
	T dval;
};



struct Particle {
public:
	int m_id;
	PARTICLETYPE m_type;
	Ch_ <part_prec_3> m_position;  //part_prec_3
	Ch_ <part_prec_3> m_velocity;  //part_prec_3

	float absVelocity;           //part_prec
	part_prec m_SmR;                 //part_prec
	int nrOfNeighbours = 0;


	Particle(int id, PARTICLETYPE type, glm::vec3 pos, glm::vec3 vel, float SmR, part_prec density, float mass)
		: m_id(id), m_type(type), m_SmR(SmR), m_mass(mass) {
		m_position.val = pos; m_position.dval = part_prec_3(0.0);
		m_velocity.val = vel; m_velocity.dval = part_prec_3(0.0);
		m_density.val = density; m_density.dval = part_prec(0.0);
		//if (m_type == PARTICLETYPE::REAL)
		p_art_water();
		//else
			//p_gas();
		setDVisc();
	}

	Particle(Particle* part): Particle(part->m_id, part->m_type, part->m_position.val, part->m_velocity.val, part->m_SmR, part->m_density.val, part->m_mass){
		this->m_position.dval = part->m_position.dval;
		this->m_velocity.dval = part->m_velocity.dval;
		this->m_density.dval = part->m_density.dval;
		this->InitPolygonNormal = part->InitPolygonNormal;
		this->m_pressure.val = part->m_pressure.val;
		this->m_pressure.dval = part->m_pressure.dval;
		this->m_DVisc = part->m_DVisc;
		//this->VirtualCounterpart = part->VirtualCounterpart;
		//
		//
		//for (auto vcn : this->VirtualCounterpartNormals) {
		//	part->VirtualCounterpartNormals.push_back(vcn);
		//}
		//for (auto vcf : this->VirtualCounterpartFlags) {
		//	part->VirtualCounterpartFlags.push_back(vcf);
		//}
		//for (auto vcd : this->VC_DistanceToBoundary) {
		//	part->VC_DistanceToBoundary.push_back(vcd);
		//}
		//for (auto vcv : this->VC_BoundaryVelocity) {
		//	part->VC_BoundaryVelocity.push_back(vcv);
		//}

	}

	~Particle() {
		//VirtualCounterpartNormals.resize(0);
		//VirtualCounterpartFlags.resize(0);
	}

	glm::vec3 m_color{ 0.9f, 0.9f, 0.9f };

	//float correctiveKernel;

	glm::vec3 oldAcceleration;
	cd_prec m_mass;                //part_prec
	Ch_<cd_prec> m_density;        //part_prec
	Ch_<float> m_pressure;       //part_prec

	float m_SoundVelocity;
	float m_Temperature = 357.f;
	float m_DVisc;
	//std::vector<glm::vec3> m_eps;







	//part_prec distance_before_teleportation=0.0;
	//part_prec_3 normal_of_periodic_boundary;
	//part_prec_3 position_of_periodic_boundary;
	//ƒÀﬂ —Œ«ƒ¿Õ»ﬂ ¬»–“”¿À‹Õ€’ ◊¿—“»÷
	//bool VirtualCounterpart = false;
	//std::vector<part_prec_3> VirtualCounterpartNormals;
	//std::vector<bool> VirtualCounterpartFlags;
	//std::vector<part_prec> VC_DistanceToBoundary;
	//std::vector<bool> periodicBoundary;
	//std::vector<part_prec_3> VC_DistanceToPeriodic;
	//std::vector<part_prec_3> VC_Before_PeriodicVelocity;
	//std::vector<part_prec_3> VC_Before_PeriodicPosition;
	//std::vector<part_prec_3> VC_BoundaryVelocity;

	—D_Boundary* particle_boundary;
	void assignToBoundary(—D_Boundary* cd_boundaty) { particle_boundary = cd_boundaty; }

	float dx(Particle* other) { return this->m_position.val.x - other->m_position.val.x; }
	float dy(Particle* other) { return this->m_position.val.y - other->m_position.val.y; }
	float dz(Particle* other) {	return this->m_position.val.z - other->m_position.val.z; }
	float distance(Particle* other) { return std::sqrt(pow(this->dx(other), 2) + pow(this->dy(other), 2) + pow(this->dz(other), 2)); }



	float dVx(Particle* other) { return this->m_velocity.val.x - other->m_velocity.val.x; }
	float dVy(Particle* other) { return this->m_velocity.val.y - other->m_velocity.val.y; }
	float dVz(Particle* other) { return this->m_velocity.val.z - other->m_velocity.val.z; }

	part_prec_3 dV(Particle* other) { return { this->dVx(other),this->dVy(other),this->dVz(other) }; }


	float dx(glm::vec3 pos) { return this->m_position.val.x - pos.x; }
	float dy(glm::vec3 pos) { return this->m_position.val.y - pos.y; }
	float dz(glm::vec3 pos) { return this->m_position.val.z - pos.z; }
	float distance(glm::vec3 pos) { return std::sqrt(pow(this->dx(pos), 2) + pow(this->dy(pos), 2) + pow(this->dz(pos), 2)); }
	

	void p_art_water() {
		float gamma = 7.0f;
		float Dens0 = 1000;
		float b = 1.013E05f;
		setSoundVel(1480.f);
		//p = b * ((m_dens / Dens0)**gamma - 1)
		//m_pressure.val = (b*pow((m_density.val / Dens0), gamma) - 1);
		m_pressure.val = Dens0 * pow(m_SoundVelocity,2) / gamma * (pow(m_density.val / Dens0, gamma) - 1.f);
		//std::cout << "pressure = " << Dens0 << "*" << pow(m_SoundVelocity, 2) << "/" << gamma << "*(" << "(" << m_density.val / Dens0 << ")^" << gamma << "-" << 1.f << ")" << "=" << m_pressure.val << "\n";
	}
	void p_gas() {
		float gamma = 1.4f;
		m_pressure.val = ((gamma - 1)*m_density.val*m_Temperature);
		setSoundVel(sqrt((gamma - 1)*m_Temperature));

	}

	void addNeighbour() {
		nrOfNeighbours++;
	}
	void refreshNeighbours() {
		nrOfNeighbours = 0;
	}

private:
	void setSoundVel(float val) {
		m_SoundVelocity = val;
	}


	void setDVisc() {
		if (m_type == PARTICLETYPE::REAL)
			m_DVisc = 0;
		if (m_type == PARTICLETYPE::VIRTUAL)
			m_DVisc = 0;
		if (m_type == PARTICLETYPE::BOUNDARY)
			m_DVisc = 1.0E-03;


		m_DVisc = 1000E-06;
	}






public: 
	glm::vec3 InitPolygonNormal;// = glm::vec3(0.f);
	glm::vec3 GetNormal() { return this->InitPolygonNormal; }
};



