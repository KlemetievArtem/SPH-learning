#pragma once



class BoundaryData {
private:
	bool VirtualCounterpart = false;
	std::vector<part_prec_3> VirtualCounterpartNormals;
	std::vector<bool> VirtualCounterpartFlags;
	std::vector<part_prec> VC_DistanceToBoundary;
	std::vector<bool> periodicBoundary;

	std::vector<part_prec_3> VC_DistanceToPeriodic;
	std::vector<part_prec_3> VC_Before_PeriodicVelocity;
	std::vector<part_prec_3> VC_Before_PeriodicPosition;
	std::vector<part_prec_3> VC_BoundaryVelocity;
public:
	BoundaryData(){ 
		assert("BoundaryData" && 0);
	}
	BoundaryData(part_prec_3 normal, part_prec distance, part_prec_3 Velocity) {
		VirtualCounterpartNormals.push_back(normal);
		VirtualCounterpartFlags.push_back(true);
		VC_DistanceToBoundary.push_back(distance);
		//VC_BoundaryVelocity.push_back(SPH.Particles[neighbour_p->m_id]->m_velocity.val);
		VC_BoundaryVelocity.push_back(Velocity);
	}
	void furtherInitialization(bool condition, part_prec_3 distanceToPeriodic_part, part_prec_3 virtualPosition, part_prec_3 virtualVelocity) {
		if (condition) {
			periodicBoundary.push_back(condition);
			VC_DistanceToPeriodic.push_back(distanceToPeriodic_part);
			VC_Before_PeriodicPosition.push_back(virtualPosition);
			VC_Before_PeriodicVelocity.push_back(virtualVelocity);
		}
		else {
			periodicBoundary.push_back(condition);
		}
	}

	void addingBoundary(part_prec_3 normal, part_prec distance, part_prec_3 Velocity, bool condition, part_prec_3 distanceToPeriodic_part, part_prec_3 virtualPosition, part_prec_3 virtualVelocity) {
		VirtualCounterpartNormals.push_back(normal);
		VirtualCounterpartFlags.push_back(true);
		VC_DistanceToBoundary.push_back(distance);
		VC_BoundaryVelocity.push_back(Velocity);
		if (condition) {
			periodicBoundary.push_back(condition);
			VC_DistanceToPeriodic.push_back(distanceToPeriodic_part);
			VC_Before_PeriodicPosition.push_back(virtualPosition);
			VC_Before_PeriodicVelocity.push_back(virtualVelocity);
		}
		else {
			periodicBoundary.push_back(condition);
		}
	
	}




	bool firstBoundaryisPeriodic() { return periodicBoundary[0]; }

	part_prec_3 firstBoundaryVelocity() { return VC_BoundaryVelocity[0]; }
	part_prec_3 firstDistanceVector() { return VirtualCounterpartNormals[0] * VC_DistanceToBoundary[0]; }
	part_prec_3 firstNormal() { return VirtualCounterpartNormals[0]; }
	part_prec_3 firstDistanceToPeriodic() { return VC_DistanceToPeriodic[0]; }
	part_prec_3 firstPeriodicBoundaryVelocity() { return VC_Before_PeriodicVelocity[0]; }

	void hasVirtualCounterpart() { VirtualCounterpart = true; }
	void VirtualCounterpartReset() { VirtualCounterpart = false; }
	bool virtualCounterpartCheck() { return VirtualCounterpart; }
	int NrOfNormals() { return VirtualCounterpartNormals.size(); }
	void resetAll() {
		VirtualCounterpartNormals.resize(0);
		VirtualCounterpartFlags.resize(0);
		VC_DistanceToBoundary.resize(0);
		VC_BoundaryVelocity.resize(0);
		VirtualCounterpart = false;

		periodicBoundary.resize(0);
		VC_DistanceToPeriodic.resize(0);
		VC_Before_PeriodicPosition.resize(0);
		VC_Before_PeriodicVelocity.resize(0);
	}
	bool newVirtualCounterpart(part_prec_3& normal, float angle_cos) {
		bool differentNormal;
		for (int i = 0;i < VirtualCounterpartNormals.size();i++) {
			if (abs(VirtualCounterpartNormals[i].x*normal.x +VirtualCounterpartNormals[i].y*normal.y + VirtualCounterpartNormals[i].z*normal.z) / (sqrt(pow(VirtualCounterpartNormals[i].x, 2) + pow(VirtualCounterpartNormals[i].y, 2) + pow(VirtualCounterpartNormals[i].z, 2))*sqrt(pow(normal.x, 2) + pow(normal.y, 2) + pow(normal.z, 2))) >= angle_cos) {
				differentNormal = true;
			}
			else {
				differentNormal = false;
			}
		}
		return differentNormal;


		//for (int i = 0;i < SPH.Particles[real_p->m_id]->VirtualCounterpartNormals.size();i++) {
		//	if (abs(SPH.Particles[real_p->m_id]->VirtualCounterpartNormals[i].x*normal.x + SPH.Particles[real_p->m_id]->VirtualCounterpartNormals[i].y*normal.y + SPH.Particles[real_p->m_id]->VirtualCounterpartNormals[i].z*normal.z) / (sqrt(pow(SPH.Particles[real_p->m_id]->VirtualCounterpartNormals[i].x, 2) + pow(SPH.Particles[real_p->m_id]->VirtualCounterpartNormals[i].y, 2) + pow(SPH.Particles[real_p->m_id]->VirtualCounterpartNormals[i].z, 2))*sqrt(pow(normal.x, 2) + pow(normal.y, 2) + pow(normal.z, 2))) >= 0.9f) {
		//		SPH.Particles[real_p->m_id]->VirtualCounterpart = SPH.Particles[real_p->m_id]->VirtualCounterpartFlags[i];
		//	}	
		//}

	}

};

class TeleportationData {
private:
	part_prec distance_before_teleportation;
	part_prec_3 normal_of_periodic_boundary;
	part_prec_3 position_of_periodic_boundary;
public:
	TeleportationData() { assert("TeleportationData" && 0); }
	TeleportationData(part_prec distance, part_prec_3 normal, part_prec_3 distanceToPeriodic_part) {
		distance_before_teleportation = distance;
		normal_of_periodic_boundary = normal;
		position_of_periodic_boundary = distanceToPeriodic_part;
	}
	part_prec distanceToBoundary() { return distance_before_teleportation; }
	part_prec_3 boundaryNormal(){ return normal_of_periodic_boundary; }
	part_prec_3 teleportingDistance() { return position_of_periodic_boundary; }
};





class BoundaryMentor {
public:
	BoundaryMentor() {}
	std::map<int, BoundaryData> BoundaryLinks;
	std::map<int, TeleportationData> TeleportationLinks;


	void resetBoundaryData() {
		for (auto& link : BoundaryLinks) {
			link.second.resetAll();
		}
		BoundaryLinks.clear();
	}

	void resetTeleportationData() {
		TeleportationLinks.clear();
	}



	bool boundarySecBar1() {
		if (BoundaryLinks.size() != 0)
			return true;
		return false;
	}
	bool boundarySecBar2(int part_id) {
		if (BoundaryLinks.count(part_id) != 0){
			return true;
		}
		return false;
	}




};