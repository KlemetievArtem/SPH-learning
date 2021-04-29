#pragma once

class uniformTreeGrid_2D {
	uniformTreeGrid_2D* leaves[4];
	uniformTreeGrid_2D* const parent;
	cd_prec const m_xmin, m_xmax, m_ymin, m_ymax;
	std::vector<Particle*> particles_inside;
	std::set<uniformTreeGrid_2D*> neigbours;
public:
	int CheckLevel() {
		int retVal = 0;
		uniformTreeGrid_2D* current_parent_ptr = parent;
		while (current_parent_ptr != nullptr) {
			retVal++;
			current_parent_ptr = current_parent_ptr->parent;
			//std::cout << current_parent_ptr << " ";
		}
		return retVal;
		//std::cout << "\n";
	}
	uniformTreeGrid_2D(cd_prec xmin, cd_prec xmax, cd_prec ymin, cd_prec ymax, uniformTreeGrid_2D* parent_ptr = nullptr) :
		m_xmin(xmin), m_xmax(xmax), m_ymin(ymin), m_ymax(ymax),parent(parent_ptr) {
		//std::cout << "(" << xmin << "," << xmax << "," << ymin << "," << ymax << "," << id << ") ";
		//if(parent_ptr != nullptr){
		//	std::cout << m_id << "  " << parent_ptr->m_id << "\n";
		//}
		//else {
		//	std::cout << m_id << "  " << "_" << "\n";
		//}
		int EXTRNAL_PARAMETER_FOR_METHOD = 3;
		if (CheckLevel() < EXTRNAL_PARAMETER_FOR_METHOD) {
			leaves[0] = new uniformTreeGrid_2D(m_xmin, (m_xmin+m_xmax) / 2.0, m_ymin, (m_ymin+m_ymax) / 2.0,  this);
			leaves[1] = new uniformTreeGrid_2D((m_xmin + m_xmax) / 2.0, m_xmax, m_ymin, (m_ymin + m_ymax) / 2.0,  this);
			leaves[2] = new uniformTreeGrid_2D(m_xmin, (m_xmin + m_xmax) / 2.0, (m_ymin + m_ymax) / 2.0, m_ymax,  this);
			leaves[3] = new uniformTreeGrid_2D((m_xmin + m_xmax) / 2.0, m_xmax, (m_ymin + m_ymax) / 2.0, m_ymax,  this);
		}

	}

	void FindNeighboursCells(std::set<uniformTreeGrid_2D*>& allpossibleNB) {
		if (leaves[0] == nullptr) {
			allpossibleNB.insert(this);
		}
		else {
			for (auto l : leaves) {
				l->FindNeighboursCells(allpossibleNB);
			}
		}
	}

	bool CellsAreNeighbours(uniformTreeGrid_2D* p_nb, uniformTreeGrid_2D* pp_nb) {
		int numberOfSimmilarX = 0;
		int numberOfSimmilarY = 0;
		cd_prec pXvals[2] = { p_nb->m_xmin,p_nb->m_xmax };
		cd_prec pYvals[2] = { p_nb->m_ymin,p_nb->m_ymax };
		cd_prec ppXvals[2] = { pp_nb->m_xmin,pp_nb->m_xmax };
		cd_prec ppYvals[2] = { pp_nb->m_ymin,pp_nb->m_ymax };
		for (auto x : pXvals) {
			for (auto xx : ppXvals) {
				if (x == xx) numberOfSimmilarX++;
			}
		}
		for (auto y : pYvals) {
			for (auto yy : ppYvals) {
				if (y == yy) numberOfSimmilarY++;
			}
		}
		if ((numberOfSimmilarX + numberOfSimmilarY >= 3)or((numberOfSimmilarX==1)and(numberOfSimmilarY==1))) {
			return true;
		}


		return false;
	}
	
	void AssignNeighbourCells(std::set<uniformTreeGrid_2D*>& allpossibleNB) {
		for (auto p_nb : allpossibleNB) {
			for (auto pp_nb : allpossibleNB) {
				if (CellsAreNeighbours(p_nb,pp_nb)) {
					p_nb->neigbours.insert(pp_nb);
				}
			}
		}
	}
	
	void addPoint(Particle* part_ptr) {
		if (leaves[0] == nullptr) {
			//std::cout << "particle was added to " << this->m_id << " " << this->parent->m_id << "\n";
			this->particles_inside.push_back(part_ptr); 
		}
		else {
			for (auto l : leaves) {
				if (((part_ptr->m_position.val.y < l->m_ymax) and (part_ptr->m_position.val.y > l->m_ymin)) and ((part_ptr->m_position.val.x < l->m_xmax) and (part_ptr->m_position.val.x > l->m_xmin))) {
					l->addPoint(part_ptr);
				}
			}
		}
	}
	void CreatingParticlePairs(std::vector<ParticlePair*>& PartPairs, part_prec KernelCoef, DIMENSIONS dim) {
		if (leaves[0] == nullptr) {
			for (auto p : particles_inside) {
				std::set<Particle*> setOfNb;
				for (auto n : this->neigbours) {
					for (auto pp : n->particles_inside) {
						if (pp->m_id > p->m_id) {
							if (((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY))) continue; // Не учитываем пары Граница-Граница 
							if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
						}
					}
				}
				for (auto nb_part : setOfNb) { PartPairs.push_back(new ParticlePair(p, nb_part, dim)); } //CREATING PAIRS
			}
		}
		else {
			for (auto l : leaves) {
				l->CreatingParticlePairs(PartPairs, KernelCoef, dim);
			}
		}
	}

	void addingVirtualParticles(std::vector<ParticlePair*>& PartPairs, part_prec KernelCoef, DIMENSIONS dim) {
		if (leaves[0] == nullptr) {
			for (auto p : particles_inside) {
				std::set<Particle*> setOfNb;
				for (auto n : this->neigbours) {
					for (auto pp : n->particles_inside) {
						if (pp->m_id > p->m_id) {
							if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue;
							if ((pp->m_type == PARTICLETYPE::VIRTUAL) or (p->m_type == PARTICLETYPE::VIRTUAL)) {
								if ((pp->m_type == PARTICLETYPE::VIRTUAL) and (p->m_type == PARTICLETYPE::VIRTUAL)) continue;
								if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
							}
						}
					}
				}
				for (auto nb_part : setOfNb) { PartPairs.push_back(new ParticlePair(p, nb_part, dim)); } //CREATING PAIRS
			}
		}
		else {
			for (auto l : leaves) {
				//std::cout << "                    |              \n";
				l->addingVirtualParticles(PartPairs, KernelCoef, dim);
			}
		}
	}
	void clear() {
		if (leaves[0] == nullptr) {
			particles_inside.clear();
		}
		else {
			for (auto l : leaves) {
				l->clear();
			}
		}
	}

	~uniformTreeGrid_2D() {
		//leaves[0] = nullptr;
		//leaves[1] = nullptr;
		//leaves[2] = nullptr;
		//leaves[3] = nullptr;
		//clear();
		delete leaves[0];
		delete leaves[1];
		delete leaves[2];
		delete leaves[3];
	}
};




class uniformTreeGrid_3D {
	uniformTreeGrid_3D* leaves[8];
	uniformTreeGrid_3D* const parent;
	cd_prec const m_xmin, m_xmax, m_ymin, m_ymax, m_zmin, m_zmax;
	std::vector<Particle*> particles_inside;

	std::set<uniformTreeGrid_3D*> neigbours;
public:
	int CheckLevel() {
		int retVal = 0;
		uniformTreeGrid_3D* current_parent_ptr = parent;
		while (current_parent_ptr != nullptr) {
			retVal++;
			current_parent_ptr = current_parent_ptr->parent;
			//std::cout << current_parent_ptr << " ";
		}
		return retVal;
		//std::cout << "\n";
	}
	uniformTreeGrid_3D(cd_prec xmin, cd_prec xmax, cd_prec ymin, cd_prec ymax, cd_prec zmin, cd_prec zmax, uniformTreeGrid_3D* parent_ptr = nullptr) :
		m_xmin(xmin), m_xmax(xmax), m_ymin(ymin), m_ymax(ymax), m_zmin(zmin), m_zmax(zmax),  parent(parent_ptr) {
		//std::cout << "(" << xmin << "," << xmax << "," << ymin << "," << ymax << "," << id << ") ";
		//if (parent_ptr != nullptr) {
		//	std::cout << m_id << "  " << parent_ptr->m_id << "\n";
		//}
		//else {
		//	std::cout << m_id << "  " << "_" << "\n";
		//}

		int EXTRNAL_PARAMETER_FOR_METHOD = 3;
		if (CheckLevel() < EXTRNAL_PARAMETER_FOR_METHOD) {
			leaves[0] = new uniformTreeGrid_3D(m_xmin, (m_xmin + m_xmax) / 2.0, m_ymin, (m_ymin + m_ymax) / 2.0, m_zmin, (m_zmin + m_zmax) / 2.0, this);
			leaves[1] = new uniformTreeGrid_3D((m_xmin + m_xmax) / 2.0, m_xmax, m_ymin, (m_ymin + m_ymax) / 2.0, m_zmin, (m_zmin + m_zmax) / 2.0, this);
			leaves[2] = new uniformTreeGrid_3D(m_xmin, (m_xmin + m_xmax) / 2.0, (m_ymin + m_ymax) / 2.0, m_ymax, m_zmin, (m_zmin + m_zmax) / 2.0, this);
			leaves[3] = new uniformTreeGrid_3D((m_xmin + m_xmax) / 2.0, m_xmax, (m_ymin + m_ymax) / 2.0, m_ymax, m_zmin, (m_zmin + m_zmax) / 2.0, this);
			leaves[4] = new uniformTreeGrid_3D(m_xmin, (m_xmin + m_xmax) / 2.0, m_ymin, (m_ymin + m_ymax) / 2.0, (m_zmin + m_zmax) / 2.0, m_zmax, this);
			leaves[5] = new uniformTreeGrid_3D((m_xmin + m_xmax) / 2.0, m_xmax, m_ymin, (m_ymin + m_ymax) / 2.0, (m_zmin + m_zmax) / 2.0, m_zmax, this);
			leaves[6] = new uniformTreeGrid_3D(m_xmin, (m_xmin + m_xmax) / 2.0, (m_ymin + m_ymax) / 2.0, m_ymax, (m_zmin + m_zmax) / 2.0, m_zmax, this);
			leaves[7] = new uniformTreeGrid_3D((m_xmin + m_xmax) / 2.0, m_xmax, (m_ymin + m_ymax) / 2.0, m_ymax, (m_zmin + m_zmax) / 2.0, m_zmax, this);
		}

	}

	void FindNeighboursCells(std::set<uniformTreeGrid_3D*>& allpossibleNB) {
		if (leaves[0] == nullptr) {
			allpossibleNB.insert(this);
		}
		else {
			for (auto l : leaves) {
				l->FindNeighboursCells(allpossibleNB);
			}
		}
	}

	bool CellsAreNeighbours(uniformTreeGrid_3D* p_nb, uniformTreeGrid_3D* pp_nb) {
		int numberOfSimmilarX = 0;
		int numberOfSimmilarY = 0;
		int numberOfSimmilarZ = 0;
		cd_prec pXvals[2] = { p_nb->m_xmin,p_nb->m_xmax };
		cd_prec pYvals[2] = { p_nb->m_ymin,p_nb->m_ymax };
		cd_prec pZvals[2] = { p_nb->m_zmin,p_nb->m_zmax };
		cd_prec ppXvals[2] = { pp_nb->m_xmin,pp_nb->m_xmax };
		cd_prec ppYvals[2] = { pp_nb->m_ymin,pp_nb->m_ymax };
		cd_prec ppZvals[2] = { pp_nb->m_zmin,pp_nb->m_zmax };
		for (auto x : pXvals) {
			for (auto xx : ppXvals) {
				if (x == xx) numberOfSimmilarX++;
			}
		}
		for (auto y : pYvals) {
			for (auto yy : ppYvals) {
				if (y == yy) numberOfSimmilarY++;
			}
		}
		for (auto z : pZvals) {
			for (auto zz : ppZvals) {
				if (z == zz) numberOfSimmilarZ++;
			}
		}
		if (((numberOfSimmilarX==2)and(numberOfSimmilarY == 2))or((numberOfSimmilarX == 2) and (numberOfSimmilarZ == 2))or((numberOfSimmilarY == 2) and (numberOfSimmilarZ == 2))
			or((numberOfSimmilarX == 1) and (numberOfSimmilarY == 1) and (numberOfSimmilarZ == 1)) or ((numberOfSimmilarX == 2) and (numberOfSimmilarY == 1) and (numberOfSimmilarZ == 1))
			or ((numberOfSimmilarX == 1) and (numberOfSimmilarY == 2) and (numberOfSimmilarZ == 1)) or ((numberOfSimmilarX == 1) and (numberOfSimmilarY == 1) and (numberOfSimmilarZ == 2))) {
			return true;
		}


		return false;
	}

	void AssignNeighbourCells(std::set<uniformTreeGrid_3D*>& allpossibleNB) {
		for (auto p_nb : allpossibleNB) {
			for (auto pp_nb : allpossibleNB) {
				if (CellsAreNeighbours(p_nb, pp_nb)) {
					p_nb->neigbours.insert(pp_nb);
				}
			}
		}

	}

	void addPoint(Particle* part_ptr) {
		if (leaves[0] == nullptr) {
			//std::cout << "particle was added to " << this->m_id << " " << this->parent->m_id << "\n";
			this->particles_inside.push_back(part_ptr);
		}
		else {
			for (auto l : leaves) {
				if (((part_ptr->m_position.val.y < l->m_ymax) and (part_ptr->m_position.val.y > l->m_ymin)) and ((part_ptr->m_position.val.x < l->m_xmax) and (part_ptr->m_position.val.x > l->m_xmin))) {
					l->addPoint(part_ptr);
				}
			}
		}
	}
	void CreatingParticlePairs(std::vector<ParticlePair*>& PartPairs, part_prec KernelCoef, DIMENSIONS dim) {
		if (leaves[0] == nullptr) {
			for (auto p : particles_inside) {
				std::set<Particle*> setOfNb;
				for (auto n : this->neigbours) {
					for (auto pp : n->particles_inside) {
						if (pp->m_id > p->m_id) {
							if (((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY))) continue; // Не учитываем пары Граница-Граница 
							if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
						}
					}
				}
				for (auto nb_part : setOfNb) { PartPairs.push_back(new ParticlePair(p, nb_part, dim)); } //CREATING PAIRS
			}
		}
		else {
			for (auto l : leaves) {
				l->CreatingParticlePairs(PartPairs, KernelCoef, dim);
			}
		}
	}

	void addingVirtualParticles(std::vector<ParticlePair*>& PartPairs, part_prec KernelCoef, DIMENSIONS dim) {
		if (leaves[0] == nullptr) {
			for (auto p : particles_inside) {
				std::set<Particle*> setOfNb;
				for (auto n : this->neigbours) {
					for (auto pp : n->particles_inside) {
						if (pp->m_id > p->m_id) {
							if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue;
							if ((pp->m_type == PARTICLETYPE::VIRTUAL) or (p->m_type == PARTICLETYPE::VIRTUAL)) {
								if ((pp->m_type == PARTICLETYPE::VIRTUAL) and (p->m_type == PARTICLETYPE::VIRTUAL)) continue;
								if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
							}
						}
					}
				}
				for (auto nb_part : setOfNb) { PartPairs.push_back(new ParticlePair(p, nb_part, dim)); } //CREATING PAIRS
			}
		}
		else {
			for (auto l : leaves) {
				//std::cout << "                    |              \n";
				l->addingVirtualParticles(PartPairs, KernelCoef, dim);
			}
		}
	}


	void clear() {
		if (leaves[0] == nullptr) {
			particles_inside.clear();
		}
		else {
			for (auto l : leaves) {
				l->clear();
			}
		}
	}

	~uniformTreeGrid_3D() {
		//leaves[0] = nullptr;
		//leaves[1] = nullptr;
		//leaves[2] = nullptr;
		//leaves[3] = nullptr;
		//clear();
		delete leaves[0];
		delete leaves[1];
		delete leaves[2];
		delete leaves[3];
		delete leaves[4];
		delete leaves[5];
		delete leaves[6];
		delete leaves[7];
	}
};
