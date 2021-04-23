#pragma once



enum UGIdentifier {
	ID_0_0,
	ID_1_0,
	ID_0_1,
	ID_1_1,
	PARENT
};

class uniformTreeGrid_2D {
	uniformTreeGrid_2D* leaves[4];
	uniformTreeGrid_2D* parent;
	UGIdentifier m_id;
	cd_prec m_xmin, m_xmax, m_ymin, m_ymax;
	std::vector<Particle*> particles_inside;
public:
	uniformTreeGrid_2D* getParentptr() {
		return parent;
	}
	int CheckLevel() {
		int retVal = 0;
		uniformTreeGrid_2D* current_parent_ptr = parent;
		while (current_parent_ptr != nullptr) {
			retVal++;
			current_parent_ptr = current_parent_ptr->getParentptr();
			//std::cout << current_parent_ptr << " ";
		}
		return retVal;
		//std::cout << "\n";
	}
	uniformTreeGrid_2D(cd_prec xmin, cd_prec xmax, cd_prec ymin, cd_prec ymax, UGIdentifier id = PARENT, uniformTreeGrid_2D* parent_ptr = nullptr) {
		std::cout << "(" << xmin << "," << xmax << "," << ymin << "," << ymax << "," << id << ") ";
		m_xmin = xmin;
		m_xmax = xmax;
		m_ymin = ymin;
		m_ymax = ymax;
		m_id = id;
		parent = parent_ptr;
		std::cout << this << " " << parent << "\n";
		if (CheckLevel() < 3) {
			leaves[0] = new uniformTreeGrid_2D(m_xmin, (m_xmin+m_xmax) / 2.0, m_ymin, (m_ymin+m_ymax) / 2.0, ID_0_0, this);
			leaves[1] = new uniformTreeGrid_2D((m_xmin + m_xmax) / 2.0, m_xmax, m_ymin, (m_ymin + m_ymax) / 2.0, ID_1_0, this);
			leaves[2] = new uniformTreeGrid_2D(m_xmin, (m_xmin + m_xmax) / 2.0, (m_ymin + m_ymax) / 2.0, m_ymax, ID_0_1, this);
			leaves[3] = new uniformTreeGrid_2D((m_xmin + m_xmax) / 2.0, m_xmax, (m_ymin + m_ymax) / 2.0, m_ymax, ID_1_1, this);
		}

	}
	void addPoint(Particle* part_ptr) {
		if (leaves[0] == nullptr) {
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
  		//std::cout << "(" << m_xmin << "," << m_xmax << "," << m_ymin << "," << m_ymax << "," << m_id << ")\n";
		if (leaves[0] == nullptr) {
			//std::cout << "(" << m_xmin << "," << m_xmax << "," << m_ymin << "," << m_ymax << "," << m_id << ")\n";
			//std::cout << this->particles_inside.size() << "\n";
			if (this->parent->parent == nullptr) {
				switch (this->m_id)	{
				case(ID_0_0):
					for (auto p : particles_inside) {
						std::set<Particle*> setOfNb;
						for (auto pp : particles_inside) {
							if (pp->m_id > p->m_id) {
								if (((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY))) continue; // Не учитываем пары Граница-Граница  
								//if (((pp->m_type == PARTICLETYPE::VIRTUAL) and (p->m_type == PARTICLETYPE::VIRTUAL))) continue; // Не учитываем пары Граница-Граница  
								//if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) PartPairs.push_back(new ParticlePair(p, pp, dim)); //CREATING PAIRS
								if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
							}
						}
						for (auto pp : this->parent->leaves[1]->particles_inside) {
							if (pp->m_id > p->m_id) {
								if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница  
								//if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) PartPairs.push_back(new ParticlePair(p, pp, dim)); //CREATING PAIRS
								if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
							}
						}
						for (auto pp : this->parent->leaves[2]->particles_inside) {
							if (pp->m_id > p->m_id) {
								if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница  
								//if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) PartPairs.push_back(new ParticlePair(p, pp, dim)); //CREATING PAIRS
								if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
							}
						}
						for (auto pp : this->parent->leaves[3]->particles_inside) {
							if (pp->m_id > p->m_id) {
								if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница  
								//if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) PartPairs.push_back(new ParticlePair(p, pp, dim)); //CREATING PAIRS
								if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
							}
						}
						for (auto nb_part : setOfNb) {
							if (nb_part->m_id > p->m_id) {
								//std::cout << p->m_id << "    " << nb_part->m_id << "\n";
								PartPairs.push_back(new ParticlePair(p, nb_part, dim)); //CREATING PAIRS
							}
						}
					}
					for (auto p : this->parent->leaves[1]->particles_inside) {
						std::set<Particle*> setOfNb;
						for (auto pp : this->parent->leaves[1]->particles_inside) {
							if (pp->m_id > p->m_id) {
								if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница  
								//if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) PartPairs.push_back(new ParticlePair(p, pp, dim)); //CREATING PAIRS
								if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
							}
						}
						for (auto nb_part : setOfNb) {
							if (nb_part->m_id > p->m_id) {
								//std::cout << p->m_id << "    " << nb_part->m_id << "\n";
								PartPairs.push_back(new ParticlePair(p, nb_part, dim)); //CREATING PAIRS
							}
						}
					}
					for (auto p : this->parent->leaves[2]->particles_inside) {
						std::set<Particle*> setOfNb;
						for (auto pp : this->parent->leaves[2]->particles_inside) {
							if (pp->m_id > p->m_id) {
								if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница  
								//if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) PartPairs.push_back(new ParticlePair(p, pp, dim)); //CREATING PAIRS
								if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
							}
						}
						for (auto nb_part : setOfNb) {
							if (nb_part->m_id > p->m_id) {
								//std::cout << p->m_id << "    " << nb_part->m_id << "\n";
								PartPairs.push_back(new ParticlePair(p, nb_part, dim)); //CREATING PAIRS
							}
						}
					}
					for (auto p : this->parent->leaves[3]->particles_inside) {
						std::set<Particle*> setOfNb;
						for (auto pp : this->parent->leaves[3]->particles_inside) {
							if (pp->m_id > p->m_id) {
								if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница  
								//if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) PartPairs.push_back(new ParticlePair(p, pp, dim)); //CREATING PAIRS
								if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
							}
						}
						for (auto nb_part : setOfNb) {
							if (nb_part->m_id > p->m_id) {
								//std::cout << p->m_id << "    " << nb_part->m_id << "\n";
								PartPairs.push_back(new ParticlePair(p, nb_part, dim)); //CREATING PAIRS
							}
						}
					}

					break;
				default:
					break;
				}
			}
			else {
				for (auto p : particles_inside) {
					std::set<Particle*> setOfNb;
					uniformTreeGrid_2D* First_gen_ptr = parent;
					while (First_gen_ptr->getParentptr() != nullptr) { First_gen_ptr = First_gen_ptr->getParentptr(); } // ищем верхний UG
					switch (this->m_id) {
					case(ID_0_0):
						for (auto pp : particles_inside) {
							if(pp->m_id > p->m_id){
								if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница  
								//if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) PartPairs.push_back(new ParticlePair(p, pp, dim)); //CREATING PAIRS
								if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
							}
						}
						for (auto pp : this->parent->leaves[1]->particles_inside) {
							if (pp->m_id > p->m_id) {
								if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница  
								//if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) PartPairs.push_back(new ParticlePair(p, pp, dim)); //CREATING PAIRS
								if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
							}
						}
						for (auto pp : this->parent->leaves[2]->particles_inside) {
							if (pp->m_id > p->m_id) {
								if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница  
								//if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) PartPairs.push_back(new ParticlePair(p, pp, dim)); //CREATING PAIRS
								if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
							}
						}
						for (auto pp : this->parent->leaves[3]->particles_inside) {
							if (pp->m_id > p->m_id) {
								if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница  
								//if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) PartPairs.push_back(new ParticlePair(p, pp, dim)); //CREATING PAIRS
								if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
							}
						}
						break;
					case(ID_0_1):
						for (auto pp : particles_inside) {
							if (pp->m_id > p->m_id) {
								if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница  
								//if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) PartPairs.push_back(new ParticlePair(p, pp, dim)); //CREATING PAIRS
								if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
							}
						}
						if ((First_gen_ptr->m_id = ID_0_0) or (First_gen_ptr->m_id = ID_1_0)) {
							for (auto pp : this->parent->parent->leaves[1]->leaves[0]->particles_inside) {
								if (pp->m_id > p->m_id) {
									if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница  
									//if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) PartPairs.push_back(new ParticlePair(p, pp, dim)); //CREATING PAIRS
									if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
								}
							}
							for (auto pp : this->parent->parent->leaves[1]->leaves[2]->particles_inside) {
								if (pp->m_id > p->m_id) {
									if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница  
									//if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) PartPairs.push_back(new ParticlePair(p, pp, dim)); //CREATING PAIRS
									if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
								}
							}
						}
						for (auto pp : this->parent->leaves[3]->particles_inside) {
							if (pp->m_id > p->m_id) {
								if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница  
								//if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) PartPairs.push_back(new ParticlePair(p, pp, dim)); //CREATING PAIRS
								if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
							}
						}

						break;
					case(ID_1_0):
						for (auto pp : particles_inside) {
							if (pp->m_id > p->m_id) {
								if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница  
								//if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) PartPairs.push_back(new ParticlePair(p, pp, dim)); //CREATING PAIRS
								if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
							}
						}
						if ((First_gen_ptr->m_id = ID_0_0) or (First_gen_ptr->m_id = ID_0_1)) {
							for (auto pp : this->parent->parent->leaves[2]->leaves[0]->particles_inside) {
								if (pp->m_id > p->m_id) {
									if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница  
									//if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) PartPairs.push_back(new ParticlePair(p, pp, dim)); //CREATING PAIRS
									if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
								}
							}
							for (auto pp : this->parent->parent->leaves[2]->leaves[1]->particles_inside) {
								if (pp->m_id > p->m_id) {
									if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница  
									//if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) PartPairs.push_back(new ParticlePair(p, pp, dim)); //CREATING PAIRS
									if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
								}
							}
						}
						for (auto pp : this->parent->leaves[3]->particles_inside) {
							if (pp->m_id > p->m_id) {
								if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница  
								//if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) PartPairs.push_back(new ParticlePair(p, pp, dim)); //CREATING PAIRS
								if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
							}
						}
						break;
					case(ID_1_1):
						for (auto pp : particles_inside) {
							if (pp->m_id > p->m_id) {
								if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница  
								//if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) PartPairs.push_back(new ParticlePair(p, pp, dim)); //CREATING PAIRS
								if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
							}
						}
						if ((First_gen_ptr->m_id = ID_0_0) or (First_gen_ptr->m_id = ID_1_0)) {
							for (auto pp : this->parent->parent->leaves[1]->leaves[2]->particles_inside) {
								if (pp->m_id > p->m_id) {
									if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница 
									//if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) PartPairs.push_back(new ParticlePair(p, pp, dim)); //CREATING PAIRS
									if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
								}
							}
						}
						if ((First_gen_ptr->m_id = ID_0_0) or (First_gen_ptr->m_id = ID_0_1)) {
							for (auto pp : this->parent->parent->leaves[2]->leaves[1]->particles_inside) {
								if (pp->m_id > p->m_id) {
									if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница  
									//if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) PartPairs.push_back(new ParticlePair(p, pp, dim)); //CREATING PAIRS
									if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
								}
							}
						}
						if ((First_gen_ptr->m_id = ID_0_0)){
							for (auto pp : this->parent->parent->leaves[3]->leaves[0]->particles_inside) {
								if (pp->m_id > p->m_id) {
									if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница  
									//if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) PartPairs.push_back(new ParticlePair(p, pp, dim)); //CREATING PAIRS
									if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
								}
							}
						}
						break;
					default:
						break;
					}
					for (auto nb_part : setOfNb) {
						if (nb_part->m_id > p->m_id) {
							//std::cout << p->m_id << "    " << nbpart->m_id << "\n";
							PartPairs.push_back(new ParticlePair(p, nb_part, dim)); //CREATING PAIRS
						}
					}
				}
			}

			//std::cout << "SPH_CD::neighbourSearch::PARTICLE_PAIRS_WERE_CREATED::" << PartPairs.size() << "\n";
			//switch (this->m_id)	{
			//case(ID_0_0):
			//	for (auto p : particles_inside) {
			//		for (auto pp : particles_inside) {
			//			if (p == pp) continue;
			//			if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница  
			//			if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) PartPairs.push_back(new ParticlePair(p, pp, dim)); //CREATING PAIRS
			//		}
			//		for (auto pp : this->parent->leaves[1]->particles_inside) {
			//			if (p == pp) continue;
			//			if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница  
			//			if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) PartPairs.push_back(new ParticlePair(p, pp, dim)); //CREATING PAIRS
			//		}
			//		for (auto pp : this->parent->leaves[2]->particles_inside) {
			//			if (p == pp) continue;
			//			if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница  
			//			if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) PartPairs.push_back(new ParticlePair(p, pp, dim)); //CREATING PAIRS
			//		}
			//		for (auto pp : this->parent->leaves[3]->particles_inside) {
			//			if (p == pp) continue;
			//			if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница  
			//			if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) PartPairs.push_back(new ParticlePair(p, pp, dim)); //CREATING PAIRS
			//		}
			//	}
			//	break;
			//case(ID_0_1):
			//	if (this->parent->parent == nullptr) { }
			//	else {
			//		for (auto p : particles_inside) {
			//			for (auto pp : particles_inside) {
			//				if (p == pp) continue;
			//				if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница  
			//				if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) PartPairs.push_back(new ParticlePair(p, pp, dim)); //CREATING PAIRS
			//			}
			//			for (auto pp : this->parent->parent->leaves[1]->leaves[0]->particles_inside) {
			//				if (p == pp) continue;
			//				if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница  
			//				if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) PartPairs.push_back(new ParticlePair(p, pp, dim)); //CREATING PAIRS
			//			}
			//			for (auto pp : this->parent->parent->leaves[1]->leaves[2]->particles_inside) {
			//				if (p == pp) continue;
			//				if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница  
			//				if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) PartPairs.push_back(new ParticlePair(p, pp, dim)); //CREATING PAIRS
			//			}
			//			for (auto pp : this->parent->leaves[3]->particles_inside) {
			//				if (p == pp) continue;
			//				if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница  
			//				if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) PartPairs.push_back(new ParticlePair(p, pp, dim)); //CREATING PAIRS
			//			}
			//		}
			//	}
			//	break;
			//case(ID_1_0):
			//	if (this->parent->parent == nullptr) {}
			//	else {
			//		for (auto p : particles_inside) {
			//			for (auto pp : particles_inside) {
			//				if (p == pp) continue;
			//				if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница  
			//				if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) PartPairs.push_back(new ParticlePair(p, pp, dim)); //CREATING PAIRS
			//			}
			//			for (auto pp : this->parent->parent->leaves[2]->leaves[0]->particles_inside) {
			//				if (p == pp) continue;
			//				if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница  
			//				if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) PartPairs.push_back(new ParticlePair(p, pp, dim)); //CREATING PAIRS
			//			}
			//			for (auto pp : this->parent->parent->leaves[2]->leaves[1]->particles_inside) {
			//				if (p == pp) continue;
			//				if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница  
			//				if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) PartPairs.push_back(new ParticlePair(p, pp, dim)); //CREATING PAIRS
			//			}
			//			for (auto pp : this->parent->leaves[3]->particles_inside) {
			//				if (p == pp) continue;
			//				if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница  
			//				if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) PartPairs.push_back(new ParticlePair(p, pp, dim)); //CREATING PAIRS
			//			}
			//		}
			//	}
			//	break;
			//case(ID_1_1):
			//	if (this->parent->parent == nullptr) {}
			//	else {
			//		for (auto p : particles_inside) {
			//			for (auto pp : particles_inside) {
			//				if (p == pp) continue;
			//				if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница  
			//				if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) PartPairs.push_back(new ParticlePair(p, pp, dim)); //CREATING PAIRS
			//			}
			//			for (auto pp : this->parent->parent->leaves[1]->leaves[2]->particles_inside) {
			//				if (p == pp) continue;
			//				if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница  
			//				if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) PartPairs.push_back(new ParticlePair(p, pp, dim)); //CREATING PAIRS
			//			}
			//			for (auto pp : this->parent->parent->leaves[2]->leaves[1]->particles_inside) {
			//				if (p == pp) continue;
			//				if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница  
			//				if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) PartPairs.push_back(new ParticlePair(p, pp, dim)); //CREATING PAIRS
			//			}
			//			for (auto pp : this->parent->parent->leaves[3]->leaves[0]->particles_inside) {
			//				if (p == pp) continue;
			//				if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница  
			//				if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) PartPairs.push_back(new ParticlePair(p, pp, dim)); //CREATING PAIRS
			//			}
			//		}
			//	}
			//	break;
			//default:
			//	break;
			//}

		}
		else {
			for (auto l : leaves) {
				//std::cout << "                    |              \n";
				l->CreatingParticlePairs(PartPairs, KernelCoef, dim);
			}
		}
	}

	void addingVirtualParticles(std::vector<ParticlePair*>& PartPairs, part_prec KernelCoef, DIMENSIONS dim) {
		//std::cout << "(" << m_xmin << "," << m_xmax << "," << m_ymin << "," << m_ymax << "," << m_id << ")\n";
		if (leaves[0] == nullptr) {

			if (this->parent->parent == nullptr) {
				switch (this->m_id) {
				case(ID_0_0):
					for (auto p : particles_inside) {
						std::set<Particle*> setOfNb;
						for (auto pp : particles_inside) {
							if (pp->m_id > p->m_id) {
								if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue;
								if ((pp->m_type == PARTICLETYPE::VIRTUAL) or (p->m_type == PARTICLETYPE::VIRTUAL)) {
									if ((pp->m_type == PARTICLETYPE::VIRTUAL) and (p->m_type == PARTICLETYPE::VIRTUAL)) continue;
									if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
								}
							}
						}
						for (auto pp : this->parent->leaves[1]->particles_inside) {
							if (pp->m_id > p->m_id) {
								if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue;
								if ((pp->m_type == PARTICLETYPE::VIRTUAL) or (p->m_type == PARTICLETYPE::VIRTUAL)) {
									if ((pp->m_type == PARTICLETYPE::VIRTUAL) and (p->m_type == PARTICLETYPE::VIRTUAL)) continue;
									if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
								}
							}
						}
						for (auto pp : this->parent->leaves[2]->particles_inside) {
							if (pp->m_id > p->m_id) {
								if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue;
								if ((pp->m_type == PARTICLETYPE::VIRTUAL) or (p->m_type == PARTICLETYPE::VIRTUAL)) {
									if ((pp->m_type == PARTICLETYPE::VIRTUAL) and (p->m_type == PARTICLETYPE::VIRTUAL)) continue;
									if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
								}
							}
						}
						for (auto pp : this->parent->leaves[3]->particles_inside) {
							if (pp->m_id > p->m_id) {
								if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue;
								if ((pp->m_type == PARTICLETYPE::VIRTUAL) or (p->m_type == PARTICLETYPE::VIRTUAL)) {
									if ((pp->m_type == PARTICLETYPE::VIRTUAL) and (p->m_type == PARTICLETYPE::VIRTUAL)) continue;
									if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
								}
							}
						}
						for (Particle* nbpart : setOfNb) {
							if (nbpart->m_id > p->m_id) {
								//std::cout << p->m_id << "    " << nbpart->m_id << "\n";
								PartPairs.push_back(new ParticlePair(p, nbpart, dim)); //CREATING PAIRS
							}
						}
					}
					for (auto p : this->parent->leaves[1]->particles_inside) {
						std::set<Particle*> setOfNb;
						for (auto pp : this->parent->leaves[1]->particles_inside) {
							if (pp->m_id > p->m_id) {
								if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница  
								if ((pp->m_type == PARTICLETYPE::VIRTUAL) or (p->m_type == PARTICLETYPE::VIRTUAL)) {
									if ((pp->m_type == PARTICLETYPE::VIRTUAL) and (p->m_type == PARTICLETYPE::VIRTUAL)) continue;
									if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
								}
							}
						}
						for (auto nb_part : setOfNb) {
							if (nb_part->m_id > p->m_id) {
								//std::cout << p->m_id << "    " << nb_part->m_id << "\n";
								PartPairs.push_back(new ParticlePair(p, nb_part, dim)); //CREATING PAIRS
							}
						}
					}
					for (auto p : this->parent->leaves[2]->particles_inside) {
						std::set<Particle*> setOfNb;
						for (auto pp : this->parent->leaves[2]->particles_inside) {
							if (pp->m_id > p->m_id) {
								if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница  
								if ((pp->m_type == PARTICLETYPE::VIRTUAL) or (p->m_type == PARTICLETYPE::VIRTUAL)) {
									if ((pp->m_type == PARTICLETYPE::VIRTUAL) and (p->m_type == PARTICLETYPE::VIRTUAL)) continue;
									if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
								}
							}
						}
						for (auto nb_part : setOfNb) {
							if (nb_part->m_id > p->m_id) {
								//std::cout << p->m_id << "    " << nb_part->m_id << "\n";
								PartPairs.push_back(new ParticlePair(p, nb_part, dim)); //CREATING PAIRS
							}
						}
					}
					for (auto p : this->parent->leaves[3]->particles_inside) {
						std::set<Particle*> setOfNb;
						for (auto pp : this->parent->leaves[3]->particles_inside) {
							if (pp->m_id > p->m_id) {
								if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница  
								if ((pp->m_type == PARTICLETYPE::VIRTUAL) or (p->m_type == PARTICLETYPE::VIRTUAL)) {
									if ((pp->m_type == PARTICLETYPE::VIRTUAL) and (p->m_type == PARTICLETYPE::VIRTUAL)) continue;
									if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
								}
							}
						}
						for (auto nb_part : setOfNb) {
							if (nb_part->m_id > p->m_id) {
								//std::cout << p->m_id << "    " << nb_part->m_id << "\n";
								PartPairs.push_back(new ParticlePair(p, nb_part, dim)); //CREATING PAIRS
							}
						}
					}
					break;
				default:
					break;
				}
			}
			else {
				for (auto p : particles_inside) {
					std::set<Particle*> setOfNb;
					switch (this->m_id) {
					case(ID_0_0):
						for (auto pp : particles_inside) {
							if (pp->m_id > p->m_id) {
								if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue;
								if ((pp->m_type == PARTICLETYPE::VIRTUAL) or (p->m_type == PARTICLETYPE::VIRTUAL)) {
									if ((pp->m_type == PARTICLETYPE::VIRTUAL) and (p->m_type == PARTICLETYPE::VIRTUAL)) continue;
									if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
								}
							}
						}
						for (auto pp : this->parent->leaves[1]->particles_inside) {
							if (pp->m_id > p->m_id) {
								if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue;
								if ((pp->m_type == PARTICLETYPE::VIRTUAL) or (p->m_type == PARTICLETYPE::VIRTUAL)) {
									if ((pp->m_type == PARTICLETYPE::VIRTUAL) and (p->m_type == PARTICLETYPE::VIRTUAL)) continue;
									if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
								}
							}
						}
						for (auto pp : this->parent->leaves[2]->particles_inside) {
							if (pp->m_id > p->m_id) {
								if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue;
								if ((pp->m_type == PARTICLETYPE::VIRTUAL) or (p->m_type == PARTICLETYPE::VIRTUAL)) {
									if ((pp->m_type == PARTICLETYPE::VIRTUAL) and (p->m_type == PARTICLETYPE::VIRTUAL)) continue;
									if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
								}
							}
						}
						for (auto pp : this->parent->leaves[3]->particles_inside) {
							if (pp->m_id > p->m_id) {
								if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue;
								if ((pp->m_type == PARTICLETYPE::VIRTUAL) or (p->m_type == PARTICLETYPE::VIRTUAL)) {
									if ((pp->m_type == PARTICLETYPE::VIRTUAL) and (p->m_type == PARTICLETYPE::VIRTUAL)) continue;
									if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
								}
							}
						}
						break;
					case(ID_0_1):
						for (auto pp : particles_inside) {
							if (pp->m_id > p->m_id) {
								if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue;
								if ((pp->m_type == PARTICLETYPE::VIRTUAL) or (p->m_type == PARTICLETYPE::VIRTUAL)) {
									if ((pp->m_type == PARTICLETYPE::VIRTUAL) and (p->m_type == PARTICLETYPE::VIRTUAL)) continue;
									if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
								}
							}
						}
						for (auto pp : this->parent->parent->leaves[1]->leaves[0]->particles_inside) {
							if (pp->m_id > p->m_id) {
								if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue;
								if ((pp->m_type == PARTICLETYPE::VIRTUAL) or (p->m_type == PARTICLETYPE::VIRTUAL)) {
									if ((pp->m_type == PARTICLETYPE::VIRTUAL) and (p->m_type == PARTICLETYPE::VIRTUAL)) continue;
									if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
								}
							}
						}
						for (auto pp : this->parent->parent->leaves[1]->leaves[2]->particles_inside) {
							if (pp->m_id > p->m_id) {
								if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue;
								if ((pp->m_type == PARTICLETYPE::VIRTUAL) or (p->m_type == PARTICLETYPE::VIRTUAL)) {
									if ((pp->m_type == PARTICLETYPE::VIRTUAL) and (p->m_type == PARTICLETYPE::VIRTUAL)) continue;
									if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
								}
							}
						}
						for (auto pp : this->parent->leaves[3]->particles_inside) {
							if (pp->m_id > p->m_id) {
								if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue;
								if ((pp->m_type == PARTICLETYPE::VIRTUAL) or (p->m_type == PARTICLETYPE::VIRTUAL)) {
									if ((pp->m_type == PARTICLETYPE::VIRTUAL) and (p->m_type == PARTICLETYPE::VIRTUAL)) continue;
									if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
								}
							}
						}
						break;
					case(ID_1_0):
						for (auto pp : particles_inside) {
							if (pp->m_id > p->m_id) {
								if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue;
								if ((pp->m_type == PARTICLETYPE::VIRTUAL) or (p->m_type == PARTICLETYPE::VIRTUAL)) {
									if ((pp->m_type == PARTICLETYPE::VIRTUAL) and (p->m_type == PARTICLETYPE::VIRTUAL)) continue;
									if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
								}
							}
						}
						for (auto pp : this->parent->parent->leaves[2]->leaves[0]->particles_inside) {
							if (pp->m_id > p->m_id) {
								if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue;
								if ((pp->m_type == PARTICLETYPE::VIRTUAL) or (p->m_type == PARTICLETYPE::VIRTUAL)) {
									if ((pp->m_type == PARTICLETYPE::VIRTUAL) and (p->m_type == PARTICLETYPE::VIRTUAL)) continue;
									if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
								}
							}
						}
						for (auto pp : this->parent->parent->leaves[2]->leaves[1]->particles_inside) {
							if (pp->m_id > p->m_id) {
								if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue;
								if ((pp->m_type == PARTICLETYPE::VIRTUAL) or (p->m_type == PARTICLETYPE::VIRTUAL)) {
									if ((pp->m_type == PARTICLETYPE::VIRTUAL) and (p->m_type == PARTICLETYPE::VIRTUAL)) continue;
									if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
								}
							}
						}
						for (auto pp : this->parent->leaves[3]->particles_inside) {
							if (pp->m_id > p->m_id) {
								if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue;
								if ((pp->m_type == PARTICLETYPE::VIRTUAL) or (p->m_type == PARTICLETYPE::VIRTUAL)) {
									if ((pp->m_type == PARTICLETYPE::VIRTUAL) and (p->m_type == PARTICLETYPE::VIRTUAL)) continue;
									if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
								}
							}
						}
						break;
					case(ID_1_1):
						for (auto pp : particles_inside) {
							if (pp->m_id > p->m_id) {
								if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue;
								if ((pp->m_type == PARTICLETYPE::VIRTUAL) or (p->m_type == PARTICLETYPE::VIRTUAL)) {
									if ((pp->m_type == PARTICLETYPE::VIRTUAL) and (p->m_type == PARTICLETYPE::VIRTUAL)) continue;
									if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
								}
							}
						}
						for (auto pp : this->parent->parent->leaves[1]->leaves[2]->particles_inside) {
							if (pp->m_id > p->m_id) {
								if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue;
								if ((pp->m_type == PARTICLETYPE::VIRTUAL) or (p->m_type == PARTICLETYPE::VIRTUAL)) {
									if ((pp->m_type == PARTICLETYPE::VIRTUAL) and (p->m_type == PARTICLETYPE::VIRTUAL)) continue;
									if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
								}
							}
						}
						for (auto pp : this->parent->parent->leaves[2]->leaves[1]->particles_inside) {
							if (pp->m_id > p->m_id) {
								if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue;
								if ((pp->m_type == PARTICLETYPE::VIRTUAL) or (p->m_type == PARTICLETYPE::VIRTUAL)) {
									if ((pp->m_type == PARTICLETYPE::VIRTUAL) and (p->m_type == PARTICLETYPE::VIRTUAL)) continue;
									if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
								}
							}
						}
						for (auto pp : this->parent->parent->leaves[3]->leaves[0]->particles_inside) {
							if (pp->m_id > p->m_id) {
								if ((pp->m_type == PARTICLETYPE::BOUNDARY) and (p->m_type == PARTICLETYPE::BOUNDARY)) continue;
								if ((pp->m_type == PARTICLETYPE::VIRTUAL) or (p->m_type == PARTICLETYPE::VIRTUAL)) {
									if ((pp->m_type == PARTICLETYPE::VIRTUAL) and (p->m_type == PARTICLETYPE::VIRTUAL)) continue;
									if (p->distance(pp) <= KernelCoef * (p->m_SmR + pp->m_SmR) / 2.0) setOfNb.insert(pp);
								}
							}
						}
						break;
					default:
						break;
					}
					for (Particle* nbpart : setOfNb) {
						if (nbpart->m_id > p->m_id){
							//std::cout << p->m_id << "    " << nbpart->m_id << "\n";
							PartPairs.push_back(new ParticlePair(p, nbpart, dim)); //CREATING PAIRS
						}
					}
				}
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
		delete leaves[0];
		delete leaves[1];
		delete leaves[2];
		delete leaves[3];
	}
};