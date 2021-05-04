#include "ComputationalDomain.h"

part_prec InitialSmR = 0.0125;
part_prec Dens0 = 1.0;

float getRandomNumber(float min, float max) {
	static const double fraction = 1.0 / (static_cast<double>(RAND_MAX) + 1.0);
	return (rand()*fraction*(max - min) + min);
}

glm::mat3x3 Xrot(float radians) {
	glm::mat3x3 Retval;
	Retval[0][0] = 1; Retval[0][1] = 0;            Retval[0][2] = 0;
	Retval[1][0] = 0; Retval[1][1] = cos(radians); Retval[1][2] = -sin(radians);
	Retval[2][0] = 0; Retval[2][1] = sin(radians); Retval[2][2] = cos(radians);
	return Retval;
};
glm::mat3x3 Yrot(float radians) {
	glm::mat3x3 Retval;
	Retval[0][0] = cos(radians);     Retval[0][1] = 0; Retval[0][2] = sin(radians);
	Retval[1][0] = 0;                Retval[1][1] = 1; Retval[1][2] = 0;
	Retval[2][0] = 0 - sin(radians); Retval[2][1] = 0; Retval[2][2] = cos(radians);
	return Retval;
};
glm::mat3x3 Zrot(float radians) {
	glm::mat3x3 Retval;
	Retval[0][0] = cos(radians); Retval[0][1] = -sin(radians); Retval[0][2] = 0;
	Retval[1][0] = sin(radians); Retval[1][1] = cos(radians);  Retval[1][2] = 0;
	Retval[2][0] = 0;            Retval[2][1] = 0;             Retval[2][2] = 1;
	return Retval;
};

void SPH_CD::Initilization() {
	if (computationalDomain_mode == MODE_CD::CD_DEBUG)	{
		std::cout << "SPH_CD::Initilization\n";
	}
	
	int Nz = 1;
	int Ny = int(sqrt(m_options.nrOfParticles[PARTICLETYPE::REAL]));
	int Nx = int(sqrt(m_options.nrOfParticles[PARTICLETYPE::REAL]));

	float dx = (getXmax() - getXmin()) / static_cast<float>(Nx);
	float dy = (getYmax() - getYmin()) / static_cast<float>(Ny);
	float dz = (getZmax() - getZmin()) / static_cast<float>(Nz);

	dx = (1.0) / static_cast<float>(Nx);
	dy = (1.0) / static_cast<float>(Ny);
	dz = (1.0) / static_cast<float>(Nz);

	m_options.average_dim_steps.x = dx;
	m_options.average_dim_steps.y = dy;
	m_options.average_dim_steps.z = dz;

	float Xstart = getXmin() + dx / 2;
	float Ystart = getYmin() + dy / 2;
	float Zstart = getZmin();

	//Xstart = 0.25f;
	//Ystart = 0.25f;


	int nx = 0;
	int ny = 0;
	int nz = 0;
	
	part_prec part_x;
	part_prec part_y;
	part_prec part_z;
	for (int i = 0;i < m_options.nrOfParticles[PARTICLETYPE::REAL];i++) {

		switch (m_options.distributionREAL) {
		case(RANDOM):
			part_x = getRandomNumber(getXmin(), getXmax());
			if (nrOfDim > D1) part_y = getRandomNumber(getYmin(), getYmax());
			else part_y = getRandomNumber(0.f, 0.f);

			if (nrOfDim > D2) part_z = getRandomNumber(getZmin(), getZmax());
			else part_z = getRandomNumber(0.f, 0.f);
			break;
		case(UNIFORM):
			Zstart = getZmin();

			part_x = Xstart + (nx) * dx;
			if (nrOfDim > D1) part_y = Ystart + (ny ) * dy;
			else part_y = 0.f;
			if (nrOfDim > D2) part_z = Zstart + (nz) * dz;
			else part_z = 0.f;
			nx++;
			if (nx == Nx) {
				nx = 0;
				ny++;
				if (ny == Ny) {
					ny = 0;
					nz++;
					if (nz == Nz) {
						m_options.nrOfParticles[PARTICLETYPE::REAL] = i+1;
						break;
					}
				}
			}
			break;
		default:
			break;
		}


		float smR = InitialSmR;
		cd_prec density = Dens0;
		cd_prec mass = density*m_options.average_dim_steps.x*m_options.average_dim_steps.y*m_options.average_dim_steps.z;

		SPH.Particles.push_back(new Particle(i, PARTICLETYPE::REAL,glm::vec3(part_x, part_y, part_z), glm::vec3(0.0f, 0.f, 0.f), smR, density, mass));
		//std::cout << part_x << ", " << part_y << ", " << part_z << "\n";
	}


	std::cout << "SPH_CD::Initilization::REAL_PARTICLES_WERE_CREATED::" << m_options.nrOfParticles[PARTICLETYPE::REAL] << "\n";

	glm::vec3 AverageNormal;
	int nrOfParticlesPerBoundary = static_cast<int>(m_options.nrOfParticles[PARTICLETYPE::BOUNDARY]/Boundaries.size());
	int BoundaryParticleCount =0;
	for (auto*& i : Boundaries) {
		int nrOfPolygons = static_cast<int>(i->ReturnMesh()->getNumberOfIndice()) / 3;

		glm::vec3 OXYnormal = { 0.f,0.f,1.f };
		float PlaneCoefA;
		float PlaneCoefB;
		float PlaneCoefC;
		float PlaneCoefD;

		glm::vec3 Plane1;
		glm::vec3 Plane2;
		glm::vec3 Plane3;

		for (int p = 0;p < nrOfPolygons;p++) {
			//WE WANT TO ROTATE EVERY POLYGON IN A WAY, THAT TWO OF HIS SIDES COULD BE ANALISED AS LINEAR FUNCTIONS 
			GLuint* indexarray = i->ReturnMesh()->getIndexArray();
			int p1 = indexarray[p * 3 + 0];
			int p2 = indexarray[p * 3 + 1];
			int p3 = indexarray[p * 3 + 2];
			//std::cout << p1 << ", " << p2 << ", " << p3 << "\n";
			Vertex* vertexArray = i->ReturnMesh()->getVertexArray();
			glm::vec3 pos1 = vertexArray[p1].position;
			glm::vec3 pos2 = vertexArray[p2].position;
			glm::vec3 pos3 = vertexArray[p3].position;

			AverageNormal = (vertexArray[p1].normal + vertexArray[p2].normal + vertexArray[p3].normal);
			AverageNormal[0] /= 3;
			AverageNormal[1] /= 3;
			AverageNormal[2] /= 3;


			//std::cout << "pos1:{" << pos1.x << ", " << pos1.y << ", " << pos1.z << "},pos2:{" << pos2.x << ", " << pos2.y << ", " << pos2.z << "},pos3{" << pos3.x << ", " << pos3.y << ", " << pos3.z << "}\n";

			Plane1 = pos1;
			Plane2 = pos2;
			Plane3 = pos3;
			//Зная три точки, можно определить функцию плоскости
			PlaneCoefA = (Plane2.y - Plane1.y)*(Plane3.z - Plane1.z) - (Plane3.y - Plane1.y)*(Plane2.z - Plane1.z);
			PlaneCoefB = -(Plane2.x - Plane1.x)*(Plane3.z - Plane1.z) - (Plane3.x - Plane1.x)*(Plane2.z - Plane1.z);
			PlaneCoefC = (Plane2.x - Plane1.x)*(Plane3.y - Plane1.y) - (Plane3.x - Plane1.x)*(Plane2.y - Plane1.y);
			PlaneCoefD = -Plane1.x*PlaneCoefA - Plane1.y*PlaneCoefB - Plane1.z*PlaneCoefC;
			//std::cout << PlaneCoefA << "*x + " << PlaneCoefB << "*y + " << PlaneCoefC << "*z + " << PlaneCoefD << " =  0\n";
			//Далее определим, линию пересечение полученной плоскости с плоскостью OXY  (z=0):
			// x*PlaneCoefA + y*PlaneCoefB + PlaneCoefD = 0
			//Определим пересечение этой линии с каждой линией, из которых состоит полигон

			glm::vec3 interP;
			std::vector<glm::vec3> PolyLine0 = { pos1 ,pos2 ,pos3 };
			std::vector<glm::vec3> PolyLine1 = { pos2 ,pos3 ,pos1 };

			std::vector<glm::vec3> interPoints;


			float t;
			float a;
			for (int side = 0;side < 3;side++) {
				if(PlaneCoefA != 0){
					a = -PolyLine0[side].z / (PolyLine1[side].z - PolyLine0[side].z);
					t = a * (PolyLine1[side].y - PolyLine0[side].y) + PolyLine0[side].y;
					if (-PlaneCoefD / PlaneCoefA - t * PlaneCoefB / PlaneCoefA == a*(PolyLine1[side].x - PolyLine0[side].x) + PolyLine0[side].x){

					interP.x = -PlaneCoefD / PlaneCoefA - t * PlaneCoefB / PlaneCoefA;
					interP.y = t;
					interP.z = 0.f;
					interPoints.push_back(interP);
					}
				}
				else {
					a = -PolyLine0[side].z / (PolyLine1[side].z - PolyLine0[side].z);
					t = a * (PolyLine1[side].x - PolyLine0[side].x) + PolyLine0[side].x;
					if (-PlaneCoefD / PlaneCoefB - t * PlaneCoefA / PlaneCoefB == a * (PolyLine1[side].y - PolyLine0[side].y) + PolyLine0[side].y) {

						interP.x = t;
						interP.y = -PlaneCoefD / PlaneCoefB - t * PlaneCoefA / PlaneCoefB;
						interP.z = 0.f;
						interPoints.push_back(interP);
					}
				}
			}
			//std::cout << "cross_product: {" << glm::cross(OXYnormal, AverageNormal).x << ", " << glm::cross(OXYnormal, AverageNormal).y << ", " << glm::cross(OXYnormal, AverageNormal).z << "}\n";

			//std::cout << "GLOBAL COORDINATS\n";
			//std::cout << "pos1:{" << pos1.x << ", " << pos1.y << ", " << pos1.z << "},pos2:{" << pos2.x << ", " << pos2.y << ", " << pos2.z << "},pos3{" << pos3.x << ", " << pos3.y << ", " << pos3.z << "}\n";
			//std::cout << "Moving\n";
			glm::vec3 newCoordOrigin = pos1;
			pos1 += -newCoordOrigin;
			pos2 += -newCoordOrigin;
			pos3 += -newCoordOrigin;


			//std::cout << "pos1:{" << pos1.x << ", " << pos1.y << ", " << pos1.z << "},pos2:{" << pos2.x << ", " << pos2.y << ", " << pos2.z << "},pos3{" << pos3.x << ", " << pos3.y << ", " << pos3.z << "}\n";

			float fi_x, fi_y, fi_z = 0.f;
			bool rotationFlags[] = {false,false, false};
			// Z-ROTATION
			if ((pos2.y != 0) or (pos2.x < 0)) {
				fi_z = (asin(pos2.y / sqrt(pow(pos2.y, 2) + pow(pos2.x, 2))));
				//std::cout << glm::degrees(fi_z) << "\n";
				//if (pos2.x*pos2.y < 0) {
				//	fi_z = -fi_z;
				//}
				if ((pos2.x < 0))
					fi_z = M_PI - fi_z;
				//std::cout << glm::degrees(fi_z) << "\n";
				pos2 = Zrot(fi_z) * pos2;
				pos3 = Zrot(fi_z) * pos3;

				rotationFlags[2] = true;
			}
			//else { std::cout << "Z-rotation was skipped\n"; }
			//std::cout << "pos2.y->0   pos2:{" << pos2.x << ", " << pos2.y << ", " << pos2.z << "}\n";

			// Y-ROTATION
			if ((pos2.z != 0) or (pos2.x < 0)) {
				fi_y = (asin(pos2.z / sqrt(pow(pos2.x, 2) + pow(pos2.z, 2))));
				//std::cout << glm::degrees(fi_y) << "\n";
				//if (pos2.x*pos2.z < 0) {
				//	fi_y = -fi_y;
				//}
				if ((pos2.x < 0))
					fi_y += M_PI;
				//std::cout << glm::degrees(-fi_y) << "\n";
				pos2 = Yrot(-fi_y) * pos2;
				pos3 = Yrot(-fi_y) * pos3;

				rotationFlags[1] = true;
			}
			//else{ std::cout << "Y-rotation was skipped\n" ;}
			//std::cout << "pos2.z->0   pos2:{" << pos2.x << ", " << pos2.y << ", " << pos2.z << "}\n";

			// X-ROTATION
			if ((pos3.z != 0) or (pos3.y < 0)) {
				fi_x = (asin(pos3.z / sqrt(pow(pos3.z, 2) + pow(pos3.y, 2))));
				//std::cout << glm::degrees(fi_x) << "\n";
				if (pos3.y*pos3.z < 0) {
					fi_x = -fi_x;
				}
				if ((pos3.y < 0))
					fi_x += M_PI;
				//std::cout << glm::degrees(fi_x) << "\n";
				pos2 = Xrot(fi_x) * pos2;
				pos3 = Xrot(fi_x) * pos3;

				rotationFlags[0] = true;

			}

			//std::cout << "pos1:{" << pos1.x << ", " << pos1.y << ", " << pos1.z << "},pos2:{" << pos2.x << ", " << pos2.y << ", " << pos2.z << "},pos3{" << pos3.x << ", " << pos3.y << ", " << pos3.z << "}\n";

			//std::cout << "interPoints[0]" << interPoints[0].x << ", " << interPoints[0].y << ", " << interPoints[0].z << "\n";
			//std::cout << "interPoints[1]" << interPoints[1].x << ", " << interPoints[1].y << ", " << interPoints[1].z << "\n";
			for (auto &ip : interPoints) {
				//std::cout << "ip[] " << ip.x << ", " << ip.y << ", " << ip.z << "\n";
				ip += -newCoordOrigin;
				if (rotationFlags[2] == true)
					ip = Zrot(fi_z) * ip;
				if (rotationFlags[1] == true)
					ip = Yrot(-fi_y) * ip;
				if (rotationFlags[0] == true)
					ip = Xrot(fi_x) * ip;
				//std::cout << "ipnew[] " << ip.x << ", " << ip.y << ", " << ip.z << "\n";
			}
			//std::cout << "interPointsNew[0]" << interPoints[0].x << ", " << interPoints[0].y << ", " << interPoints[0].z << "\n";
			//std::cout << "interPointsNew[1]" << interPoints[1].x << ", " << interPoints[1].y << ", " << interPoints[1].z << "\n";

			//else { std::cout << "X-rotation was skipped\n"; }
			//std::cout << "pos3.z->0   pos3:{" << pos3.x << ", " << pos3.y << ", " << pos3.z << "}\n";


			//std::cout << "Rotation: "<< glm::degrees(fi_x) << ", " << glm::degrees(-fi_y) << ", " << glm::degrees(fi_z) << "\n" ;

			//std::cout << "rotating\n";
			//std::cout << "pos1:{" << pos1.x << ", " << pos1.y << ", " << pos1.z << "},pos2:{" << pos2.x << ", " << pos2.y << ", " << pos2.z << "},pos3{" << pos3.x << ", " << pos3.y << ", " << pos3.z << "}\n";
			LinearFunc LF1_3 = LinearFuncCoefficients(glm::vec2(pos1.x, pos1.y), glm::vec2(pos3.x, pos3.y));
			LinearFunc LF2_3 = LinearFuncCoefficients(glm::vec2(pos2.x, pos2.y), glm::vec2(pos3.x, pos3.y));

			LinearFunc LFIntersection = LinearFuncCoefficients(glm::vec2(interPoints[0].x, interPoints[0].y), glm::vec2(interPoints[1].x, interPoints[1].y));
			//std::cout << "LFIntersection(x) = " << LFIntersection.k << "*x + " << LFIntersection.b<<"\n";

			glm::vec3 mincoord = pos1;
			glm::vec3 maxcoord = pos1;
			
			// MAX AND MIN VALUES PER POLYGON
			for (int dim = 0;dim < 3;dim++) {
				if (mincoord[dim] > pos2[dim])
					mincoord[dim] = pos2[dim];
				if (mincoord[dim] > pos3[dim])
					mincoord[dim] = pos3[dim];
				if (maxcoord[dim] < pos2[dim])
					maxcoord[dim] = pos2[dim];
				if (maxcoord[dim] < pos3[dim])
					maxcoord[dim] = pos3[dim];
			}

			part_prec_3 rand_pos;
			//std::cout << "nrOfParticlesPerPolygon: " << static_cast<int>(nrOfParticlesPerBoundary / nrOfPolygons) <<" \n\n";


			float partB_x = 0.f;
			float partB_y = 0.f;
			float partB_z = 0.f;

			int BNx = int(sqrt(static_cast<int>(nrOfParticlesPerBoundary / nrOfPolygons)));
			int BNy = int(sqrt(static_cast<int>(nrOfParticlesPerBoundary / nrOfPolygons)));

							
			float Bdx = (maxcoord[0] - mincoord[0]) / static_cast<float>(BNx);
			float Bdy = (maxcoord[1] - mincoord[1]) / static_cast<float>(BNy);

			//FUTURE WARNING
			if (Bdx > Bdy*2)
				Bdx /= 2;
			else if (Bdy > Bdx*2)
				Bdy /= 2;
			//FUTURE WARNING



			float min2D_x = interPoints[0].x;
			float max2D_x = interPoints[0].x;
			float min2D_y = interPoints[0].y;
			float max2D_y = interPoints[0].y;

			//std::cout << "interPointsNew[0]" << interPoints[0].x << ", " << interPoints[0].y << ", " << interPoints[0].z << "\n";
			//std::cout << "interPointsNew[1]" << interPoints[1].x << ", " << interPoints[1].y << ", " << interPoints[1].z << "\n";
			//std::cout << "test" << min2D_x << ", " << max2D_x << "\n";
			if (nrOfDim == D2) {
				if (min2D_x > interPoints[1].x) {
					min2D_x = interPoints[1].x;
				}
				if (max2D_x < interPoints[1].x) {
					max2D_x = interPoints[1].x;
				}
				if (min2D_y > interPoints[1].y) {
					min2D_y = interPoints[1].y;
				}
				if (max2D_y < interPoints[1].y) {
					max2D_y = interPoints[1].y;
				}
			}
			//std::cout << "test_x " << min2D_x << ", " << max2D_x << "\n";
			//std::cout << "test_y " << min2D_y << ", " << max2D_y << "\n";

			if (nrOfDim == D2) {
				BNx = static_cast<int>(nrOfParticlesPerBoundary / nrOfPolygons);
				//std::cout << BNx << "\n";

				BNy = BNx;
				Bdx = (max2D_x - min2D_x) / (static_cast<float>(BNx) );
				Bdy = (max2D_y - min2D_y) / (static_cast<float>(BNy) );
				
			}




			int Bnx = 0;
			int Bny = 0;
			


			for (int j = 0;j < static_cast<int>(nrOfParticlesPerBoundary / nrOfPolygons);) {
				switch (m_options.distributionBOUNDARY) {
				case(RANDOM):
					if (nrOfDim == D2){
						partB_z = 0.f;
						partB_x = getRandomNumber(min2D_x,max2D_x);

						if (abs(min2D_x - max2D_x) < 10E-06){
							partB_y = getRandomNumber(min2D_y, max2D_y);
							//std::cout << "getRandomNumber(" << min2D_y << ", " << max2D_y << ") = " << partB_y << "\n";
						}
						else {
							partB_y = LFIntersection.k*partB_x + LFIntersection.b;
						}
						//std::cout << partB_x << ", " << partB_y << "\n";
					}
					if(nrOfDim == D3) {
						partB_x = getRandomNumber(mincoord[0], maxcoord[0]);
						partB_y = getRandomNumber(mincoord[1], maxcoord[1]);
						partB_z = 0.f;
					}
					break;
				case(UNIFORM):
					if (nrOfDim == D3) {
						Xstart = (partB_y - LF1_3.b) / LF1_3.k;
						Ystart = 0.f;
						partB_x = Xstart + Bnx * Bdx;
						partB_y = Ystart + Bny * Bdy;
						partB_z = 0;

						Bnx++;
						if (partB_x > (partB_y - LF2_3.b) / LF2_3.k) {
							Bnx = 0;
							Bny++;
							if (partB_y > maxcoord[1]) {
								partB_y = 0.f;
								j = static_cast<int>(nrOfParticlesPerBoundary / nrOfPolygons);
								break;
							}
						}
					}
					if (nrOfDim == D2) {
						partB_z = 0.f;
						Xstart = min2D_x + Bdx / 2.f;
						Ystart = min2D_y + Bdy / 2.f;
						partB_x = Xstart + Bnx * Bdx;


						if (abs(min2D_x - max2D_x) < 10E-06) {
							partB_y = Ystart + Bnx * Bdy;
							//std::cout << "getRandomNumber(" << min2D_y << ", " << max2D_y << ") = " << partB_y << "\n";
						}
						else {
							partB_y = LFIntersection.k*partB_x + LFIntersection.b;
						}
						//std::cout << partB_x << ", " << partB_y << "\n";
						Bnx++;
					}

					break;
				default:
					break;
				}

				rand_pos = { partB_x ,partB_y , partB_z };





				float smR = InitialSmR;
				cd_prec density = Dens0;
				cd_prec mass = density * m_options.average_dim_steps.x*m_options.average_dim_steps.y*m_options.average_dim_steps.z;


				//std::cout << rotationFlags[0] << rotationFlags[1] << rotationFlags[2] << "\n";
				if ((rand_pos.x < pos3.x) and (rand_pos.x > pos1.x)) {
					//std::cout << "[" << pos1.x << " , "<< pos3.x <<"]\n";
					//std::cout << "LF1_3(x) = " << LF1_3.k << "*x + " << LF1_3.b<<"\n";
					if (rand_pos.y < LF1_3.k*rand_pos.x + LF1_3.b) {
						
						if (rotationFlags[0])
							rand_pos = Xrot(-fi_x)* rand_pos;
						if (rotationFlags[1])
							rand_pos = Yrot(fi_y)* rand_pos;
						if (rotationFlags[2])
							rand_pos = Zrot(-fi_z)* rand_pos;
						rand_pos += newCoordOrigin;

						//std::cout << i->getBCval()[BCPARAMETR::VELOCITY_X] << "\n";
						glm::vec3 bc_velocity;

						if (i->isPeriodic()) {
							bc_velocity = glm::vec3(0.0, 0.0, 0.0);
						}else{
							bc_velocity = glm::vec3(i->getBCvals()[BCPARAMETER::VELOCITY_X], i->getBCvals()[BCPARAMETER::VELOCITY_Y], i->getBCvals()[BCPARAMETER::VELOCITY_Z]);
						}
						SPH.Particles.push_back(new Particle(m_options.nrOfParticles[PARTICLETYPE::REAL] + BoundaryParticleCount, PARTICLETYPE::BOUNDARY, rand_pos, bc_velocity, smR, density, mass));
						SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] + BoundaryParticleCount]->InitPolygonNormal = AverageNormal;
						SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] + BoundaryParticleCount]->assignToBoundary(i);

						//std::cout << m_options.nrOfParticles[PARTICLETYPE::REAL] + BoundaryParticleCount << "  " << rand_pos.x << ", " << rand_pos.y << ", " << rand_pos.z << "\n";
						j++; BoundaryParticleCount++;
					}
				}
				else if ((rand_pos.x > pos3.x) and (rand_pos.x < pos2.x)) {
					//std::cout << "[" << pos3.x << " , " << pos2.x << "]\n";
					//std::cout << "LF2_3(x) = " << LF2_3.k << "*x + " << LF2_3.b << "\n";
					if (rand_pos.y < LF2_3.k*rand_pos.x + LF2_3.b) {
						
						if (rotationFlags[0])
							rand_pos = Xrot(-fi_x)* rand_pos;
						if (rotationFlags[1])
							rand_pos = Yrot(fi_y)* rand_pos;
						if (rotationFlags[2])
							rand_pos = Zrot(-fi_z)* rand_pos;
						rand_pos += newCoordOrigin;
						

						glm::vec3 bc_velocity;
						if (i->isPeriodic()) {
							bc_velocity = glm::vec3(0.0, 0.0, 0.0);
						}
						else {
							bc_velocity = glm::vec3(i->getBCvals()[BCPARAMETER::VELOCITY_X], i->getBCvals()[BCPARAMETER::VELOCITY_Y], i->getBCvals()[BCPARAMETER::VELOCITY_Z]);
						}


						SPH.Particles.push_back(new Particle(m_options.nrOfParticles[PARTICLETYPE::REAL] + BoundaryParticleCount, PARTICLETYPE::BOUNDARY, rand_pos, bc_velocity, smR, density, mass));
						SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] + BoundaryParticleCount]->InitPolygonNormal = AverageNormal;
						SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] + BoundaryParticleCount]->assignToBoundary(i);
						//std::cout << m_options.nrOfParticles[PARTICLETYPE::REAL] + BoundaryParticleCount << "  " << rand_pos.x << ", " << rand_pos.y << ", " << rand_pos.z << "\n";
						j++; BoundaryParticleCount++;

					}
				}

			}
			//safeDelete(&vertexArray);
			vertexArray = nullptr;
			indexarray = nullptr;
			//if (vertexArray) delete vertexArray;
			//if (indexarray) delete indexarray;
			//std::cout << p << "   " << BoundaryParticleCount << "\n";
		}
	}
	m_options.nrOfParticles[PARTICLETYPE::BOUNDARY] = BoundaryParticleCount;
	std::cout << "SPH_CD::Initilization::BOUNDARY_PARTICLES_WERE_CREATED::"<< BoundaryParticleCount <<"\n";
	if (computationalDomain_mode == MODE_CD::CD_DEBUG) {
		std::cout << "	id:		first:" << SPH.Particles[0]->m_id << "	";
		std::cout << "	middle:" << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] / 2]->m_id << "	";
		std::cout << "last:" << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL]-1]->m_id << "\n";

		std::cout << "	position:	{" << SPH.Particles[0]->m_position.val.x << "," << SPH.Particles[0]->m_position.val.y << "," << SPH.Particles[0]->m_position.val.z << "} ";
		std::cout << "	{" << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] / 2]->m_position.val.x << "," << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] / 2]->m_position.val.y << "," << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] / 2]->m_position.val.z << "} ";
		std::cout << "	{" << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1]->m_position.val.x << "," << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1]->m_position.val.y << "," << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1]->m_position.val.z << "}\n";

		std::cout << "	velocity:	{" << SPH.Particles[0]->m_velocity.val.x << "," << SPH.Particles[0]->m_velocity.val.y << "," << SPH.Particles[0]->m_velocity.val.z << "} ";
		std::cout << "	{" << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] / 2]->m_velocity.val.x << "," << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] / 2]->m_velocity.val.y << "," << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] / 2]->m_velocity.val.z << "} ";
		std::cout << "	{" << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1]->m_velocity.val.x << "," << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1]->m_velocity.val.y << "," << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1]->m_velocity.val.z << "}\n";

		std::cout << "	density:	" << SPH.Particles[0]->m_density.val << "	";
		std::cout << "	" << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] / 2]->m_density.val << "	";
		std::cout << "	" << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1]->m_density.val << "\n";

		std::cout << "	mass:		" << SPH.Particles[0]->m_mass << "	";
		std::cout << "	" << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] / 2]->m_mass << "	";
		std::cout << "	" << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1]->m_mass << "\n";
	}
	std::cout << "\n";
	if (m_options.NBSAlg == NEIGBOURS_SEARCH_ALGORITHM::UNIFORM_GRID) {
		uniformGridPartitionInitialization();
	}

}

#define PARTICLE_IS(type) (i->m_type == PARTICLETYPE::type)
#define PARTICLE_IS_NOT(type) (i->m_type != PARTICLETYPE::type)
#define ANY_PARTICLE (true)
// NEIGBOUR PARTICLE IS 


void SPH_CD::PRB_refresh() {
	if (computationalDomain_mode == MODE_CD::CD_DEBUG) {
		std::cout << "SPH_CD::PRB_refresh\n";
	}
	PRB.clear();
	for (auto*& i : SPH.Particles) {
		//if (PARTICLE_IS(REAL)) {
		if (ANY_PARTICLE) {
			PRB.push_back({ i->m_position.val, i->m_color, static_cast<part_prec>(2.0*i->m_SmR/ InitialSmR* 0.025 *0.125)});
		}
	}
}

void SPH_CD::UpdateRendering(std::vector<Model*>* models, Texture* tex, Texture* tex_specualar, std::vector<Material*>* materials) {
	if (computationalDomain_mode == MODE_CD::CD_DEBUG) {
		std::cout << "SPH_CD::UpdateRendering\n";
	}
	//РЕНДЕРИНГ ВИРТУАЛЬНЫХ ЧАСТИЦ
	//for (auto*& i : (*models)[ModelId]->meshes) {
	//	delete i;
	//}
	//(*models)[ModelId]->meshes.resize(0);
	//for (int i = 0;i < PRB.size();i++) {
	//	//std::cout << i << "\n";
	//	(*models)[ModelId]->meshes.push_back(new Mesh(&Quad(glm::vec3(-0.5f), Quad::QuadNormal(Quad::Z, Quad::PLUS), 1.f, 1.f, PRB[i].color), PRB[i].position, glm::vec3(0.f), glm::vec3(0.f), glm::vec3(PRB[i].size)));
	//}
	//m_options.firstCycle = false;
	//std::cout << "models: " << models->size() << "\n";
	Mesh* mesh;
	Material* material;
	//for (int mat = 0;mat < materials->size(); mat++) {
	//	if(mat >= MaterialId)
	//		delete (*materials)[mat];
	//}
	//materials->resize(MaterialId);
	if (m_options.firstCycle == true) {
		m_options.firstCycle = false;
		//std::cout << "before: " << models->size() << "\n";
		//for (auto*& i : (*models)[ModelId]->meshes) {
		//	delete i;
		//}
		//(*models)[ModelId]->meshes.resize(0);
		for (int i = 0;i < PRB.size();i++) {
			//std::cout << i << "\n";
			//(*models)[ModelId]->meshes.push_back(new Mesh(&Quad(glm::vec3(-0.5f), Quad::QuadNormal(Quad::Z, Quad::PLUS), 1.f, 1.f, PRB[i].color), PRB[i].position, glm::vec3(0.f), glm::vec3(0.f), glm::vec3(PRB[i].size)));
			mesh = (new Mesh(&Quad(glm::vec3(-0.5f), Quad::QuadNormal(Quad::Z, Quad::PLUS), 1.f, 1.f, glm::vec3(0.3f)), PRB[i].position, glm::vec3(0.f), glm::vec3(0.f), glm::vec3(PRB[i].size)));
			material = (new Material(PRB[i].color, glm::vec3(0.9f), glm::vec3(1.f), 0, 1));
			(*models).push_back(new Model(PRB[i].position, material, tex, tex_specualar, mesh));
			(*models)[ModelId + i]->getMaterial()->ChangeLighting(PRB[i].color, glm::vec3(0.9f), glm::vec3(1.f));
			(*models)[ModelId + i]->moveTo(PRB[i].position);
			(*models)[ModelId + i]->scaleUpTo(glm::vec3(PRB[i].size));
			delete material;
			delete mesh;
		}
		//std::cout << "after: " << models->size() << "\n";

		assignNewModels(ModelId + m_options.nrOfParticles[REAL] + m_options.nrOfParticles[BOUNDARY]);
	}
	else {
		//std::cout << "before: " << models->size() << "\n";
		for (int m = 0;m < models->size();m++) {
			if (m >= CD_ModelId) {
				delete (*models)[m];
			}
		}
		models->resize(CD_ModelId);
		//std::cout << "mid: " << models->size() << "\n";
		for (int i = 0;i < PRB.size();i++) {
			if (i < m_options.nrOfParticles[REAL] + m_options.nrOfParticles[BOUNDARY]) {
				//(*models)[ModelId]->meshes[i]->changeColorTo(PRB[i].color);
				//(*models)[ModelId]->meshes[i]->moveTo(PRB[i].position);
				//(*models)[ModelId]->meshes[i]->changeScaleTo(glm::vec3(PRB[i].size));
				//(*models)[ModelId]->meshes[i]->update();
				(*models)[ModelId + i]->getMaterial()->ChangeLighting(PRB[i].color, glm::vec3(0.9f), glm::vec3(1.f));
				(*models)[ModelId + i]->moveTo(PRB[i].position);
				(*models)[ModelId + i]->scaleUpTo(glm::vec3(PRB[i].size));
			}
			else {
				//std::cout << PRB[i].position.x << " " << PRB[i].position.y << " " << PRB[i].position.z << "\n";
				//(*models)[ModelId]->meshes.push_back(new Mesh(&Quad(glm::vec3(-0.5f), Quad::QuadNormal(Quad::Z, Quad::PLUS), 1.f, 1.f, PRB[i].color), PRB[i].position, glm::vec3(0.f), glm::vec3(0.f), glm::vec3(PRB[i].size)));

				mesh = (new Mesh(&Quad(glm::vec3(-0.5f), Quad::QuadNormal(Quad::Z, Quad::PLUS), 1.f, 1.f, glm::vec3(0.3f)), PRB[i].position, glm::vec3(0.f), glm::vec3(0.f), glm::vec3(PRB[i].size)));
				//materials->push_back(new Material(PRB[i].color, glm::vec3(0.9f), glm::vec3(1.f), 0, 1));
				//(*models).push_back(new Model(PRB[i].position, (*materials)[materials->size() - 1], tex, tex_specualar, mesh));
				material = (new Material(PRB[i].color, glm::vec3(0.9f), glm::vec3(1.f), 0, 1));
				(*models).push_back(new Model(PRB[i].position, material, tex, tex_specualar, mesh));
				(*models)[ModelId + i]->getMaterial()->ChangeLighting(PRB[i].color, glm::vec3(0.9f), glm::vec3(1.f));
				(*models)[ModelId + i]->moveTo(PRB[i].position);
				(*models)[ModelId + i]->scaleUpTo(glm::vec3(PRB[i].size));
				delete material;
				delete mesh;
			}
		}
		//std::cout << "after: " << models->size() << "\n";
	}
}



void SPH_CD::DeletingVirtualParticles() {
	if (computationalDomain_mode == MODE_CD::CD_DEBUG) {
		std::cout << "SPH_CD::DeletingVirtualParticles\n";
	}
	deletingParticlePairs();
	if (m_options.boundary_handling == RENORMALIZATION) {

	}
	else {
		if (m_options.nrOfParticles[PARTICLETYPE::VIRTUAL] != 0) {
			for (int i = 0; i < m_options.nrOfParticles[PARTICLETYPE::VIRTUAL];i++) {
				delete SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] + m_options.nrOfParticles[PARTICLETYPE::BOUNDARY] + i];
			}
			SPH.Particles.resize(m_options.nrOfParticles[PARTICLETYPE::REAL] + m_options.nrOfParticles[PARTICLETYPE::BOUNDARY]);
			m_options.nrOfParticles[PARTICLETYPE::VIRTUAL] = 0;
		}
	}
}
void SPH_CD::AfterRendering(std::vector<Model*>* models) {
	//for (auto*& i : SPH.Particles) {
	//	i->refreshNeighbours();
	//}
}



void SPH_CD::timeStep_thread(cd_prec dt, std::atomic<bool>& dataReadyForRender, std::atomic<bool>& dataIsRendering) {
	if (computationalDomain_mode == MODE_CD::CD_DEBUG) {
		std::cout << "SPH_CD::timeStep\n";
	}
	//Saving values from previous step to render

	{
		DeletingVirtualParticles();
		BM.resetTeleportationData();
		RecalcPrep();
		//ColoringBtType();
		Coloring();
		dataReadyForRender.store(false);
		while (dataIsRendering) {}
		//std::cout << "SPH_CD::timeStepEnd::dataReadyForRender:" << dataReadyForRender << "\n";
		PRB_refresh();
		dataReadyForRender.store(true);
		//std::cout << "SPH_CD::timeStepEnd::dataReadyForRender:" << dataReadyForRender << "\n";
	}
	//Reseting needed values
	for (auto*& i : SPH.Particles) {
		i->refreshNeighbours();
		i->m_velocity.dval = part_prec_3(0.f);
	}

	if (m_options.firstCycle == false) {
		if (computationalDomain_mode == MODE_CD::CD_DEBUG) {
			std::cout << "SPH_CD::timeIntegration\n";
		}

		switch (m_options.timeIntegrationScheme) {
		case(EXPLICIT):
			////Полностью явная схема
			//DensityCalculation();
			//PressurePartCalculation();
			//ViscosityPartCalculation();
			//i->m_position.dval = i->m_velocity.val;
			//i->m_velocity.val += i->m_velocity.dval * dt;
			//i->m_position.val += i->m_position.dval * dt;
			//i->m_density.val += i->m_density.dval * dt;
		{
			DensityAndVelocityRecalculation();
			ExternalForces();
			deletingParticlePairs();
			for (auto*& i : SPH.Particles) {
				if (i->m_type == PARTICLETYPE::REAL) {
					//TELEPORTING
					part_prec_3 pos_prediction = i->m_velocity.val *static_cast<part_prec>(dt);
					if ((BM.TeleportationLinks.count(i->m_id) == 1) and (BM.TeleportationLinks[i->m_id].distanceToBoundary() != 0.0) and ((sqrt(pow(pos_prediction.x, 2) + pow(pos_prediction.y, 2) + pow(pos_prediction.z, 2))) != 0.0)) {
						double n_vel_cos_angle = ((pos_prediction.x*(-BM.TeleportationLinks[i->m_id].boundaryNormal().x)) + (pos_prediction.y*(-BM.TeleportationLinks[i->m_id].boundaryNormal().y)) + (pos_prediction.z*(-BM.TeleportationLinks[i->m_id].boundaryNormal().z))) / (sqrt(pow(pos_prediction.x, 2) + pow(pos_prediction.y, 2) + pow(pos_prediction.z, 2))*sqrt(pow(-BM.TeleportationLinks[i->m_id].boundaryNormal().x, 2) + pow(-BM.TeleportationLinks[i->m_id].boundaryNormal().y, 2) + pow(-BM.TeleportationLinks[i->m_id].boundaryNormal().z, 2)));
						if ((sqrt(pow(pos_prediction.x, 2) + pow(pos_prediction.y, 2) + pow(pos_prediction.z, 2))*n_vel_cos_angle > BM.TeleportationLinks[i->m_id].distanceToBoundary()) and (n_vel_cos_angle > 0.0))
							pos_prediction += BM.TeleportationLinks[i->m_id].teleportingDistance();
					}
					i->m_velocity.val += i->m_velocity.dval *static_cast<part_prec>(dt);
					i->m_position.val += pos_prediction;
					switch (m_options.densityChangeUsing) {
					case(VAL):
						i->m_density.val = i->m_density.dval;
						break;
					case(DVAL_DT):
						i->m_density.val += i->m_density.dval * dt;
						break;
					}
					i->p_art_water();
				}
			}
		}
		break;
		case(IMPLICIT):
			////Полностью неявная схема
			//while(condition){
			//	PressurePartCalculation();
			//	ViscosityPartCalculation();
			//	i->m_velocity.val += i->m_velocity.dval * dt;
			//	i->m_position.dval = i->m_velocity.val;
			//	i->m_position.val += i->m_position.dval * dt;
			//	DensityCalculation();
			//	i->m_density.val += i->m_density.dval * dt;
			//}
		{
			while (false) {
				VelocityRecalculation();
				ExternalForces();
				deletingParticlePairs();
				for (auto*& i : SPH.Particles) {
					if (i->m_type == PARTICLETYPE::REAL) {
						i->m_velocity.val += i->m_velocity.dval *static_cast<part_prec>(dt);
						//TELEPORTING
						part_prec_3 pos_prediction = i->m_velocity.val *static_cast<part_prec>(dt);
						if ((BM.TeleportationLinks.count(i->m_id) == 1) and (BM.TeleportationLinks[i->m_id].distanceToBoundary() != 0.0) and ((sqrt(pow(pos_prediction.x, 2) + pow(pos_prediction.y, 2) + pow(pos_prediction.z, 2))) != 0.0)) {
							double n_vel_cos_angle = ((pos_prediction.x*(-BM.TeleportationLinks[i->m_id].boundaryNormal().x)) + (pos_prediction.y*(-BM.TeleportationLinks[i->m_id].boundaryNormal().y)) + (pos_prediction.z*(-BM.TeleportationLinks[i->m_id].boundaryNormal().z))) / (sqrt(pow(pos_prediction.x, 2) + pow(pos_prediction.y, 2) + pow(pos_prediction.z, 2))*sqrt(pow(-BM.TeleportationLinks[i->m_id].boundaryNormal().x, 2) + pow(-BM.TeleportationLinks[i->m_id].boundaryNormal().y, 2) + pow(-BM.TeleportationLinks[i->m_id].boundaryNormal().z, 2)));
							if ((sqrt(pow(pos_prediction.x, 2) + pow(pos_prediction.y, 2) + pow(pos_prediction.z, 2))*n_vel_cos_angle > BM.TeleportationLinks[i->m_id].distanceToBoundary()) and (n_vel_cos_angle > 0.0))
								pos_prediction += BM.TeleportationLinks[i->m_id].teleportingDistance();
						}
						i->m_position.val += pos_prediction;
						i->m_velocity.dval = part_prec_3(0.f);
					}
				}
				DeletingVirtualParticles();
				BM.resetTeleportationData();
				RecalcPrep();
				DensityRecalculation();
				deletingParticlePairs();
				for (auto*& i : SPH.Particles) {
					if (i->m_type == PARTICLETYPE::REAL) {
						switch (m_options.densityChangeUsing) {
						case(VAL):
							i->m_density.val = i->m_density.dval;
							break;
						case(DVAL_DT):
							i->m_density.val += i->m_density.dval * dt;
							break;
						}
						i->p_art_water();
					}
				}
				DeletingVirtualParticles();
				BM.resetTeleportationData();
				RecalcPrep();
			}
		}
		break;
		case(SEMI_IMPICIT):
			////Полунеявная схема
			//PressurePartCalculation();
			//ViscosityPartCalculation();
			//i->m_velocity.val += i->m_velocity.dval * dt;
			//i->m_position.dval = i->m_velocity.val;
			//i->m_position.val += i->m_position.dval * dt;
			//DensityCalculation();
			//i->m_density.val += i->m_density.dval * dt;
		{
			VelocityRecalculation();
			ExternalForces();
			deletingParticlePairs();
			for (auto*& i : SPH.Particles) {
				if (i->m_type == PARTICLETYPE::REAL) {
					i->m_velocity.val += i->m_velocity.dval*static_cast<part_prec>(dt);
					//TELEPORTING
					part_prec_3 pos_prediction = i->m_velocity.val *static_cast<part_prec>(dt);
					if ((BM.TeleportationLinks.count(i->m_id) == 1) and (BM.TeleportationLinks[i->m_id].distanceToBoundary() != 0.0) and ((sqrt(pow(pos_prediction.x, 2) + pow(pos_prediction.y, 2) + pow(pos_prediction.z, 2))) != 0.0)) {
						double n_vel_cos_angle = ((pos_prediction.x*(-BM.TeleportationLinks[i->m_id].boundaryNormal().x)) + (pos_prediction.y*(-BM.TeleportationLinks[i->m_id].boundaryNormal().y)) + (pos_prediction.z*(-BM.TeleportationLinks[i->m_id].boundaryNormal().z))) / (sqrt(pow(pos_prediction.x, 2) + pow(pos_prediction.y, 2) + pow(pos_prediction.z, 2))*sqrt(pow(-BM.TeleportationLinks[i->m_id].boundaryNormal().x, 2) + pow(-BM.TeleportationLinks[i->m_id].boundaryNormal().y, 2) + pow(-BM.TeleportationLinks[i->m_id].boundaryNormal().z, 2)));
						if ((sqrt(pow(pos_prediction.x, 2) + pow(pos_prediction.y, 2) + pow(pos_prediction.z, 2))*n_vel_cos_angle > BM.TeleportationLinks[i->m_id].distanceToBoundary()) and (n_vel_cos_angle > 0.0))
							pos_prediction += BM.TeleportationLinks[i->m_id].teleportingDistance();
					}
					i->m_position.val += pos_prediction;
				}
			}
			DeletingVirtualParticles();
			BM.resetTeleportationData();
			RecalcPrep();
			DensityRecalculation();
			deletingParticlePairs();
			for (auto*& i : SPH.Particles) {
				if (i->m_type == PARTICLETYPE::REAL) {
					switch (m_options.densityChangeUsing) {
					case(VAL):
						i->m_density.val = i->m_density.dval;
						break;
					case(DVAL_DT):
						i->m_density.val += i->m_density.dval * dt;
						break;
					}
					i->p_art_water();
				}
			}
		}
		break;
		case(SECOND_ORDER_SCEME):
			////схема второго порядка
			//PressurePartCalculation();
			//ViscosityPartCalculation();
			//DensityCalculation();
			//i->m_position.dval = i->m_velocity.val;
			//i->m_position.val += i->m_position.dval * dt/2.0;
			// part_prec old_rho = i->m_density.val;
			//i->m_density.val += i->m_density.dval * dt/2.0;
			//DensityCalculation();
			//PressurePartCalculation();
			//ViscosityPartCalculation();
			//i->m_velocity.val += i->m_velocity.dval * dt;
			//i->m_position.dval = i->m_velocity.val;
			//i->m_position.val += i->m_position.dval * dt / 2.0;
			//	part_prec Eps_rho = -(i->m_density.dval/i->m_density.val)*dt
			//i->m_density.val  = old_rho*(2-Eps_rho)/(2+Eps_rho)
		{
			deletingParticlePairs();
			for (auto*& i : SPH.Particles) {
				if (i->m_type == PARTICLETYPE::REAL) {
					//TELEPORTING
					part_prec_3 pos_prediction = i->m_velocity.val *static_cast<part_prec>(dt / 2.0);
					if ((BM.TeleportationLinks.count(i->m_id) == 1) and (BM.TeleportationLinks[i->m_id].distanceToBoundary() != 0.0) and ((sqrt(pow(pos_prediction.x, 2) + pow(pos_prediction.y, 2) + pow(pos_prediction.z, 2))) != 0.0)) {
						double n_vel_cos_angle = ((pos_prediction.x*(-BM.TeleportationLinks[i->m_id].boundaryNormal().x)) + (pos_prediction.y*(-BM.TeleportationLinks[i->m_id].boundaryNormal().y)) + (pos_prediction.z*(-BM.TeleportationLinks[i->m_id].boundaryNormal().z))) / (sqrt(pow(pos_prediction.x, 2) + pow(pos_prediction.y, 2) + pow(pos_prediction.z, 2))*sqrt(pow(-BM.TeleportationLinks[i->m_id].boundaryNormal().x, 2) + pow(-BM.TeleportationLinks[i->m_id].boundaryNormal().y, 2) + pow(-BM.TeleportationLinks[i->m_id].boundaryNormal().z, 2)));
						if ((sqrt(pow(pos_prediction.x, 2) + pow(pos_prediction.y, 2) + pow(pos_prediction.z, 2))*n_vel_cos_angle > BM.TeleportationLinks[i->m_id].distanceToBoundary()) and (n_vel_cos_angle > 0.0))
							pos_prediction += BM.TeleportationLinks[i->m_id].teleportingDistance();
					}
					i->m_position.val += pos_prediction;
				}
			}
			BM.resetTeleportationData();
			DeletingVirtualParticles();
			RecalcPrep();
			DensityRecalculation();
			std::vector<part_prec> old_rho;
			old_rho.resize(SPH.Particles.size(), 0.f);
			VelocityRecalculation();
			ExternalForces();
			deletingParticlePairs();
			for (auto*& i : SPH.Particles) {
				if (i->m_type == PARTICLETYPE::REAL) {
					i->m_velocity.val += i->m_velocity.dval *static_cast<part_prec>(dt);
					//TELEPORTING
					part_prec_3 pos_prediction = i->m_velocity.val *static_cast<part_prec>(dt / 2.0);
					if ((BM.TeleportationLinks.count(i->m_id) == 1) and (BM.TeleportationLinks[i->m_id].distanceToBoundary() != 0.0) and ((sqrt(pow(pos_prediction.x, 2) + pow(pos_prediction.y, 2) + pow(pos_prediction.z, 2))) != 0.0)) {
						double n_vel_cos_angle = ((pos_prediction.x*(-BM.TeleportationLinks[i->m_id].boundaryNormal().x)) + (pos_prediction.y*(-BM.TeleportationLinks[i->m_id].boundaryNormal().y)) + (pos_prediction.z*(-BM.TeleportationLinks[i->m_id].boundaryNormal().z))) / (sqrt(pow(pos_prediction.x, 2) + pow(pos_prediction.y, 2) + pow(pos_prediction.z, 2))*sqrt(pow(-BM.TeleportationLinks[i->m_id].boundaryNormal().x, 2) + pow(-BM.TeleportationLinks[i->m_id].boundaryNormal().y, 2) + pow(-BM.TeleportationLinks[i->m_id].boundaryNormal().z, 2)));
						if ((sqrt(pow(pos_prediction.x, 2) + pow(pos_prediction.y, 2) + pow(pos_prediction.z, 2))*n_vel_cos_angle > BM.TeleportationLinks[i->m_id].distanceToBoundary()) and (n_vel_cos_angle > 0.0))
							pos_prediction += BM.TeleportationLinks[i->m_id].teleportingDistance();
					}
					i->m_position.val += pos_prediction;
					old_rho[i->m_id] = i->m_density.val;
					i->m_density.val += i->m_density.dval * dt / 2.0;
				}
			}
			DeletingVirtualParticles();
			RecalcPrep();
			DensityRecalculation();
			std::vector<part_prec> Eps_rho;
			Eps_rho.resize(SPH.Particles.size(), 0.f);
			deletingParticlePairs();
			for (auto*& i : SPH.Particles) {
				if (i->m_type == PARTICLETYPE::REAL) {
					Eps_rho[i->m_id] = -(i->m_density.dval / i->m_density.val)*dt;
					i->m_density.val = old_rho[i->m_id] * (2.0 - Eps_rho[i->m_id]) / (2.0 + Eps_rho[i->m_id]);
					i->p_art_water();
				}
			}
		}
		break;
		case(NONE):
			break;
		default:
			break;
		}

		if (computationalDomain_mode == MODE_CD::CD_DEBUG) {
			std::cout << "	id:		first:" << SPH.Particles[0]->m_id << "	";
			std::cout << "	middle:" << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] / 2]->m_id << "	";
			std::cout << "last:" << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1]->m_id << "\n";

			std::cout << "	position:	{" << SPH.Particles[0]->m_position.val.x << "," << SPH.Particles[0]->m_position.val.y << "," << SPH.Particles[0]->m_position.val.z << "} ";
			std::cout << "	{" << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] / 2]->m_position.val.x << "," << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] / 2]->m_position.val.y << "," << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] / 2]->m_position.val.z << "} ";
			std::cout << "	{" << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1]->m_position.val.x << "," << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1]->m_position.val.y << "," << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1]->m_position.val.z << "}\n";

			std::cout << "	velocity:	{" << SPH.Particles[0]->m_velocity.val.x << "," << SPH.Particles[0]->m_velocity.val.y << "," << SPH.Particles[0]->m_velocity.val.z << "} ";
			std::cout << "	{" << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] / 2]->m_velocity.val.x << "," << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] / 2]->m_velocity.val.y << "," << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] / 2]->m_velocity.val.z << "} ";
			std::cout << "	{" << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1]->m_velocity.val.x << "," << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1]->m_velocity.val.y << "," << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1]->m_velocity.val.z << "}\n";

			std::cout << "	density:	" << SPH.Particles[0]->m_density.val << "	";
			std::cout << "	" << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] / 2]->m_density.val << "	";
			std::cout << "	" << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1]->m_density.val << "\n";

			std::cout << "	mass:		" << SPH.Particles[0]->m_mass << "	";
			std::cout << "	" << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] / 2]->m_mass << "	";
			std::cout << "	" << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1]->m_mass << "\n";
		}
	}
	else {
	}


	//Coloring();
	//RecalcPrep();
	//deletingParticlePairs();
	//PRB_refresh();


	//if (computationalDomain_mode == MODE_CD::CD_DEBUG) { std::cout << "SPH_CD::Deleting_Pairs\n"; }
	//deletingVirtualParticles();
	// ТЕБЯ МОЖНО И ПО РАНЬШЕ

	SaveMaxVelocity();

	if (computationalDomain_mode == MODE_CD::CD_DEBUG) { std::cout << "SPH_CD::timeStepEnd\n"; }
	timeStepEnd();

	if (computationalDomain_mode == MODE_CD::CD_DEBUG) {
		std::cout << "Current Time: " << getCurrentTime() << "\n";
		std::cout << "current timeStep: " << getDeltaTime() << "   minimal time scale: " << getStatMinTime() << "\n";
	}
	if (getStatMinTime() < getDeltaTime()) { setDeltaTime(getStatMinTime()); }
	else { setDeltaTime(getInitialDeltaTime()); }


}




void SPH_CD::timeStep(cd_prec dt) {
	if (computationalDomain_mode == MODE_CD::CD_DEBUG) {
		std::cout << "SPH_CD::timeStep\n";
	}

	//Saving values from previous step to render
	{
		DeletingVirtualParticles();
		BM.resetTeleportationData();
		RecalcPrep();
		//ColoringBtType();
		Coloring();
		PRB_refresh();
	}
	//Reseting needed values
	for (auto*& i : SPH.Particles) {
		i->m_velocity.dval = part_prec_3(0.f);
	}




	if (m_options.firstCycle == false) {
		if (computationalDomain_mode == MODE_CD::CD_DEBUG) {
			std::cout << "SPH_CD::timeIntegration\n";
		}

		switch (m_options.timeIntegrationScheme) {
		case(EXPLICIT): 
		////Полностью явная схема
		//DensityCalculation();
		//PressurePartCalculation();
		//ViscosityPartCalculation();
		//i->m_position.dval = i->m_velocity.val;
		//i->m_velocity.val += i->m_velocity.dval * dt;
		//i->m_position.val += i->m_position.dval * dt;
		//i->m_density.val += i->m_density.dval * dt;
		{
			DensityAndVelocityRecalculation();
			ExternalForces();
			deletingParticlePairs();
			for (auto*& i : SPH.Particles) {
				if (i->m_type == PARTICLETYPE::REAL) {
					//TELEPORTING
					part_prec_3 pos_prediction = i->m_velocity.val *static_cast<part_prec>(dt);
					if ((BM.TeleportationLinks.count(i->m_id) == 1) and (BM.TeleportationLinks[i->m_id].distanceToBoundary() != 0.0) and ((sqrt(pow(pos_prediction.x, 2) + pow(pos_prediction.y, 2) + pow(pos_prediction.z, 2))) != 0.0)) {
						double n_vel_cos_angle = ((pos_prediction.x*(-BM.TeleportationLinks[i->m_id].boundaryNormal().x)) + (pos_prediction.y*(-BM.TeleportationLinks[i->m_id].boundaryNormal().y)) + (pos_prediction.z*(-BM.TeleportationLinks[i->m_id].boundaryNormal().z))) / (sqrt(pow(pos_prediction.x, 2) + pow(pos_prediction.y, 2) + pow(pos_prediction.z, 2))*sqrt(pow(-BM.TeleportationLinks[i->m_id].boundaryNormal().x, 2) + pow(-BM.TeleportationLinks[i->m_id].boundaryNormal().y, 2) + pow(-BM.TeleportationLinks[i->m_id].boundaryNormal().z, 2)));
						if ((sqrt(pow(pos_prediction.x, 2) + pow(pos_prediction.y, 2) + pow(pos_prediction.z, 2))*n_vel_cos_angle > BM.TeleportationLinks[i->m_id].distanceToBoundary()) and (n_vel_cos_angle > 0.0))
							pos_prediction += BM.TeleportationLinks[i->m_id].teleportingDistance();
					}
					i->m_velocity.val += i->m_velocity.dval *static_cast<part_prec>(dt);
					i->m_position.val += pos_prediction;
					switch (m_options.densityChangeUsing) {
					case(VAL):
						i->m_density.val = i->m_density.dval;
						break;
					case(DVAL_DT):
						i->m_density.val += i->m_density.dval * dt;
						break;
					}
					i->p_art_water();
				}
			}
		}
			break;
		case(IMPLICIT):
		////Полностью неявная схема
		//while(condition){
		//	PressurePartCalculation();
		//	ViscosityPartCalculation();
		//	i->m_velocity.val += i->m_velocity.dval * dt;
		//	i->m_position.dval = i->m_velocity.val;
		//	i->m_position.val += i->m_position.dval * dt;
		//	DensityCalculation();
		//	i->m_density.val += i->m_density.dval * dt;
		//}
		{
			while (false) {
				VelocityRecalculation();
				ExternalForces();
				deletingParticlePairs();
				for (auto*& i : SPH.Particles) {
					if (i->m_type == PARTICLETYPE::REAL) {
						i->m_velocity.val += i->m_velocity.dval *static_cast<part_prec>(dt);
						//TELEPORTING
						part_prec_3 pos_prediction = i->m_velocity.val *static_cast<part_prec>(dt);
						if ((BM.TeleportationLinks.count(i->m_id) == 1) and (BM.TeleportationLinks[i->m_id].distanceToBoundary() != 0.0) and ((sqrt(pow(pos_prediction.x, 2) + pow(pos_prediction.y, 2) + pow(pos_prediction.z, 2))) != 0.0)) {
							double n_vel_cos_angle = ((pos_prediction.x*(-BM.TeleportationLinks[i->m_id].boundaryNormal().x)) + (pos_prediction.y*(-BM.TeleportationLinks[i->m_id].boundaryNormal().y)) + (pos_prediction.z*(-BM.TeleportationLinks[i->m_id].boundaryNormal().z))) / (sqrt(pow(pos_prediction.x, 2) + pow(pos_prediction.y, 2) + pow(pos_prediction.z, 2))*sqrt(pow(-BM.TeleportationLinks[i->m_id].boundaryNormal().x, 2) + pow(-BM.TeleportationLinks[i->m_id].boundaryNormal().y, 2) + pow(-BM.TeleportationLinks[i->m_id].boundaryNormal().z, 2)));
							if ((sqrt(pow(pos_prediction.x, 2) + pow(pos_prediction.y, 2) + pow(pos_prediction.z, 2))*n_vel_cos_angle > BM.TeleportationLinks[i->m_id].distanceToBoundary()) and (n_vel_cos_angle > 0.0))
								pos_prediction += BM.TeleportationLinks[i->m_id].teleportingDistance();
						}
						i->m_position.val += pos_prediction;
						i->m_velocity.dval = part_prec_3(0.f);
					}
				}
				DeletingVirtualParticles();
				BM.resetTeleportationData();
				RecalcPrep();
				DensityRecalculation();
				deletingParticlePairs();
				for (auto*& i : SPH.Particles) {
					if (i->m_type == PARTICLETYPE::REAL) {
						switch (m_options.densityChangeUsing) {
						case(VAL):
							i->m_density.val = i->m_density.dval;
							break;
						case(DVAL_DT):
							i->m_density.val += i->m_density.dval * dt;
							break;
						}
						i->p_art_water();
					}
				}
				DeletingVirtualParticles();
				BM.resetTeleportationData();
				RecalcPrep();
			}
		}
			break;
		case(SEMI_IMPICIT):
		////Полунеявная схема
		//PressurePartCalculation();
		//ViscosityPartCalculation();
		//i->m_velocity.val += i->m_velocity.dval * dt;
		//i->m_position.dval = i->m_velocity.val;
		//i->m_position.val += i->m_position.dval * dt;
		//DensityCalculation();
		//i->m_density.val += i->m_density.dval * dt;
		{
			VelocityRecalculation();
			ExternalForces();
			deletingParticlePairs();
			for (auto*& i : SPH.Particles) {
				if (i->m_type == PARTICLETYPE::REAL) {
					i->m_velocity.val += i->m_velocity.dval*static_cast<part_prec>(dt);
					//TELEPORTING
					part_prec_3 pos_prediction = i->m_velocity.val *static_cast<part_prec>(dt);
					if ((BM.TeleportationLinks.count(i->m_id) == 1) and (BM.TeleportationLinks[i->m_id].distanceToBoundary() != 0.0) and ((sqrt(pow(pos_prediction.x, 2) + pow(pos_prediction.y, 2) + pow(pos_prediction.z, 2))) != 0.0)) {
						double n_vel_cos_angle = ((pos_prediction.x*(-BM.TeleportationLinks[i->m_id].boundaryNormal().x)) + (pos_prediction.y*(-BM.TeleportationLinks[i->m_id].boundaryNormal().y)) + (pos_prediction.z*(-BM.TeleportationLinks[i->m_id].boundaryNormal().z))) / (sqrt(pow(pos_prediction.x, 2) + pow(pos_prediction.y, 2) + pow(pos_prediction.z, 2))*sqrt(pow(-BM.TeleportationLinks[i->m_id].boundaryNormal().x, 2) + pow(-BM.TeleportationLinks[i->m_id].boundaryNormal().y, 2) + pow(-BM.TeleportationLinks[i->m_id].boundaryNormal().z, 2)));
						if ((sqrt(pow(pos_prediction.x, 2) + pow(pos_prediction.y, 2) + pow(pos_prediction.z, 2))*n_vel_cos_angle > BM.TeleportationLinks[i->m_id].distanceToBoundary()) and (n_vel_cos_angle > 0.0))
							pos_prediction += BM.TeleportationLinks[i->m_id].teleportingDistance();
					}
					i->m_position.val += pos_prediction;
				}
			}
			DeletingVirtualParticles();
			BM.resetTeleportationData();
			RecalcPrep();
			DensityRecalculation();
			deletingParticlePairs();
			for (auto*& i : SPH.Particles) {
				if (i->m_type == PARTICLETYPE::REAL) {
					switch (m_options.densityChangeUsing) {
					case(VAL):
						i->m_density.val = i->m_density.dval;
						break;
					case(DVAL_DT):
						i->m_density.val += i->m_density.dval * dt;
						break;
					}
					i->p_art_water();
				}
			}
		}
			break;
		case(SECOND_ORDER_SCEME):
		////схема второго порядка
		//PressurePartCalculation();
		//ViscosityPartCalculation();
		//DensityCalculation();
		//i->m_position.dval = i->m_velocity.val;
		//i->m_position.val += i->m_position.dval * dt/2.0;
		// part_prec old_rho = i->m_density.val;
		//i->m_density.val += i->m_density.dval * dt/2.0;
		//DensityCalculation();
		//PressurePartCalculation();
		//ViscosityPartCalculation();
		//i->m_velocity.val += i->m_velocity.dval * dt;
		//i->m_position.dval = i->m_velocity.val;
		//i->m_position.val += i->m_position.dval * dt / 2.0;
		//	part_prec Eps_rho = -(i->m_density.dval/i->m_density.val)*dt
		//i->m_density.val  = old_rho*(2-Eps_rho)/(2+Eps_rho)
		{
			deletingParticlePairs();
			for (auto*& i : SPH.Particles) {
				if (i->m_type == PARTICLETYPE::REAL) {
					//TELEPORTING
					part_prec_3 pos_prediction = i->m_velocity.val *static_cast<part_prec>(dt / 2.0);
					if ((BM.TeleportationLinks.count(i->m_id) == 1) and (BM.TeleportationLinks[i->m_id].distanceToBoundary() != 0.0) and ((sqrt(pow(pos_prediction.x, 2) + pow(pos_prediction.y, 2) + pow(pos_prediction.z, 2))) != 0.0)) {
						double n_vel_cos_angle = ((pos_prediction.x*(-BM.TeleportationLinks[i->m_id].boundaryNormal().x)) + (pos_prediction.y*(-BM.TeleportationLinks[i->m_id].boundaryNormal().y)) + (pos_prediction.z*(-BM.TeleportationLinks[i->m_id].boundaryNormal().z))) / (sqrt(pow(pos_prediction.x, 2) + pow(pos_prediction.y, 2) + pow(pos_prediction.z, 2))*sqrt(pow(-BM.TeleportationLinks[i->m_id].boundaryNormal().x, 2) + pow(-BM.TeleportationLinks[i->m_id].boundaryNormal().y, 2) + pow(-BM.TeleportationLinks[i->m_id].boundaryNormal().z, 2)));
						if ((sqrt(pow(pos_prediction.x, 2) + pow(pos_prediction.y, 2) + pow(pos_prediction.z, 2))*n_vel_cos_angle > BM.TeleportationLinks[i->m_id].distanceToBoundary()) and (n_vel_cos_angle > 0.0))
							pos_prediction += BM.TeleportationLinks[i->m_id].teleportingDistance();
					}
					i->m_position.val += pos_prediction;
				}
			}
			BM.resetTeleportationData();
			DeletingVirtualParticles();
			RecalcPrep();
			DensityRecalculation();
			std::vector<part_prec> old_rho;
			old_rho.resize(SPH.Particles.size(), 0.f);
			VelocityRecalculation();
			ExternalForces();
			deletingParticlePairs();
			for (auto*& i : SPH.Particles) {
				if (i->m_type == PARTICLETYPE::REAL) {
					i->m_velocity.val += i->m_velocity.dval *static_cast<part_prec>(dt) ;
					//TELEPORTING
					part_prec_3 pos_prediction = i->m_velocity.val *static_cast<part_prec>(dt/2.0);
					if ((BM.TeleportationLinks.count(i->m_id) == 1) and (BM.TeleportationLinks[i->m_id].distanceToBoundary() != 0.0) and ((sqrt(pow(pos_prediction.x, 2) + pow(pos_prediction.y, 2) + pow(pos_prediction.z, 2))) != 0.0)) {
						double n_vel_cos_angle = ((pos_prediction.x*(-BM.TeleportationLinks[i->m_id].boundaryNormal().x)) + (pos_prediction.y*(-BM.TeleportationLinks[i->m_id].boundaryNormal().y)) + (pos_prediction.z*(-BM.TeleportationLinks[i->m_id].boundaryNormal().z))) / (sqrt(pow(pos_prediction.x, 2) + pow(pos_prediction.y, 2) + pow(pos_prediction.z, 2))*sqrt(pow(-BM.TeleportationLinks[i->m_id].boundaryNormal().x, 2) + pow(-BM.TeleportationLinks[i->m_id].boundaryNormal().y, 2) + pow(-BM.TeleportationLinks[i->m_id].boundaryNormal().z, 2)));
						if ((sqrt(pow(pos_prediction.x, 2) + pow(pos_prediction.y, 2) + pow(pos_prediction.z, 2))*n_vel_cos_angle > BM.TeleportationLinks[i->m_id].distanceToBoundary()) and (n_vel_cos_angle > 0.0))
							pos_prediction += BM.TeleportationLinks[i->m_id].teleportingDistance();
					}
					i->m_position.val += pos_prediction;
					old_rho[i->m_id] = i->m_density.val;
					i->m_density.val += i->m_density.dval * dt / 2.0;
				}
			}
			DeletingVirtualParticles();
			RecalcPrep();
			DensityRecalculation();
			std::vector<part_prec> Eps_rho;
			Eps_rho.resize(SPH.Particles.size(), 0.f);
			deletingParticlePairs();
			for (auto*& i : SPH.Particles) {
				if (i->m_type == PARTICLETYPE::REAL) {
					Eps_rho[i->m_id] = -(i->m_density.dval / i->m_density.val)*dt;
					i->m_density.val = old_rho[i->m_id] * (2.0 - Eps_rho[i->m_id]) / (2.0 + Eps_rho[i->m_id]);
					i->p_art_water();
				}
			}
		}
			break;
		case(NONE):
			break;
		default:
			break;
		}

		if (computationalDomain_mode == MODE_CD::CD_DEBUG) {
			std::cout << "	id:		first:" << SPH.Particles[0]->m_id << "	";
			std::cout << "	middle:" << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] / 2]->m_id << "	";
			std::cout << "last:" << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1]->m_id << "\n";

			std::cout << "	position:	{" << SPH.Particles[0]->m_position.val.x << "," << SPH.Particles[0]->m_position.val.y << "," << SPH.Particles[0]->m_position.val.z << "} ";
			std::cout << "	{" << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] / 2]->m_position.val.x << "," << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] / 2]->m_position.val.y << "," << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] / 2]->m_position.val.z << "} ";
			std::cout << "	{" << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1]->m_position.val.x << "," << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1]->m_position.val.y << "," << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1]->m_position.val.z << "}\n";

			std::cout << "	velocity:	{" << SPH.Particles[0]->m_velocity.val.x << "," << SPH.Particles[0]->m_velocity.val.y << "," << SPH.Particles[0]->m_velocity.val.z << "} ";
			std::cout << "	{" << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] / 2]->m_velocity.val.x << "," << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] / 2]->m_velocity.val.y << "," << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] / 2]->m_velocity.val.z << "} ";
			std::cout << "	{" << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1]->m_velocity.val.x << "," << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1]->m_velocity.val.y << "," << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1]->m_velocity.val.z << "}\n";

			std::cout << "	density:	" << SPH.Particles[0]->m_density.val << "	";
			std::cout << "	" << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] / 2]->m_density.val << "	";
			std::cout << "	" << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1]->m_density.val << "\n";

			std::cout << "	mass:		" << SPH.Particles[0]->m_mass << "	";
			std::cout << "	" << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] / 2]->m_mass << "	";
			std::cout << "	" << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1]->m_mass << "\n";
		}
	}
	else {
		m_options.firstCycle = false;
	}


	//Coloring();
	//RecalcPrep();
	//deletingParticlePairs();
	//PRB_refresh();


	//if (computationalDomain_mode == MODE_CD::CD_DEBUG) { std::cout << "SPH_CD::Deleting_Pairs\n"; }
	//deletingVirtualParticles();
	// ТЕБЯ МОЖНО И ПО РАНЬШЕ
	
	SaveMaxVelocity();

	if (computationalDomain_mode == MODE_CD::CD_DEBUG) { std::cout << "SPH_CD::timeStepEnd\n"; }
	timeStepEnd();
	std::cout << "Current Time: " << getCurrentTime() << "\n";
	std::cout << "current timeStep: " << getDeltaTime() << "   minimal time scale: " << getStatMinTime() << "\n";
	if (getStatMinTime() < getDeltaTime()) { setDeltaTime(getStatMinTime()); }
	else { setDeltaTime(getInitialDeltaTime()); }

}

bool SPH_CD::chekParticlePairsFor(Particle* p1, Particle* p2) {
	for (auto pair : SPH.ParticlePairs) {
		if (((pair->Mpart_ptr->m_id == p1->m_id) and (pair->NBpart_ptr->m_id == p2->m_id)) or ((pair->Mpart_ptr->m_id == p2->m_id) and (pair->NBpart_ptr->m_id == p1->m_id))) {
			return true;
		}
	}
	return false;
}

void SPH_CD::neighbourSearch() {
	if (computationalDomain_mode == MODE_CD::CD_DEBUG) { std::cout << "SPH_CD::neighbourSearch::PARTICLE_PAIRS_WERE_CREATED::" << SPH.ParticlePairs.size() << " --> "; }
	BM.resetBoundaryData();
	switch (m_options.NBSAlg) {
	case(NEIGBOURS_SEARCH_ALGORITHM::DIRECT):
		// Точно правильно
		for (int i = 0;i < SPH.Particles.size();i++) {  //  SPH.Particles.size()-1 ет.к. у последнего уже нету соседей
			Particle* i_p = SPH.Particles[i];
			//std::cout << i_p->m_id << "_th particle has neighbours: \n";
			for (int j = i + 1;j < SPH.Particles.size();j++) { //  SPH.Particles.size()-1 ет.к. у последнего уже нету соседей
				Particle* j_p = SPH.Particles[j];
				if (i_p == j_p) continue;
				if ((j_p->m_type == PARTICLETYPE::BOUNDARY) and (i_p->m_type == PARTICLETYPE::BOUNDARY)) continue; // Не учитываем пары Граница-Граница  
				if (i_p->distance(j_p) <= m_options.smoothingKernelLengthCoefficient * (i_p->m_SmR + j_p->m_SmR) / 2.0) SPH.ParticlePairs.push_back(new ParticlePair(i_p, j_p, nrOfDim)); //CREATING PAIRS
				j_p = nullptr;
			}
			//delete j_p;
			i_p = nullptr;
		}
		//delete i_p;
		break;
	case(NEIGBOURS_SEARCH_ALGORITHM::UNIFORM_GRID):
		for (int i = 0;i < SPH.Particles.size();i++) {  //  SPH.Particles.size()-1 ет.к. у последнего уже нету соседей
			Particle* i_p = SPH.Particles[i];
			UG->addPoint(i_p);
		}
		UG->CreatingParticlePairs(SPH.ParticlePairs, m_options.smoothingKernelLengthCoefficient, nrOfDim);
		//UG->clear();
		break;
	default:
		break;
	}

	// Mass recalculation
	if (computationalDomain_mode == MODE_CD::CD_DEBUG) { std::cout << "SPH_CD::neighbourSearch::PARTICLE_PAIRS_WERE_CREATED::" << SPH.ParticlePairs.size() << "\n"; }
	if (m_options.firstCycle == true) {
		float factor;
		float smR = m_options.smoothingKernelLengthCoefficient * InitialSmR;
		if(nrOfDim == D1) factor = (1.f / smR);
		else if (nrOfDim == D2) factor = (15.f / (7.f * M_PI*pow(smR, 2)));
		else if (nrOfDim == D3) factor = (3.f / (2.f * M_PI*pow(smR, 3)));

		float averageNrOfNeighbours = 0;
		int maxNrOfNeighbours = 0;
		int chosenParticleId;

		for (auto*& i : SPH.Particles) {
			if (i->m_type == REAL) {
				averageNrOfNeighbours += i->nrOfNeighbours;
				if (i->nrOfNeighbours > maxNrOfNeighbours) {
					maxNrOfNeighbours = i->nrOfNeighbours;
					chosenParticleId = i->m_id;
				}
			}
		}
		averageNrOfNeighbours = averageNrOfNeighbours / m_options.nrOfParticles[REAL];
		//averageNrOfNeighbours = 20.0;

		float aveW;
		aveW = factor / 4.f;
		aveW = factor * (2.f / 3.f ); //что это
		int counter = 0;
		float aveDistance = 0.f;
		for (auto*& i : SPH.ParticlePairs) {
			if ((i->Mpart_ptr->m_id == chosenParticleId) or (i->NBpart_ptr->m_type == REAL)) {
				aveW += i->W(i->Mpart_ptr);
				aveDistance += i->getPairDistance();
				counter++;
			}
			if ((i->NBpart_ptr->m_id == chosenParticleId) or (i->Mpart_ptr->m_type == REAL)) {
				aveW += i->W(i->NBpart_ptr);
				aveDistance += i->getPairDistance();
				counter++;
			}
		}
		aveW = aveW / counter;
		cd_prec density = Dens0;
		cd_prec mass = density / averageNrOfNeighbours / aveW;


		std::cout << "Rho/(dx*dy) = " << density / m_options.nrOfParticles[REAL] << "\n";
		std::cout <<"averageNrOfNeighbours =  " << averageNrOfNeighbours << ", New mass = Rho/n/aveW =  " << mass << "\n";
		if (computationalDomain_mode == MODE_CD::CD_DEBUG) { std::cout << "				Changing mass from " << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] / 2]->m_mass << " to " << mass << ".\n"; }

		for (auto*& i : SPH.Particles) { i->m_mass = mass; }
		//choose
		for (auto*& i : SPH.Particles) { i->m_mass = i->m_density.val*m_options.average_dim_steps.x*m_options.average_dim_steps.y*m_options.average_dim_steps.z; }
		//deletingParticlePairs();
		

	}
}


void SPH_CD::uniformGridPartitionInitialization() {
	if (computationalDomain_mode == MODE_CD::CD_DEBUG) { std::cout << "SPH_CD::uniformGridPartitionInitialization\n"; }
	UG = new uniformTreeGrid_2D(getXmin() - 2.0*m_options.smoothingKernelLengthCoefficient * InitialSmR, getXmax() + 2.0* m_options.smoothingKernelLengthCoefficient * InitialSmR, getYmin() - 2.0*m_options.smoothingKernelLengthCoefficient * InitialSmR, getYmax() + 2.0*m_options.smoothingKernelLengthCoefficient * InitialSmR);
	std::set<uniformTreeGrid_2D*> allpossibleNB;
	UG->FindNeighboursCells(allpossibleNB);
	UG->AssignNeighbourCells(allpossibleNB);
}



//  LOGIC CUALITY OF LIFE (SPH.ParticlePairs-loop)
// MAIN PARTICLE IS 
#define PPL_MAIN_PARTICLE_IS(type) (i->Mpart_ptr->m_type == PARTICLETYPE::type)
#define PPL_MAIN_PARTICLE_IS_NOT(type) (i->Mpart_ptr->m_type != PARTICLETYPE::type)
	// NEIGBOUR PARTICLE IS 
#define PPL_NEIGHBOUR_PARTICLE_IS(type) (i->NBpart_ptr->m_type == PARTICLETYPE::type)
#define PPL_NEIGHBOUR_PARTICLE_IS_NOT(type) (i->NBpart_ptr->m_type != PARTICLETYPE::type)


#define PPL_BOTH_PARTICLES_ARE(type) ((i->NBpart_ptr->m_type == PARTICLETYPE::type) and (i->Mpart_ptr->m_type == PARTICLETYPE::type))
#define PPL_BOTH_PARTICLES_ARE_NOT(type) ((i->NBpart_ptr->m_type != PARTICLETYPE::type) and (i->Mpart_ptr->m_type != PARTICLETYPE::type))

#define PPL_ONE_OF_PARTICLE_IS_1_OTHER_IS_2(type1,type2) (((i->NBpart_ptr->m_type == PARTICLETYPE::type1) and (i->Mpart_ptr->m_type == PARTICLETYPE::type2)) or ((i->NBpart_ptr->m_type == PARTICLETYPE::type2) and (i->Mpart_ptr->m_type == PARTICLETYPE::type1)))


void SPH_CD::creatingVirtualParticles() {
	if (computationalDomain_mode == MODE_CD::CD_DEBUG) {
    		std::cout << "SPH_CD::creatingVirtualParticles\n";
	}
	std::vector<glm::vec3> PB;
	PB.resize(m_options.nrOfParticles[PARTICLETYPE::REAL]+ m_options.nrOfParticles[PARTICLETYPE::BOUNDARY],glm::vec3(0.f));

	float Min_smR = 10;
	for (auto*& i : SPH.Particles) {
		if (i->m_SmR < Min_smR)
			Min_smR = i->m_SmR;
	}
	setStatMinScale(Min_smR);
	Min_smR *= m_options.smoothingKernelLengthCoefficient;
	int VirtualParticleCount = 0;

	for (auto*& i : SPH.ParticlePairs) {
		//ОПРЕДЕЛЕНИЕ ЧАСТИЦ, КОТОРЫЕ ДОЛЖНЫ БЫТЬ ОТЗЕРКАЛЕНЫ, ДЛЯ ЗАДАНИЯ ИМИ Г.У.
		//ОПРЕДЕЛИМ УРАВНЕНИЕ ПЛОСКОСТИ, В КОТОРОЙ ЛЕЖИТ ПОЛИГОН, ЗНАЯ ТОЧКУ, ПРИНАДЛЕЖАЩУЮ ЭТОЙ ПЛОСКОСТИ И ЕЁ НОРМАЛЬ
		//POS_BOUND - точка на плоскости;
		//NORMAL - нормаль плоскости,
		//ТОГДА УР-ИЕ ПЛОСКОСТИ ИМЕЕТ ВИД
		// NORMAL.X*(X-POS_BOUND.X) + NORMAL.Y*(Y-POS_BOUND.Y) + NORMAL.Z*(Z-POS_BOUND.Z) = 0
		// POS_REAL - точка рядом с плоскостью
		// КАНОНИЧЕСКОЕ УРАВНЕНИЕ ПРЯМОЙ, КОТОРАЯ ПРОХОДИТ ЧЕРЕЗ ТОЧКУ РЯДОМ С ПЛОСКОСТЬЮ И ИМЕЕТ НАПРАВЛЯЮЩИЙ ВЕКТОР - НОРМАЛЬ
		// (X-POS_REAL.X)/NORMAL.X = (Y-POS_REAL.Y)/NORMAL.Y = (Z-POS_REAL.Z)/NORMAL.Z = T
		// ИЗ РЕЗУЛЬТАТА ПЕРЕСЕЧЕНИЯ ЭТИХ УРАВНЕНИЙ МОЖНО НАЙТИ T:
		// T = (NORMAL.X*(POS_BOUND.X-NORMAL.X) + NORMAL.Y*(POS_BOUND.Y-NORMAL.Y) + NORMAL.Z*(POS_BOUND.Z-NORMAL.Z))/(NORMAL.X*POS_REAL.X + NORMAL.Z*POS_REAL.Z + NORMAL.Z*POS_REAL.Z)
		// А ЗНАЯ T, УЖЕ МОЖНО ОПРЕДЕЛИТЬ ПРОЕКЦИЮ ТОЧКИ НА ПЛОСКОСТЬ

		Particle* real_p = SPH.Particles[0];
		Particle* neighbour_p = SPH.Particles[0];

		if ((i->Mpart_ptr->m_type == PARTICLETYPE::REAL) and (i->NBpart_ptr->m_type == PARTICLETYPE::BOUNDARY)) {
			real_p = SPH.Particles[i->Mpart_ptr->m_id];
			neighbour_p = SPH.Particles[i->NBpart_ptr->m_id];
		}
		else if ((i->NBpart_ptr->m_type == PARTICLETYPE::REAL) and (i->Mpart_ptr->m_type == PARTICLETYPE::BOUNDARY)) {
			real_p = SPH.Particles[i->NBpart_ptr->m_id];
			neighbour_p = SPH.Particles[i->Mpart_ptr->m_id];
		}
		else {
			continue;
		}

		part_prec_3 normal = neighbour_p->InitPolygonNormal;
		part_prec_3 pos_bound = neighbour_p->m_position.val;
		part_prec_3 pos_real = real_p->m_position.val;
		part_prec T = ((normal.x*(pos_bound.x - pos_real.x)) + (normal.y*(pos_bound.y - pos_real.y)) + (normal.z*(pos_bound.z - pos_real.z))) / (normal.x*normal.x + normal.y*normal.y + normal.z*normal.z);
		part_prec_3 projected_real = glm::vec3(normal.x * T + pos_real.x, normal.y * T + pos_real.y, normal.z * T + pos_real.z);
		part_prec distance = abs(SPH.Particles[real_p->m_id]->distance(projected_real));
		//std::cout << "("<<pos_bound.x << ", " << pos_bound.y<<"), (" << pos_real.x << ", " << pos_real.y  << "), "<< distance << "\n";
		//SPH.Particles[real_p->m_id]->VirtualCounterpart = false;


		//std::cout << real_p->m_id << "    " << neighbour_p->m_id << "   " << normal.x << ",  " << normal.y << "\n";

		bool VCcondition = false;
		if (BM.boundarySecBar1()) {
			if (BM.boundarySecBar2(real_p->m_id)) {
				BM.BoundaryLinks[real_p->m_id].VirtualCounterpartReset();
				VCcondition = BM.BoundaryLinks[real_p->m_id].newVirtualCounterpart(normal, 0.9f);
			}
		}
		//std::cout << VCcondition << "\n";
		//if ((distance <= Max_smR) and !(SPH.Particles[real_p->m_id]->VirtualCounterpart)) {
		if ((distance <= 2.0*Min_smR) and !(VCcondition)) {

			part_prec_3 virtualVelocity;
			part_prec_3 virtualPosition;
			part_prec_3 distanceToPeriodic_part;

			if (SPH.Particles[neighbour_p->m_id]->particle_boundary->isPeriodic()) {
				glm::vec3 sum_MeshToLOOp_position(0.f, 0.f, 0.f);
				glm::vec3 sum_Mesh_position(0.f, 0.f, 0.f);
				for (int VI = 0;VI < SPH.Particles[neighbour_p->m_id]->particle_boundary->ReturnMeshToLoop()->getNrOfVertices();VI++) {
					sum_MeshToLOOp_position += SPH.Particles[neighbour_p->m_id]->particle_boundary->ReturnMeshToLoop()->getVertexArray()[VI].position;
					sum_Mesh_position += SPH.Particles[neighbour_p->m_id]->particle_boundary->ReturnMesh()->getVertexArray()[VI].position;
				}
				sum_MeshToLOOp_position = sum_MeshToLOOp_position / static_cast<float>(SPH.Particles[neighbour_p->m_id]->particle_boundary->ReturnMeshToLoop()->getNrOfVertices());
				sum_Mesh_position = sum_Mesh_position / static_cast<float>(SPH.Particles[neighbour_p->m_id]->particle_boundary->ReturnMesh()->getNrOfVertices());
				glm::vec3 distanceToPeriodic = sum_MeshToLOOp_position - sum_Mesh_position;
				distanceToPeriodic_part = { distanceToPeriodic.x,distanceToPeriodic.y,distanceToPeriodic.z };
				BM.TeleportationLinks.insert({ real_p->m_id,TeleportationData(distance, normal, distanceToPeriodic_part) });
				virtualVelocity = SPH.Particles[real_p->m_id]->m_velocity.val;
				virtualPosition = (projected_real + distanceToPeriodic_part + normal * distance);
			}
			else {
				virtualVelocity = static_cast<part_prec>(2.0) * SPH.Particles[neighbour_p->m_id]->m_velocity.val - SPH.Particles[real_p->m_id]->m_velocity.val;
				virtualPosition = (projected_real - normal * distance);
			}

			if(virtualPosition != pos_real){
				//std::cout << "virtualVelocity{x,y,z} = {" << virtualVelocity.x << ", " << virtualVelocity.y << ", " << virtualVelocity.z << "}\n";
				SPH.Particles.push_back(new Particle(m_options.nrOfParticles[PARTICLETYPE::REAL] + m_options.nrOfParticles[PARTICLETYPE::BOUNDARY] + VirtualParticleCount, PARTICLETYPE::VIRTUAL, virtualPosition, virtualVelocity, SPH.Particles[i->Mpart_ptr->m_id]->m_SmR, SPH.Particles[i->Mpart_ptr->m_id]->m_density.val, SPH.Particles[i->Mpart_ptr->m_id]->m_mass));
				VirtualParticleCount++;

				if (BM.boundarySecBar2(real_p->m_id)) {
					if (BM.BoundaryLinks[real_p->m_id].NrOfNormals() >= 1) {
						BM.BoundaryLinks[real_p->m_id].addingBoundary(normal, distance, SPH.Particles[neighbour_p->m_id]->m_velocity.val, SPH.Particles[neighbour_p->m_id]->particle_boundary->isPeriodic(), distanceToPeriodic_part, virtualPosition, virtualVelocity);
					}
				}
				else {
					BM.BoundaryLinks.insert({ real_p->m_id,BoundaryData(normal,distance,SPH.Particles[neighbour_p->m_id]->m_velocity.val) });
					BM.BoundaryLinks[real_p->m_id].furtherInitialization(SPH.Particles[neighbour_p->m_id]->particle_boundary->isPeriodic(), distanceToPeriodic_part, virtualPosition, virtualVelocity);
				}
				if (m_options.cornerVP) {

					if (BM.BoundaryLinks[real_p->m_id].NrOfNormals() >= 2) {
					//if (SPH.Particles[real_p->m_id]->VirtualCounterpartNormals.size() == 2) {

						
						// SPH.Particles[real_p->m_id]->periodicBoundary[0]
						if ((BM.BoundaryLinks[real_p->m_id].firstBoundaryisPeriodic() == false) and (SPH.Particles[neighbour_p->m_id]->particle_boundary->isPeriodic() == false)) {
							//virtualVelocity = 0.5*(SPH.Particles[neighbour_p->m_id]->m_velocity.val + SPH.Particles[real_p->m_id]->VC_BoundaryVelocity[0]) - (SPH.Particles[real_p->m_id]->m_velocity.val - 0.5*(SPH.Particles[neighbour_p->m_id]->m_velocity.val + SPH.Particles[real_p->m_id]->VC_BoundaryVelocity[0]));
							//virtualVelocity = 0.5*(SPH.Particles[neighbour_p->m_id]->m_velocity.val + BM.BoundaryLinks[real_p->m_id].firstBoundaryVelocity()) - (SPH.Particles[real_p->m_id]->m_velocity.val - 0.5*(SPH.Particles[neighbour_p->m_id]->m_velocity.val;
							
							//virtualVelocity = 0.5*(SPH.Particles[neighbour_p->m_id]->m_velocity.val + BM.BoundaryLinks[real_p->m_id].firstBoundaryVelocity()) - (SPH.Particles[real_p->m_id]->m_velocity.val - 0.5*(SPH.Particles[neighbour_p->m_id]->m_velocity.val + BM.BoundaryLinks[real_p->m_id].firstBoundaryVelocity()));
							virtualVelocity = -virtualVelocity;
							//virtualPosition = (projected_real - normal * distance - 2.0 * SPH.Particles[real_p->m_id]->VirtualCounterpartNormals[0] * SPH.Particles[real_p->m_id]->VC_DistanceToBoundary[0]);
							virtualPosition = (projected_real - normal * distance - static_cast<part_prec>(2.0) * BM.BoundaryLinks[real_p->m_id].firstDistanceVector());
							//std::cout << "From not periodic to not periodic\n";
							}
						else if ((BM.BoundaryLinks[real_p->m_id].firstBoundaryisPeriodic() == true) and (SPH.Particles[neighbour_p->m_id]->particle_boundary->isPeriodic() == true)) {
								glm::vec3 sum_MeshToLOOp_position(0.f, 0.f, 0.f);
								glm::vec3 sum_Mesh_position(0.f, 0.f, 0.f);
								for (int VI = 0;VI < SPH.Particles[neighbour_p->m_id]->particle_boundary->ReturnMeshToLoop()->getNrOfVertices();VI++) {
									sum_MeshToLOOp_position += SPH.Particles[neighbour_p->m_id]->particle_boundary->ReturnMeshToLoop()->getVertexArray()[VI].position;
									sum_Mesh_position += SPH.Particles[neighbour_p->m_id]->particle_boundary->ReturnMesh()->getVertexArray()[VI].position;
								}
								sum_MeshToLOOp_position = sum_MeshToLOOp_position / static_cast<float>(SPH.Particles[neighbour_p->m_id]->particle_boundary->ReturnMeshToLoop()->getNrOfVertices());
								sum_Mesh_position = sum_Mesh_position / static_cast<float>(SPH.Particles[neighbour_p->m_id]->particle_boundary->ReturnMesh()->getNrOfVertices());

								glm::vec3 distanceToPeriodic = sum_MeshToLOOp_position - sum_Mesh_position;
								part_prec_3 distanceToPeriodic_part = { distanceToPeriodic.x,distanceToPeriodic.y,distanceToPeriodic.z };

								//virtualPosition += SPH.Particles[real_p->m_id]->VirtualCounterpartNormals[0] * SPH.Particles[real_p->m_id]->VC_DistanceToBoundary[0];
								virtualPosition += BM.BoundaryLinks[real_p->m_id].firstDistanceVector();

								//std::cout << "From periodic to periodic\n";
							}
						else if ((BM.BoundaryLinks[real_p->m_id].firstBoundaryisPeriodic() == true) and (SPH.Particles[neighbour_p->m_id]->particle_boundary->isPeriodic() == false)) {
								//virtualPosition += SPH.Particles[real_p->m_id]->VC_DistanceToPeriodic[0];
								virtualPosition += BM.BoundaryLinks[real_p->m_id].firstDistanceToPeriodic();
								//std::cout << "From periodic to not periodic\n";
							}
						else if ((BM.BoundaryLinks[real_p->m_id].firstBoundaryisPeriodic() == false) and (SPH.Particles[neighbour_p->m_id]->particle_boundary->isPeriodic() == true)) {
							glm::vec3 sum_MeshToLOOp_position(0.f, 0.f, 0.f);
							glm::vec3 sum_Mesh_position(0.f, 0.f, 0.f);
							for (int VI = 0;VI < SPH.Particles[neighbour_p->m_id]->particle_boundary->ReturnMeshToLoop()->getNrOfVertices();VI++) {
								sum_MeshToLOOp_position += SPH.Particles[neighbour_p->m_id]->particle_boundary->ReturnMeshToLoop()->getVertexArray()[VI].position;
								sum_Mesh_position += SPH.Particles[neighbour_p->m_id]->particle_boundary->ReturnMesh()->getVertexArray()[VI].position;
							}
							sum_MeshToLOOp_position = sum_MeshToLOOp_position / static_cast<float>(SPH.Particles[neighbour_p->m_id]->particle_boundary->ReturnMeshToLoop()->getNrOfVertices());
							sum_Mesh_position = sum_Mesh_position / static_cast<float>(SPH.Particles[neighbour_p->m_id]->particle_boundary->ReturnMesh()->getNrOfVertices());

							glm::vec3 distanceToPeriodic = sum_MeshToLOOp_position - sum_Mesh_position;
							distanceToPeriodic_part = { distanceToPeriodic.x,distanceToPeriodic.y,distanceToPeriodic.z };
							BM.TeleportationLinks[real_p->m_id] = TeleportationData(distance, normal, distanceToPeriodic_part);
							//virtualVelocity = SPH.Particles[real_p->m_id]->VC_Before_PeriodicVelocity[0];
							virtualVelocity = BM.BoundaryLinks[real_p->m_id].firstPeriodicBoundaryVelocity();
							//virtualPosition -= 2.0* SPH.Particles[real_p->m_id]->VirtualCounterpartNormals[0] * SPH.Particles[real_p->m_id]->VC_DistanceToBoundary[0];
							virtualPosition -= static_cast<part_prec>(2.0)* BM.BoundaryLinks[real_p->m_id].firstDistanceVector();
							//std::cout << "From not periodic to periodic\n";
						}
						SPH.Particles.push_back(new Particle(m_options.nrOfParticles[PARTICLETYPE::REAL] + m_options.nrOfParticles[PARTICLETYPE::BOUNDARY] + VirtualParticleCount, PARTICLETYPE::VIRTUAL, virtualPosition, virtualVelocity, SPH.Particles[i->Mpart_ptr->m_id]->m_SmR, SPH.Particles[i->Mpart_ptr->m_id]->m_density.val, SPH.Particles[i->Mpart_ptr->m_id]->m_mass));
						//if((virtualPosition.x > 0.5) and (virtualPosition.y > 0.5))
						//	std::cout << m_options.nrOfParticles[PARTICLETYPE::REAL] + m_options.nrOfParticles[PARTICLETYPE::BOUNDARY] + VirtualParticleCount << " ";
						VirtualParticleCount++;
					}
				}
				part_prec_3 Xij = { SPH.Particles[real_p->m_id]->dx(projected_real),SPH.Particles[real_p->m_id]->dy(projected_real),SPH.Particles[real_p->m_id]->dz(projected_real) };
				for (int d = 0;d < nrOfDim;d++) {
					if (BH_r/ distance <= 1) {
						PB[real_p->m_id][d] = BH_D * (pow(BH_r / distance, 12) - pow(BH_r / distance, 4))*Xij[d] / pow(distance, 2);
					}
					else {
						PB[real_p->m_id][d] = 0.f;
					}
				}
			}
		}
	}
	//std::cout << "\n COUNT: " << count << "\n";
	for (int i = 0;i < PB.size();i++) {
		//SPH.Particles[i]->m_velocity.dval += PB[i];
	}

	m_options.nrOfParticles[PARTICLETYPE::VIRTUAL] = VirtualParticleCount;
	if (computationalDomain_mode == MODE_CD::CD_DEBUG) {
		std::cout << "SPH_CD::CreatingVirtualParticles::VIRTUAL_PARTICLES_WERE_CREATED::" << VirtualParticleCount << "\n";
	}
}


void SPH_CD::addingVirtualTOPairs() {
	if (computationalDomain_mode == MODE_CD::CD_DEBUG) {
		std::cout << "SPH_CD::addingVirtualTOPairs\n";
	}
	int novirtpairs = SPH.ParticlePairs.size();
	switch (m_options.NBSAlg) {
	case(NEIGBOURS_SEARCH_ALGORITHM::DIRECT):
		// Пытаемся избежать лишних пар
		for (int i = 0;i < SPH.Particles.size();i++) {  //  SPH.Particles.size()-1 ет.к. у последнего уже нету соседей
			Particle* const i_p = SPH.Particles[i];
			//std::cout <<"\n"<< i_p->m_id << "_th particle has neighbours: ";
			//std::cout << m_options.nrOfParticles[PARTICLETYPE::REAL] + m_options.nrOfParticles[PARTICLETYPE::BOUNDARY];
			for (int j = i + 1;j < SPH.Particles.size();j++) { //  SPH.Particles.size()-1 ет.к. у последнего уже нету соседей
			//for (int j = m_options.nrOfParticles[PARTICLETYPE::REAL] + m_options.nrOfParticles[PARTICLETYPE::BOUNDARY] + i;j < SPH.Particles.size();j++) { //  SPH.Particles.size()-1 ет.к. у последнего уже нету соседей
				Particle* const j_p = SPH.Particles[j];
				if (i_p == j_p)
					continue;
				if ((j_p->m_type == PARTICLETYPE::BOUNDARY) and (i_p->m_type == PARTICLETYPE::BOUNDARY)) // Не учитываем пары Граница-Граница  
					continue;
				if ((j_p->m_type == PARTICLETYPE::VIRTUAL) or (i_p->m_type == PARTICLETYPE::VIRTUAL)) {
					if ((j_p->m_type == PARTICLETYPE::VIRTUAL) and (i_p->m_type == PARTICLETYPE::VIRTUAL))
						continue;
					if (i_p->distance(j_p) <= m_options.smoothingKernelLengthCoefficient * i_p->m_SmR) {
						//std::cout << i_p->m_id << "    " << j_p->m_id << "\n";
						SPH.ParticlePairs.push_back(new ParticlePair(i_p, j_p, nrOfDim));
					}
				}
				//j_p = nullptr;
			}
			//i_p = nullptr;
		}
		break;
	case(NEIGBOURS_SEARCH_ALGORITHM::UNIFORM_GRID):
		UG->clear();
		for (int i = 0;i < SPH.Particles.size();i++) {  //  SPH.Particles.size()-1 ет.к. у последнего уже нету соседей
			Particle* i_p = SPH.Particles[i];
			UG->addPoint(i_p);
		}
 		UG->addingVirtualParticles(SPH.ParticlePairs, m_options.smoothingKernelLengthCoefficient, nrOfDim);
		break;
	default:
		break;
	}



	if (computationalDomain_mode == MODE_CD::CD_DEBUG) {
		std::cout << "SPH_CD::addingVirtualTOPairs::PARTICLE_PAIRS_WERE_CREATED::" << novirtpairs << " + " << SPH.ParticlePairs.size() - novirtpairs << " = " << SPH.ParticlePairs.size() << "\n";
	}

}


void SPH_CD::deletingParticlePairs() {
	if (computationalDomain_mode == MODE_CD::CD_DEBUG) {
		std::cout << "SPH_CD::deletingParticlePairs   ";
	}
	if (SPH.ParticlePairs.size() != 0) {
		if (computationalDomain_mode == MODE_CD::CD_DEBUG) {
			std::cout << SPH.ParticlePairs.size() << "->";
		}
		for (auto*& i : SPH.ParticlePairs) {
			safeDelete(&i);
			if (m_options.NBSAlg == NEIGBOURS_SEARCH_ALGORITHM::UNIFORM_GRID) {
				UG->clear();
			}
		}
		SPH.ParticlePairs.resize(0);
		if (computationalDomain_mode == MODE_CD::CD_DEBUG) {
			std::cout << SPH.ParticlePairs.size() << "\n";
		}
	}
	else {
		if (computationalDomain_mode == MODE_CD::CD_DEBUG) {
			std::cout << "\n";
		}

	}
}

void SPH_CD::DensityVariation() {
	BH_r = 0.f;
	int BH_r_count = 0;

	float minDistance = (getXmax() - getXmin());

	float dvxdWx;
	std::vector<float> drhodt;
	drhodt.resize(SPH.Particles.size(),0.f);
	std::vector<ParticlePair::dsmthFunc> Funckey{ &Particle::dVx ,&Particle::dVy ,&Particle::dVz };
	for (auto*& i : SPH.ParticlePairs) {
		for (int d = 0;d < nrOfDim;d++) {
			//drhodt[i->Mpart_ptr->m_id] += (i->Mpart_ptr->m_density.val)*((i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val) * dvxdWx);
			//drhodt[i->NBpart_ptr->m_id] += (i->NBpart_ptr->m_density.val)*((i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val) * dvxdWx);

			if ((PPL_BOTH_PARTICLES_ARE(REAL)) or
				(PPL_MAIN_PARTICLE_IS(VIRTUAL) or PPL_NEIGHBOUR_PARTICLE_IS(VIRTUAL))) {
				switch (m_options.densityChangeUsing) {
				case(VAL):
					drhodt[i->Mpart_ptr->m_id] += (i->Mpart_ptr->m_mass) * i->W(i->Mpart_ptr);
					drhodt[i->NBpart_ptr->m_id] += (i->NBpart_ptr->m_mass) * i->W(i->NBpart_ptr);
					break;
				case(DVAL_DT):
					switch (m_options.densityChangeAlgorithm) {
					case(0): //мой вывод (как я понимаю)
						drhodt[i->Mpart_ptr->m_id] += (i->NBpart_ptr->m_mass) * (-i->dvfunc(Funckey[d])) * (i->dWd(static_cast<DiffAxis>(d), i->Mpart_ptr));
						drhodt[i->NBpart_ptr->m_id] += (i->Mpart_ptr->m_mass) * (i->dvfunc(Funckey[d])) * (i->dWd(static_cast<DiffAxis>(d), i->NBpart_ptr));
						break;
					case(1):
						drhodt[i->Mpart_ptr->m_id] += -(i->Mpart_ptr->m_density.val)*(i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val) *(-i->dvfunc(Funckey[d]))*(i->dWd(static_cast<DiffAxis>(d), i->Mpart_ptr));
						drhodt[i->NBpart_ptr->m_id] += -(i->NBpart_ptr->m_density.val)*(i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val) *(i->dvfunc(Funckey[d]))*(i->dWd(static_cast<DiffAxis>(d), i->NBpart_ptr));
						break;
					case(2): //Monaghan
						drhodt[i->Mpart_ptr->m_id] += (i->NBpart_ptr->m_mass)*(i->dvfunc(Funckey[d]))*(i->dWd(static_cast<DiffAxis>(d), i->Mpart_ptr));
						drhodt[i->NBpart_ptr->m_id] += (i->Mpart_ptr->m_mass)*(-i->dvfunc(Funckey[d]))*(i->dWd(static_cast<DiffAxis>(d), i->NBpart_ptr));
						break;
					default:
						break;
					}
					break;
				}
			}
		}
		if((i->Mpart_ptr->m_type == REAL)or(i->NBpart_ptr->m_type == REAL)){
			if (i->Mpart_ptr->distance(i->NBpart_ptr) < minDistance) {
				minDistance = i->Mpart_ptr->distance(i->NBpart_ptr);
			}

			BH_r += i->Mpart_ptr->distance(i->NBpart_ptr);
			BH_r_count++;
		}
	}
	BH_r = BH_r / BH_r_count;
	BH_r = 0.05f/2;

	for (auto*& i : SPH.Particles) {
		if ((i->m_type == PARTICLETYPE::REAL) or (i->m_type == PARTICLETYPE::BOUNDARY))

			i->m_density.dval = drhodt[i->m_id];

	}
	//std::cout << SPH.Particles[200]->m_density.dval << " ";

	setStatAveScale(BH_r);
}

int delta(int m, int n) {
	if (m == n)
		return 1;
	else
		return 0;
}


/*
void SPH_CD::EvaluateEps() {
	std::vector<glm::mat3x3> eps;
	eps.resize(SPH.Particles.size());
	float mrho;
	std::vector<ParticlePair::dsmthFunc> Funckey{ &Particle::dVx ,&Particle::dVy ,&Particle::dVz };
	for (auto*& i : SPH.ParticlePairs) {
		mrho = (i->NBpart_ptr)->m_mass / (i->NBpart_ptr)->m_density.val;
		for (int d = 0;d < nrOfDim;d++) {
			for (int dd = 0;dd < nrOfDim;dd++) {
				eps[i->Mpart_ptr->m_id][d][dd] += (mrho * (i->dvfunc(Funckey[d])) * (i->m_dWdp[d])) + (mrho * (i->dvfunc(Funckey[dd])) * (i->m_dWdp[dd])) - 2 / 3 * mrho*(i->dvfunc(Funckey[0]) + i->dvfunc(Funckey[1]) + i->dvfunc(Funckey[2]) *(i->m_dWdp[0] + i->m_dWdp[1] + i->m_dWdp[2])*delta(d, dd));
				eps[i->NBpart_ptr->m_id][d][dd] += (mrho * (-i->dvfunc(Funckey[d])) * (i->m_dWdp[d])) + (mrho * (-i->dvfunc(Funckey[dd])) * (i->m_dWdp[dd])) - 2 / 3 * mrho*-(i->dvfunc(Funckey[0]) + i->dvfunc(Funckey[1]) + i->dvfunc(Funckey[2]) *(i->m_dWdp[0] + i->m_dWdp[1] + i->m_dWdp[2])*delta(d, dd));
			}
		}
	}
	for (auto*& i : SPH.Particles) {
		if (i->m_type == PARTICLETYPE::REAL)
			i->m_eps = eps[i->m_id];
	}
}

void SPH_CD::InternalForces() {

	std::vector<glm::vec3> dpdx;
	dpdx.resize(SPH.Particles.size());
	std::vector<glm::vec3> dtaudx;
	dtaudx.resize(SPH.Particles.size());
	std::vector<glm::vec3> TAU;
	TAU.resize(SPH.Particles.size());
	std::vector<glm::vec3> dvdt;
	dvdt.resize(SPH.Particles.size());

	float presspart;

	for (auto*& i : SPH.ParticlePairs) {
		presspart = (i->NBpart_ptr->m_mass) *((i->Mpart_ptr->m_pressure.val)+(i->NBpart_ptr->m_pressure.val))/ ((i->Mpart_ptr->m_density.val)*(i->NBpart_ptr->m_density.val));
		for (int d = 0;d < nrOfDim;d++) {
			dpdx[i->Mpart_ptr->m_id][d] += (-presspart)*(i->m_dWdp[d]);
			dpdx[i->NBpart_ptr->m_id][d] += (-presspart)*(-i->m_dWdp[d]);

			dtaudx[i->Mpart_ptr->m_id][d] += (i->NBpart_ptr->m_mass) / ((i->Mpart_ptr->m_density.val)*(i->NBpart_ptr->m_density.val))*i->m_dWdp[d];
			dtaudx[i->NBpart_ptr->m_id][d] += (i->Mpart_ptr->m_mass) / ((i->NBpart_ptr->m_density.val)*(i->Mpart_ptr->m_density.val))*i->m_dWdp[d];

			for (int dd = 0;dd < nrOfDim;dd++) {
				TAU[i->Mpart_ptr->m_id][d] += ((i->NBpart_ptr->m_DVisc)*(i->NBpart_ptr->m_eps[d][dd]) + (i->Mpart_ptr->m_DVisc)*(i->Mpart_ptr->m_eps[d][dd]));
			}
			dtaudx[i->Mpart_ptr->m_id][d] *= TAU[i->Mpart_ptr->m_id][d];

			dvdt[i->Mpart_ptr->m_id][d] += dpdx[i->Mpart_ptr->m_id][d] + dtaudx[i->Mpart_ptr->m_id][d];
		}
	}
	for (auto*& i : SPH.Particles) {
		if(i->m_type == PARTICLETYPE::REAL)
			i->m_velocity.dval = dvdt[i->m_id];
	}
}
*/


void SPH_CD::ExternalForces() {
	if (computationalDomain_mode == MODE_CD::CD_DEBUG) {
		std::cout << "SPH_CD::ExternalForces\n";
	}
	if (true) { //gravityFlag
		part_prec_3 gravVec = { 0.0,-9.81, 0.0 };
		gravVec = { 0.0f, 0.0f, 0.0f };
		for (auto*& i : SPH.Particles) {
			i->m_velocity.dval += gravVec;
		}
	}

}

void SPH_CD::SaveMaxVelocity() {

	if (computationalDomain_mode == MODE_CD::CD_DEBUG) {
		std::cout << "SPH_CD::SaveMaxVelocity\n";
	}
	BH_D = 0.f;
	part_prec max_drho = 0.0;
	for (auto*& i : SPH.Particles) {
		i->absVelocity = sqrt(pow(i->m_velocity.val.x, 2)+ pow(i->m_velocity.val.y, 2)+ pow(i->m_velocity.val.z, 2));
		//Boundary handling parametrs assignment
		if(i->m_type == REAL){
			if (BH_D < i->absVelocity) {
				BH_D = i->absVelocity;
			}
			if (i->m_density.dval * getDeltaTime() > max_drho) {
				max_drho = i->m_density.dval * getDeltaTime();
			}
		}
	}
	

	setStatMaxVel(10.0*glm::max(BH_D,sqrt(max_drho/ Dens0)));
}

void SPH_CD::Boundary_Handling() {
	std::vector<glm::vec3> PB;
	PB.resize(SPH.Particles.size());

	std::vector<ParticlePair::dsmthFunc> Funckey{ &Particle::dx ,&Particle::dy ,&Particle::dz };
	//std::cout << BH_D << "   " << BH_r<<"\n";
	for (auto*& i : SPH.ParticlePairs) {
		//if (i->NBpart_ptr->m_type == VIRTUAL) {
			for (int d = 0;d < nrOfDim;d++) {
				//std::cout << BH_r / i->Mpart_ptr->distance(i->NBpart_ptr) << "\n";
				if (BH_r / i->Mpart_ptr->distance(i->NBpart_ptr) <= 1.f) {
					PB[i->Mpart_ptr->m_id][d] += BH_D * (pow(BH_r / i->Mpart_ptr->distance(i->NBpart_ptr), 12) - pow(BH_r / i->Mpart_ptr->distance(i->NBpart_ptr), 4))*i->dvfunc(Funckey[d]) / pow(i->Mpart_ptr->distance(i->NBpart_ptr), 2);
					std::cout << BH_D << "*(" << pow(BH_r / i->Mpart_ptr->distance(i->NBpart_ptr), 12) << " - " << pow(BH_r / i->Mpart_ptr->distance(i->NBpart_ptr), 4) << ")*" << i->dvfunc(Funckey[d])<<"/" << pow(i->Mpart_ptr->distance(i->NBpart_ptr), 2) <<"\n";
					std::cout << BH_r / i->Mpart_ptr->distance(i->NBpart_ptr) << "    " <<i->Mpart_ptr->m_id << ": " << PB[i->Mpart_ptr->m_id][0] << ", " << PB[i->Mpart_ptr->m_id][1] << ", " << PB[i->Mpart_ptr->m_id][2] << "\n";
				}
			}
		//}
	}
	for (auto*& i : SPH.Particles) {
		i->m_velocity.val += PB[i->m_id];
	}
}


//void SPH_CD::CorrectiveKernelCalculation() {
//	if (computationalDomain_mode == MODE_CD::CD_DEBUG) {
//		std::cout << "SPH_CD::CorrectiveKernelCalculation\n";
//	}
//
//	std::vector<float> corKer;
//	corKer.resize(SPH.Particles.size(),0.f);
//
//	for (auto*& i : SPH.ParticlePairs) {
//		if (i->Mpart_ptr->m_type == PARTICLETYPE::REAL and i->NBpart_ptr->m_type == PARTICLETYPE::BOUNDARY) {
//			corKer[i->Mpart_ptr->m_id] += i->NBpart_ptr->m_mass / i->NBpart_ptr->m_density.val * i->W(i->Mpart_ptr);
//		}
//		//std::cout << i->m_W << "\n";
//	}
//	for (auto*& i : SPH.Particles) {
//		if(i->m_type == PARTICLETYPE::REAL) {
//			i->correctiveKernel = corKer[i->m_id];
//		}
//		//std::cout << i->m_id << ": " << i->correctiveKernel << "\n";
//	}
//}


SPH_CD::~SPH_CD() {

	for (auto*& i : SPH.Particles) {
		delete i;
	}
	for (auto*& i : SPH.ParticlePairs) {
		delete i;
	}
	delete UG;
}

//#define COLORING_PARAMETR_NAME "abs(V)"
//#define COLORING_INPUT_PARAMETR(i) sqrt(pow(i->m_velocity.val.x,2)+pow(i->m_velocity.val.y,2)+pow(i->m_velocity.val.z,2))

#define COLORING_PARAMETR_NAME "Vx"
#define COLORING_INPUT_PARAMETR(i) i->m_velocity.val.x

//#define COLORING_PARAMETR_NAME "Vy"
//#define COLORING_INPUT_PARAMETR(i) i->m_velocity.val.y

//#define COLORING_PARAMETR_NAME "Smoothing radius"
//#define COLORING_INPUT_PARAMETR(i) i->m_SmR

//#define COLORING_PARAMETR_NAME "d(Vx)/dt "
//#define COLORING_INPUT_PARAMETR(i) i->m_velocity.dval.x

//#define COLORING_PARAMETR_NAME "d(Vy)/dt "
//#define COLORING_INPUT_PARAMETR(i) i->m_velocity.dval.y

//#define COLORING_PARAMETR_NAME "Density"
//#define COLORING_INPUT_PARAMETR(i) i->m_density.val

//#define COLORING_PARAMETR_NAME "d(Density)/dt"
//#define COLORING_INPUT_PARAMETR(i) i->m_density.dval

//#define COLORING_PARAMETR_NAME "pressure"
//#define COLORING_INPUT_PARAMETR(i) i->m_pressure.val

//#define COLORING_PARAMETR_NAME "id"
//#define COLORING_INPUT_PARAMETR(i) i->m_id

//#define COLORING_PARAMETR_NAME "nrOfNeighbours"
//#define COLORING_INPUT_PARAMETR(i) i->nrOfNeighbours

//#define COLORING_PARAMETR_NAME "mass"
//#define COLORING_INPUT_PARAMETR(i) i->m_mass




//#define COLORING_CONDITION(i) i->m_type == REAL
//#define COLORING_CONDITION(i) i->m_type == VIRTUAL
#define COLORING_CONDITION(i) (i->m_type == REAL)or(i->m_type == BOUNDARY)or(i->m_type == VIRTUAL)
//#define COLORING_CONDITION(i) (i->m_type == REAL)or(i->m_type == VIRTUAL)
//#define COLORING_CONDITION(i) (i->m_type == REAL)or(i->m_type == BOUNDARY)
void SPH_CD::Coloring() {
	if (computationalDomain_mode == MODE_CD::CD_DEBUG) {
		std::cout << "SPH_CD::Coloring\n";
	}

	float Red = 0.f;
	float Green = 0.f;
	float Blue = 0.f;
	float currentVal;

	//glm::vec3 maxcolor{ 1.f,0.0f,0.0f };
	//glm::vec3 avecolor{ 0.f,1.0f,0.0f };
	//glm::vec3 mincolor{ 0.f,0.0f,1.0f };



//#define COLORING_PARAMETR_NAME "id"
//#define COLORING_INPUT_PARAMETR(i) i->m_id

//#define COLORING_PARAMETR_NAME "nrOfNeighbours"
//#define COLORING_INPUT_PARAMETR(i) i->nrOfNeighbours

//#define COLORING_PARAMETR_NAME "mass"
//#define COLORING_INPUT_PARAMETR(i) i->m_mass

	ColoringUnknow = false;
	if (getColorParam() == "Vx") {
		changeColorParamTo("Vx");
		float maxVal = SPH.Particles[0]->m_velocity.val.x;
		float minVal = SPH.Particles[0]->m_velocity.val.x;
		for (auto*& i : SPH.Particles) {
			if (COLORING_CONDITION(i)) {
				if (i->m_velocity.val.x < minVal) { minVal = i->m_velocity.val.x; }
				if (i->m_velocity.val.x > maxVal) { maxVal = i->m_velocity.val.x; }
			}
		}
		float aveVal = (minVal + maxVal) / 2;
		m_maxVal = maxVal;
		m_minVal = minVal;

		for (auto*& i : SPH.Particles) {

			if (COLORING_CONDITION(i)) {
				currentVal = i->m_velocity.val.x;

				if ((currentVal >= minVal) and (currentVal <= aveVal)) {
					Red = 0.f;
					Green = (1.f - 0.f) / (aveVal - minVal)*(currentVal - minVal) + 0.f;
					Blue = 1.f - (1.f - 0.f) / (aveVal - minVal)*(currentVal - minVal);
				}
				else if ((currentVal > aveVal) and (currentVal <= maxVal)) {
					Red = (1.f - 0.f) / (maxVal - aveVal)*(currentVal - aveVal) + 0.f;
					Green = 1.f - (1.f - 0.f) / (maxVal - aveVal)*(currentVal - aveVal);
					Blue = 0.f;
				}
				if (minVal == maxVal) {
					Red = 0.f;
					Green = 1.f;
					Blue = 0.f;
				}
				i->m_color = glm::vec3(Red, Green, Blue);
			}
		}
	}
	else if (getColorParam() == "Vy") {
		changeColorParamTo("Vy");
		float maxVal = SPH.Particles[0]->m_velocity.val.y;
		float minVal = SPH.Particles[0]->m_velocity.val.y;
		for (auto*& i : SPH.Particles) {
			if (COLORING_CONDITION(i)) {
				if (i->m_velocity.val.y < minVal) { minVal = i->m_velocity.val.y; }
				if (i->m_velocity.val.y > maxVal) { maxVal = i->m_velocity.val.y; }
			}
		}
		float aveVal = (minVal + maxVal) / 2;
		m_maxVal = maxVal;
		m_minVal = minVal;

		for (auto*& i : SPH.Particles) {

			if (COLORING_CONDITION(i)) {
				currentVal = i->m_velocity.val.y;

				if ((currentVal >= minVal) and (currentVal <= aveVal)) {
					Red = 0.f;
					Green = (1.f - 0.f) / (aveVal - minVal)*(currentVal - minVal) + 0.f;
					Blue = 1.f - (1.f - 0.f) / (aveVal - minVal)*(currentVal - minVal);
				}
				else if ((currentVal > aveVal) and (currentVal <= maxVal)) {
					Red = (1.f - 0.f) / (maxVal - aveVal)*(currentVal - aveVal) + 0.f;
					Green = 1.f - (1.f - 0.f) / (maxVal - aveVal)*(currentVal - aveVal);
					Blue = 0.f;
				}
				if (minVal == maxVal) {
					Red = 0.f;
					Green = 1.f;
					Blue = 0.f;
				}
				i->m_color = glm::vec3(Red, Green, Blue);
			}
		}
	}
	else if (getColorParam() == "abs(V)") {
		changeColorParamTo("abs(V)");
		float maxVal = sqrt(pow(SPH.Particles[0]->m_velocity.val.x, 2) + pow(SPH.Particles[0]->m_velocity.val.y, 2) + pow(SPH.Particles[0]->m_velocity.val.z, 2));
		float minVal = sqrt(pow(SPH.Particles[0]->m_velocity.val.x, 2) + pow(SPH.Particles[0]->m_velocity.val.y, 2) + pow(SPH.Particles[0]->m_velocity.val.z, 2));
		for (auto*& i : SPH.Particles) {
			if (COLORING_CONDITION(i)) {
				if (sqrt(pow(i->m_velocity.val.x, 2) + pow(i->m_velocity.val.y, 2) + pow(i->m_velocity.val.z, 2)) < minVal) { minVal = sqrt(pow(i->m_velocity.val.x, 2) + pow(i->m_velocity.val.y, 2) + pow(i->m_velocity.val.z, 2)); }
				if (sqrt(pow(i->m_velocity.val.x, 2) + pow(i->m_velocity.val.y, 2) + pow(i->m_velocity.val.z, 2)) > maxVal) { maxVal = sqrt(pow(i->m_velocity.val.x, 2) + pow(i->m_velocity.val.y, 2) + pow(i->m_velocity.val.z, 2)); }
			}
		}
		float aveVal = (minVal + maxVal) / 2;
		m_maxVal = maxVal;
		m_minVal = minVal;

		for (auto*& i : SPH.Particles) {

			if (COLORING_CONDITION(i)) {
				currentVal = sqrt(pow(i->m_velocity.val.x, 2) + pow(i->m_velocity.val.y, 2) + pow(i->m_velocity.val.z, 2));

				if ((currentVal >= minVal) and (currentVal <= aveVal)) {
					Red = 0.f;
					Green = (1.f - 0.f) / (aveVal - minVal)*(currentVal - minVal) + 0.f;
					Blue = 1.f - (1.f - 0.f) / (aveVal - minVal)*(currentVal - minVal);
				}
				else if ((currentVal > aveVal) and (currentVal <= maxVal)) {
					Red = (1.f - 0.f) / (maxVal - aveVal)*(currentVal - aveVal) + 0.f;
					Green = 1.f - (1.f - 0.f) / (maxVal - aveVal)*(currentVal - aveVal);
					Blue = 0.f;
				}
				if (minVal == maxVal) {
					Red = 0.f;
					Green = 1.f;
					Blue = 0.f;
				}
				i->m_color = glm::vec3(Red, Green, Blue);
			}
		}
	}
	else if (getColorParam() == "Smoothing radius") {
	changeColorParamTo("Smoothing radius");
	float maxVal = SPH.Particles[0]->m_SmR;
	float minVal = SPH.Particles[0]->m_SmR;
	for (auto*& i : SPH.Particles) {
		if (COLORING_CONDITION(i)) {
			if (i->m_SmR < minVal) { minVal = i->m_SmR; }
			if (i->m_SmR > maxVal) { maxVal = i->m_SmR; }
		}
	}
	float aveVal = (minVal + maxVal) / 2;
	m_maxVal = maxVal;
	m_minVal = minVal;

	for (auto*& i : SPH.Particles) {

		if (COLORING_CONDITION(i)) {
			currentVal = i->m_SmR;

			if ((currentVal >= minVal) and (currentVal <= aveVal)) {
				Red = 0.f;
				Green = (1.f - 0.f) / (aveVal - minVal)*(currentVal - minVal) + 0.f;
				Blue = 1.f - (1.f - 0.f) / (aveVal - minVal)*(currentVal - minVal);
			}
			else if ((currentVal > aveVal) and (currentVal <= maxVal)) {
				Red = (1.f - 0.f) / (maxVal - aveVal)*(currentVal - aveVal) + 0.f;
				Green = 1.f - (1.f - 0.f) / (maxVal - aveVal)*(currentVal - aveVal);
				Blue = 0.f;
			}
			if (minVal == maxVal) {
				Red = 0.f;
				Green = 1.f;
				Blue = 0.f;
			}
			i->m_color = glm::vec3(Red, Green, Blue);
		}
	}
	}
	else if (getColorParam() == "d(Vx)/dt") {
	changeColorParamTo("d(Vx)/dt");
	float maxVal = SPH.Particles[0]->m_velocity.dval.x;
	float minVal = SPH.Particles[0]->m_velocity.dval.x;
	for (auto*& i : SPH.Particles) {
		if (COLORING_CONDITION(i)) {
			if (i->m_velocity.dval.x < minVal) { minVal = i->m_velocity.dval.x; }
			if (i->m_velocity.dval.x > maxVal) { maxVal = i->m_velocity.dval.x; }
		}
	}
	float aveVal = (minVal + maxVal) / 2;
	m_maxVal = maxVal;
	m_minVal = minVal;

	for (auto*& i : SPH.Particles) {

		if (COLORING_CONDITION(i)) {
			currentVal = i->m_velocity.dval.x;

			if ((currentVal >= minVal) and (currentVal <= aveVal)) {
				Red = 0.f;
				Green = (1.f - 0.f) / (aveVal - minVal)*(currentVal - minVal) + 0.f;
				Blue = 1.f - (1.f - 0.f) / (aveVal - minVal)*(currentVal - minVal);
			}
			else if ((currentVal > aveVal) and (currentVal <= maxVal)) {
				Red = (1.f - 0.f) / (maxVal - aveVal)*(currentVal - aveVal) + 0.f;
				Green = 1.f - (1.f - 0.f) / (maxVal - aveVal)*(currentVal - aveVal);
				Blue = 0.f;
			}
			if (minVal == maxVal) {
				Red = 0.f;
				Green = 1.f;
				Blue = 0.f;
			}
			i->m_color = glm::vec3(Red, Green, Blue);
		}
	}
	}
	else if (getColorParam() == "d(Vy)/dt") {
	changeColorParamTo("d(Vy)/dt");
	float maxVal = SPH.Particles[0]->m_velocity.dval.y;
	float minVal = SPH.Particles[0]->m_velocity.dval.y;
	for (auto*& i : SPH.Particles) {
		if (COLORING_CONDITION(i)) {
			if (i->m_velocity.dval.y < minVal) { minVal = i->m_velocity.dval.y; }
			if (i->m_velocity.dval.y > maxVal) { maxVal = i->m_velocity.dval.y; }
		}
	}
	float aveVal = (minVal + maxVal) / 2;
	m_maxVal = maxVal;
	m_minVal = minVal;

	for (auto*& i : SPH.Particles) {

		if (COLORING_CONDITION(i)) {
			currentVal = i->m_velocity.dval.y;

			if ((currentVal >= minVal) and (currentVal <= aveVal)) {
				Red = 0.f;
				Green = (1.f - 0.f) / (aveVal - minVal)*(currentVal - minVal) + 0.f;
				Blue = 1.f - (1.f - 0.f) / (aveVal - minVal)*(currentVal - minVal);
			}
			else if ((currentVal > aveVal) and (currentVal <= maxVal)) {
				Red = (1.f - 0.f) / (maxVal - aveVal)*(currentVal - aveVal) + 0.f;
				Green = 1.f - (1.f - 0.f) / (maxVal - aveVal)*(currentVal - aveVal);
				Blue = 0.f;
			}
			if (minVal == maxVal) {
				Red = 0.f;
				Green = 1.f;
				Blue = 0.f;
			}
			i->m_color = glm::vec3(Red, Green, Blue);
		}
	}
	}
	else if (getColorParam() == "Density") {
	changeColorParamTo("Density");
	float maxVal = SPH.Particles[0]->m_density.val;
	float minVal = SPH.Particles[0]->m_density.val;
	for (auto*& i : SPH.Particles) {
		if (COLORING_CONDITION(i)) {
			if (i->m_density.val < minVal) { minVal = i->m_density.val; }
			if (i->m_density.val > maxVal) { maxVal = i->m_density.val; }
		}
	}
	float aveVal = (minVal + maxVal) / 2;
	m_maxVal = maxVal;
	m_minVal = minVal;

	for (auto*& i : SPH.Particles) {

		if (COLORING_CONDITION(i)) {
			currentVal = i->m_density.val;

			if ((currentVal >= minVal) and (currentVal <= aveVal)) {
				Red = 0.f;
				Green = (1.f - 0.f) / (aveVal - minVal)*(currentVal - minVal) + 0.f;
				Blue = 1.f - (1.f - 0.f) / (aveVal - minVal)*(currentVal - minVal);
			}
			else if ((currentVal > aveVal) and (currentVal <= maxVal)) {
				Red = (1.f - 0.f) / (maxVal - aveVal)*(currentVal - aveVal) + 0.f;
				Green = 1.f - (1.f - 0.f) / (maxVal - aveVal)*(currentVal - aveVal);
				Blue = 0.f;
			}
			if (minVal == maxVal) {
				Red = 0.f;
				Green = 1.f;
				Blue = 0.f;
			}
			i->m_color = glm::vec3(Red, Green, Blue);
		}
	}
	}
	else if (getColorParam() == "d(Density)/dt") {
	changeColorParamTo("d(Density)/dt");
	float maxVal = SPH.Particles[0]->m_density.dval;
	float minVal = SPH.Particles[0]->m_density.dval;
	for (auto*& i : SPH.Particles) {
		if (COLORING_CONDITION(i)) {
			if (i->m_density.dval < minVal) { minVal = i->m_density.dval; }
			if (i->m_density.dval > maxVal) { maxVal = i->m_density.dval; }
		}
	}
	float aveVal = (minVal + maxVal) / 2;
	m_maxVal = maxVal;
	m_minVal = minVal;

	for (auto*& i : SPH.Particles) {

		if (COLORING_CONDITION(i)) {
			currentVal = i->m_density.dval;

			if ((currentVal >= minVal) and (currentVal <= aveVal)) {
				Red = 0.f;
				Green = (1.f - 0.f) / (aveVal - minVal)*(currentVal - minVal) + 0.f;
				Blue = 1.f - (1.f - 0.f) / (aveVal - minVal)*(currentVal - minVal);
			}
			else if ((currentVal > aveVal) and (currentVal <= maxVal)) {
				Red = (1.f - 0.f) / (maxVal - aveVal)*(currentVal - aveVal) + 0.f;
				Green = 1.f - (1.f - 0.f) / (maxVal - aveVal)*(currentVal - aveVal);
				Blue = 0.f;
			}
			if (minVal == maxVal) {
				Red = 0.f;
				Green = 1.f;
				Blue = 0.f;
			}
			i->m_color = glm::vec3(Red, Green, Blue);
		}
	}
	}
	else if (getColorParam() == "pressure") {
	changeColorParamTo("pressure");
	float maxVal = SPH.Particles[0]->m_pressure.val;
	float minVal = SPH.Particles[0]->m_pressure.val;
	for (auto*& i : SPH.Particles) {
		if (COLORING_CONDITION(i)) {
			if (i->m_pressure.val < minVal) { minVal = i->m_pressure.val; }
			if (i->m_pressure.val > maxVal) { maxVal = i->m_pressure.val; }
		}
	}
	float aveVal = (minVal + maxVal) / 2;
	m_maxVal = maxVal;
	m_minVal = minVal;

	for (auto*& i : SPH.Particles) {

		if (COLORING_CONDITION(i)) {
			currentVal = i->m_pressure.val;

			if ((currentVal >= minVal) and (currentVal <= aveVal)) {
				Red = 0.f;
				Green = (1.f - 0.f) / (aveVal - minVal)*(currentVal - minVal) + 0.f;
				Blue = 1.f - (1.f - 0.f) / (aveVal - minVal)*(currentVal - minVal);
			}
			else if ((currentVal > aveVal) and (currentVal <= maxVal)) {
				Red = (1.f - 0.f) / (maxVal - aveVal)*(currentVal - aveVal) + 0.f;
				Green = 1.f - (1.f - 0.f) / (maxVal - aveVal)*(currentVal - aveVal);
				Blue = 0.f;
			}
			if (minVal == maxVal) {
				Red = 0.f;
				Green = 1.f;
				Blue = 0.f;
			}
			i->m_color = glm::vec3(Red, Green, Blue);
		}
	}
	}
	else if (getColorParam() == "nrOfNeighbours") {
	changeColorParamTo("nrOfNeighbours");
	float maxVal = SPH.Particles[0]->nrOfNeighbours;
	float minVal = SPH.Particles[0]->nrOfNeighbours;
	for (auto*& i : SPH.Particles) {
		if (COLORING_CONDITION(i)) {
			if (i->nrOfNeighbours < minVal) { minVal = i->nrOfNeighbours; }
			if (i->nrOfNeighbours > maxVal) { maxVal = i->nrOfNeighbours; }
		}
	}
	float aveVal = (minVal + maxVal) / 2;
	m_maxVal = maxVal;
	m_minVal = minVal;

	for (auto*& i : SPH.Particles) {

		if (COLORING_CONDITION(i)) {
			currentVal = i->nrOfNeighbours;

			if ((currentVal >= minVal) and (currentVal <= aveVal)) {
				Red = 0.f;
				Green = (1.f - 0.f) / (aveVal - minVal)*(currentVal - minVal) + 0.f;
				Blue = 1.f - (1.f - 0.f) / (aveVal - minVal)*(currentVal - minVal);
			}
			else if ((currentVal > aveVal) and (currentVal <= maxVal)) {
				Red = (1.f - 0.f) / (maxVal - aveVal)*(currentVal - aveVal) + 0.f;
				Green = 1.f - (1.f - 0.f) / (maxVal - aveVal)*(currentVal - aveVal);
				Blue = 0.f;
			}
			if (minVal == maxVal) {
				Red = 0.f;
				Green = 1.f;
				Blue = 0.f;
			}
			i->m_color = glm::vec3(Red, Green, Blue);
		}
	}
	}
	else if (getColorParam() == "mass") {
	changeColorParamTo("mass");
	float maxVal = SPH.Particles[0]->m_mass;
	float minVal = SPH.Particles[0]->m_mass;
	for (auto*& i : SPH.Particles) {
		if (COLORING_CONDITION(i)) {
			if (i->m_mass < minVal) { minVal = i->m_mass; }
			if (i->m_mass > maxVal) { maxVal = i->m_mass; }
		}
	}
	float aveVal = (minVal + maxVal) / 2;
	m_maxVal = maxVal;
	m_minVal = minVal;

	for (auto*& i : SPH.Particles) {

		if (COLORING_CONDITION(i)) {
			currentVal = i->m_mass;

			if ((currentVal >= minVal) and (currentVal <= aveVal)) {
				Red = 0.f;
				Green = (1.f - 0.f) / (aveVal - minVal)*(currentVal - minVal) + 0.f;
				Blue = 1.f - (1.f - 0.f) / (aveVal - minVal)*(currentVal - minVal);
			}
			else if ((currentVal > aveVal) and (currentVal <= maxVal)) {
				Red = (1.f - 0.f) / (maxVal - aveVal)*(currentVal - aveVal) + 0.f;
				Green = 1.f - (1.f - 0.f) / (maxVal - aveVal)*(currentVal - aveVal);
				Blue = 0.f;
			}
			if (minVal == maxVal) {
				Red = 0.f;
				Green = 1.f;
				Blue = 0.f;
			}
			i->m_color = glm::vec3(Red, Green, Blue);
		}
	}
	}
	//else {
	//	ColoringUnknow = true;
	//	changeColorParamTo("id");
	//	float maxVal = SPH.Particles[0]->m_id;
	//	float minVal = SPH.Particles[0]->m_id;
	//	for (auto*& i : SPH.Particles) {
	//		if (COLORING_CONDITION(i)) {
	//			if (i->m_id < minVal) { minVal = i->m_id; }
	//			if (i->m_id > maxVal) { maxVal = i->m_id; }
	//		}
	//	}
	//	float aveVal = (minVal + maxVal) / 2;
	//	m_maxVal = maxVal;
	//	m_minVal = minVal;
	//
	//	for (auto*& i : SPH.Particles) {
	//
	//		if (COLORING_CONDITION(i)) {
	//			currentVal = i->m_id;
	//
	//			if ((currentVal >= minVal) and (currentVal <= aveVal)) {
	//				Red = 0.f;
	//				Green = (1.f - 0.f) / (aveVal - minVal)*(currentVal - minVal) + 0.f;
	//				Blue = 1.f - (1.f - 0.f) / (aveVal - minVal)*(currentVal - minVal);
	//			}
	//			else if ((currentVal > aveVal) and (currentVal <= maxVal)) {
	//				Red = (1.f - 0.f) / (maxVal - aveVal)*(currentVal - aveVal) + 0.f;
	//				Green = 1.f - (1.f - 0.f) / (maxVal - aveVal)*(currentVal - aveVal);
	//				Blue = 0.f;
	//			}
	//			if (minVal == maxVal) {
	//				Red = 0.f;
	//				Green = 1.f;
	//				Blue = 0.f;
	//			}
	//			i->m_color = glm::vec3(Red, Green, Blue);
	//		}
	//	}
	//}
	else {
		changeColorParamTo("nrOfNeighbours");
		float maxVal = SPH.Particles[0]->nrOfNeighbours;
		float minVal = SPH.Particles[0]->nrOfNeighbours;
		for (auto*& i : SPH.Particles) {
			if (COLORING_CONDITION(i)) {
				if (i->nrOfNeighbours < minVal) { minVal = i->nrOfNeighbours; }
				if (i->nrOfNeighbours > maxVal) { maxVal = i->nrOfNeighbours; }
			}
		}
		float aveVal = (minVal + maxVal) / 2;
		m_maxVal = maxVal;
		m_minVal = minVal;

		for (auto*& i : SPH.Particles) {

			if (COLORING_CONDITION(i)) {
				currentVal = i->nrOfNeighbours;

				if ((currentVal >= minVal) and (currentVal <= aveVal)) {
					Red = 0.f;
					Green = (1.f - 0.f) / (aveVal - minVal)*(currentVal - minVal) + 0.f;
					Blue = 1.f - (1.f - 0.f) / (aveVal - minVal)*(currentVal - minVal);
				}
				else if ((currentVal > aveVal) and (currentVal <= maxVal)) {
					Red = (1.f - 0.f) / (maxVal - aveVal)*(currentVal - aveVal) + 0.f;
					Green = 1.f - (1.f - 0.f) / (maxVal - aveVal)*(currentVal - aveVal);
					Blue = 0.f;
				}
				if (minVal == maxVal) {
					Red = 0.f;
					Green = 1.f;
					Blue = 0.f;
				}
				i->m_color = glm::vec3(Red, Green, Blue);
			}
		}
	}
}

void SPH_CD::ColoringBtType() {
	if (computationalDomain_mode == MODE_CD::CD_DEBUG) {
		std::cout << "SPH_CD::Coloring\n";
	}
	for (auto*& i : SPH.Particles) {
		if (i->m_type == REAL)
			i->m_color = glm::vec3(0.f, 0.f, 1.f);
		else if (i->m_type == BOUNDARY)
			i->m_color = glm::vec3(1.f, 0.f, 0.f);
		else if (i->m_type == VIRTUAL)
			i->m_color = glm::vec3(0.f, 1.f, 0.f);

		//std::cout << currentVal << "   " << i->m_color.x << ", " << i->m_color.y << ", " << i->m_color.z << "\n";
	}
}

void SPH_CD::punctualColorChange(int number, glm::vec3 color) {

	if ((number < m_options.nrOfParticles[PARTICLETYPE::REAL] + m_options.nrOfParticles[PARTICLETYPE::BOUNDARY] + m_options.nrOfParticles[PARTICLETYPE::VIRTUAL]) and (number >= 0)) {
		//PRB[number].color = color;
		PRB[number].size = 4.0*SPH.Particles[number]->m_SmR / InitialSmR * 0.025 *0.125;
	}
}
/*
void SPH_CD::SPH_Iterations() {
	if (computationalDomain_mode == MODE_CD::CD_DEBUG) {
		std::cout << "SPH_CD::SPH_Iterations\n";
	}
	part_prec maxI_Wij = -1.0; part_prec maxI_dWij = -1.0; part_prec maxI_d2Wij = -1.0;
	part_prec sumI_Wij = 0.0; part_prec sumI_dWij = 0.0; part_prec sumI_d2Wij = 0.0;

	part_prec maxI_Wij_dr = -1.0; part_prec maxI_dWij_dr = -1.0; part_prec maxI_d2Wij_dr = -1.0;
	part_prec sumI_Wij_dr = 0.0; part_prec sumI_dWij_dr = 0.0; part_prec sumI_d2Wij_dr = 0.0;

	part_prec maxI_Wij_dr2 = -1.0; part_prec maxI_dWij_dr2 = -1.0; part_prec maxI_d2Wij_dr2 = -1.0;
	part_prec sumI_Wij_dr2 = 0.0; part_prec sumI_dWij_dr2 = 0.0; part_prec sumI_d2Wij_dr2 = 0.0;

	int m_size = 0;
	for (int i = 1;i <= nrOfDim + 1;i++) { m_size += i; }



	//while (true) {

		//for (auto*& i : SPH.ParticlePairs) {
		//	safeDelete(&i);
		//}
		//SPH.ParticlePairs.resize(0);
		//neighbourSearch();
		//creatingVirtualParticles();
		//addingVirtualTOPairs();


		std::vector<part_prec> h;   h.resize(SPH.Particles.size(), 0.f);
		std::vector<part_prec> hx;  hx.resize(SPH.Particles.size(), 0.f);
		std::vector<part_prec> hy;  hy.resize(SPH.Particles.size(), 0.f);
		std::vector<part_prec> hxx; hxx.resize(SPH.Particles.size(), 0.f);
		std::vector<part_prec> hyy; hyy.resize(SPH.Particles.size(), 0.f);
		std::vector<part_prec> hxy; hxy.resize(SPH.Particles.size(), 0.f);


		Matrix aveIntegrals;
		aveIntegrals.resize(m_size, m_size);


		part_prec aveh = 0.0;
		part_prec avehx = 0.0;
		part_prec avehy = 0.0;
		part_prec avehxx = 0.0;
		part_prec avehyy = 0.0;
		part_prec avehxy = 0.0;

		part_prec avedrho = 0.0;


		if (m_options.SPH_iteration_algorithm == STANDART_SPH) {
			std::vector<part_prec> drho;
			std::vector<part_prec> dh;
			drho.resize(SPH.Particles.size(), 0.f);
			dh.resize(SPH.Particles.size(), 0.f);
			std::vector<part_prec_3> dv;
			dv.resize(SPH.Particles.size());
			for (auto*& i : SPH.ParticlePairs) {
				drho[i->Mpart_ptr->m_id] += (i->NBpart_ptr->m_mass) * i->W(i->Mpart_ptr);
				drho[i->NBpart_ptr->m_id] += (i->Mpart_ptr->m_mass) * i->W(i->NBpart_ptr);
				dh[i->Mpart_ptr->m_id] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val) * (i->NBpart_ptr->m_SmR) * i->W(i->Mpart_ptr);
				dh[i->NBpart_ptr->m_id] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val) * (i->Mpart_ptr->m_SmR) * i->W(i->NBpart_ptr);

				for (int d = 0;d < nrOfDim;d++) {
					dv[i->Mpart_ptr->m_id][d] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val)*(i->NBpart_ptr->m_velocity.val[d]) * i->W(i->Mpart_ptr);
					dv[i->NBpart_ptr->m_id][d] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val)*(i->Mpart_ptr->m_velocity.val[d]) * i->W(i->NBpart_ptr);
				}
			}
			for (auto*& i : SPH.Particles) {
				if (i->m_type == REAL) {
					i->m_density.val = drho[i->m_id];
					for (int d = 0;d < nrOfDim;d++) {
						i->m_velocity.val[d] = dv[i->m_id][d];
					}
					i->m_SmR = dh[i->m_id];
				}
			}
		}
		if (m_options.SPH_iteration_algorithm == CORRECTIVE_SPH) {


			std::vector<ParticlePair::dsmthFunc> FunckeydX{ &Particle::dx ,&Particle::dy ,&Particle::dz,&Particle::distance };


			   

			int dim_matrix_size = nrOfDim;

			std::vector<Matrix> FPM_matrix;
			FPM_matrix.resize(SPH.Particles.size(), Matrix()); // на удаление
			for (int i = 0;i < SPH.Particles.size();i++) { FPM_matrix[i].resize(dim_matrix_size, dim_matrix_size); }

			std::vector<std::vector<part_prec>> FPM_vecs;
			FPM_vecs.resize(SPH.Particles.size(), std::vector<part_prec>()); // на удаление
			for (auto &jj : FPM_vecs) { jj.resize(dim_matrix_size, 0.f); }

			std::vector<std::vector<part_prec>> FPM_results;
			FPM_results.resize(SPH.Particles.size(), std::vector<part_prec>()); // на удаление
			for (auto &jj : FPM_results) { jj.resize(dim_matrix_size, 0.f); }

			std::vector<FiniteParticleMethod> FPM;
			FPM.resize(SPH.Particles.size(), FiniteParticleMethod(nrOfDim));

			//std::vector<std::vector<FiniteParticleMethod>> FPM;
			//FPM.resize(nrOfDim); //nrOfParemetersToCalculate
			//for (auto i : FPM) {
			//	i.resize(SPH.Particles.size(), FiniteParticleMethod(nrOfDim));
			//}

			//dh = hi - hi' = Sum(hj*Wij*mj/rhoj)*(1/Sum(Wij*mj/rhoj)-1) - {Sum((hj-hi)*Wij,x*mj/rhoj)/Sum((xj-xi)*Wij,x*mj/rhoj)*Sum((xj-xi)*Wij*mj/rhoj)/Sum(Wij*mj/rhoj) +
			//                                                            + Sum((hj-hi)*Wij,y*mj/rhoj)/Sum((yj-yi)*Wij,y*mj/rhoj)*Sum((yj-yi)*Wij*mj/rhoj)/Sum(Wij*mj/rhoj) + 
			//                                                            + Sum((hj-hi)*Wij,z*mj/rhoj)/Sum((zj-zi)*Wij,z*mj/rhoj)*Sum((zj-zi)*Wij*mj/rhoj)/Sum(Wij*mj/rhoj) }
			std::vector<part_prec> drho; drho.resize(SPH.Particles.size(), 0.f);
			std::vector<part_prec> dh; dh.resize(SPH.Particles.size(), 0.f);

			std::vector<part_prec> ddh; ddh.resize(SPH.Particles.size(), 0.f);
			//std::vector<part_prec_3> hj_i;
			//hj_i.resize(SPH.Particles.size(), { 0.0, 0.0, 0.0 });
			std::vector<part_prec> W_denominator; W_denominator.resize(SPH.Particles.size(), 0.0);
			std::vector<part_prec_3> dr_W; dr_W.resize(SPH.Particles.size(), { 0.0, 0.0, 0.0 });
			std::vector<part_prec_3> dr2_W; dr2_W.resize(SPH.Particles.size(), { 0.0, 0.0, 0.0 });

			std::vector<part_prec_3> dWdr; dWdr.resize(SPH.Particles.size(), { 0.0, 0.0, 0.0 });
			std::vector<part_prec_3> d2Wdr2; d2Wdr2.resize(SPH.Particles.size(), { 0.0, 0.0, 0.0 });


			//std::vector<part_prec> dr_W; dr_W.resize(SPH.Particles.size(), 0.f);
			std::vector<part_prec> dr_dWdr; dr_dWdr.resize(SPH.Particles.size(), 0.f);  //   Delta_h = dh(1/W_denominator - 1) - (dr_W)/(W_denominator)*(hj_i_dWdr)/(dr_dWdr)
			std::vector<part_prec> hj_i_dWdr; hj_i_dWdr.resize(SPH.Particles.size(), 0.f);
			std::vector<part_prec> dr_d2Wdr2; dr_d2Wdr2.resize(SPH.Particles.size(), 0.f);
			//std::vector<part_prec_3> dx_Wx_denominator;
			//dx_Wx_denominator.resize(SPH.Particles.size(), { 0.0, 0.0, 0.0 });



			std::vector<part_prec> dr2_dWdr; dr2_dWdr.resize(SPH.Particles.size(), 0.f);
			std::vector<part_prec> dr2_d2Wdr2; dr2_d2Wdr2.resize(SPH.Particles.size(), 0.f);


			std::vector<part_prec_3> dv;
			dv.resize(SPH.Particles.size(), { 0.0, 0.0, 0.0 });
			for (auto*& i : SPH.ParticlePairs) {


				if (PPL_BOTH_PARTICLES_ARE(REAL) or PPL_ONE_OF_PARTICLE_IS_1_OTHER_IS_2(REAL, VIRTUAL) or PPL_ONE_OF_PARTICLE_IS_1_OTHER_IS_2(REAL, BOUNDARY)) {
					PPC_CorSPH_General(i, &FPM_matrix);
					PPC_CorSPH_vec(i, &FPM_vecs, "h");
					//if (i->Mpart_ptr->m_id == 500) std::cout << FPM_matrix[500];
					//if (i->NBpart_ptr->m_id == 500) std::cout << FPM_matrix[500];

					drho[i->Mpart_ptr->m_id] += (i->NBpart_ptr->m_mass) * i->W(i->Mpart_ptr);
					drho[i->NBpart_ptr->m_id] += (i->Mpart_ptr->m_mass) * i->W(i->NBpart_ptr);
					dh[i->Mpart_ptr->m_id] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val) * (i->NBpart_ptr->m_SmR) * i->W(i->Mpart_ptr);
					dh[i->NBpart_ptr->m_id] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val) * (i->Mpart_ptr->m_SmR) * i->W(i->NBpart_ptr);

					ddh[i->Mpart_ptr->m_id] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val) * (i->NBpart_ptr->m_SmR) * (i->dWd(R, i->Mpart_ptr));
					ddh[i->NBpart_ptr->m_id] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val) * (i->Mpart_ptr->m_SmR) * (i->dWd(R, i->NBpart_ptr));

					W_denominator[i->Mpart_ptr->m_id] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val) * i->W(i->Mpart_ptr);
					W_denominator[i->NBpart_ptr->m_id] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val) * i->W(i->NBpart_ptr);


					//std::cout << i->dWd(X, i->Mpart_ptr) << " , " << i->dWd(Y, i->Mpart_ptr) << " , " << i->dWd(R, i->Mpart_ptr) << "\n";

					//std::cout << i->dWd(R, i->Mpart_ptr)<<"\n";
					dr_dWdr[i->Mpart_ptr->m_id] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val) *(i->dvfunc(FunckeydX[R]))* (i->dWd(R, i->Mpart_ptr));
					dr_dWdr[i->NBpart_ptr->m_id] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val) *i->dvfunc(FunckeydX[R])* (i->dWd(R, i->NBpart_ptr));

					hj_i_dWdr[i->Mpart_ptr->m_id] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val) *(i->NBpart_ptr->m_SmR - i->Mpart_ptr->m_SmR)* (i->dWd(R, i->Mpart_ptr));
					hj_i_dWdr[i->NBpart_ptr->m_id] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val) *(i->Mpart_ptr->m_SmR - i->NBpart_ptr->m_SmR)* (i->dWd(R, i->NBpart_ptr));


					dr_d2Wdr2[i->Mpart_ptr->m_id] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val)*(i->dvfunc(FunckeydX[R])) * (i->d2Wd(R, R, i->Mpart_ptr));
					dr_d2Wdr2[i->NBpart_ptr->m_id] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val)*i->dvfunc(FunckeydX[R]) * (i->d2Wd(R, R, i->NBpart_ptr));




					dr2_dWdr[i->Mpart_ptr->m_id] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val) *pow(i->dvfunc(FunckeydX[R]), 2)* (i->dWd(R, i->Mpart_ptr));
					dr2_dWdr[i->NBpart_ptr->m_id] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val) *pow(i->dvfunc(FunckeydX[R]), 2)* (i->dWd(R, i->NBpart_ptr));

					dr2_d2Wdr2[i->Mpart_ptr->m_id] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val)*pow(i->dvfunc(FunckeydX[R]), 2) * (i->d2Wd(R, R, i->Mpart_ptr));
					dr2_d2Wdr2[i->NBpart_ptr->m_id] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val)*pow(i->dvfunc(FunckeydX[R]), 2) * (i->d2Wd(R, R, i->NBpart_ptr));


					for (int d = 0;d < nrOfDim;d++) {
						dv[i->Mpart_ptr->m_id][d] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val)*(i->NBpart_ptr->m_velocity.val[d]) * i->W(i->Mpart_ptr);
						dv[i->NBpart_ptr->m_id][d] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val)*(i->Mpart_ptr->m_velocity.val[d]) * i->W(i->NBpart_ptr);

						//   dx*W   dy*W   dz*W   
						dr_W[i->Mpart_ptr->m_id][d] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val)*(-(i->dvfunc(FunckeydX[d]))) * i->W(i->Mpart_ptr);
						dr_W[i->NBpart_ptr->m_id][d] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val)*((i->dvfunc(FunckeydX[d]))) * i->W(i->NBpart_ptr);

						//   dx2*W   dy2*W   dz2*W  
						dr2_W[i->Mpart_ptr->m_id] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val) *pow(i->dvfunc(FunckeydX[d]), 2)* i->W(i->Mpart_ptr);
						dr2_W[i->NBpart_ptr->m_id] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val) *pow(i->dvfunc(FunckeydX[d]), 2)* i->W(i->NBpart_ptr);

						//   dWdx   dWdy   dWdz  
						dWdr[i->Mpart_ptr->m_id][d] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val) * (i->dWd(static_cast<DiffAxis>(d), i->Mpart_ptr));
						dWdr[i->NBpart_ptr->m_id][d] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val) * (i->dWd(static_cast<DiffAxis>(d), i->NBpart_ptr));

						//   d2Wdx2   d2Wdy2   d2Wdz2  
						d2Wdr2[i->Mpart_ptr->m_id] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val) * (i->d2Wd(static_cast<DiffAxis>(d), static_cast<DiffAxis>(d), i->Mpart_ptr));
						d2Wdr2[i->NBpart_ptr->m_id] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val) * (i->d2Wd(static_cast<DiffAxis>(d), static_cast<DiffAxis>(d), i->NBpart_ptr));


						//dx_Wx_denominator[i->Mpart_ptr->m_id][d] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val)*(-(i->dvfunc(FunckeydX[d]))) * i->dWd(static_cast<DiffAxis>(d), i->Mpart_ptr);
						//dx_Wx_denominator[i->NBpart_ptr->m_id][d] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val)*((i->dvfunc(FunckeydX[d]))) * i->dWd(static_cast<DiffAxis>(d), i->NBpart_ptr);
						//hj_i[i->Mpart_ptr->m_id][d] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val)*(i->NBpart_ptr->m_SmR -i->Mpart_ptr->m_SmR) * i->dWd(static_cast<DiffAxis>(d), i->Mpart_ptr);
						//hj_i[i->NBpart_ptr->m_id][d] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val)*(i->Mpart_ptr->m_SmR - i->NBpart_ptr->m_SmR) * i->dWd(static_cast<DiffAxis>(d), i->NBpart_ptr);
					}
				}
			}
			for (auto*& i : SPH.Particles) {
				if (i->m_type == REAL) {


					aveIntegrals(0, 0) += W_denominator[i->m_id]; aveIntegrals(0, 1) += dWdr[i->m_id][0]; aveIntegrals(0, 2) += dWdr[i->m_id][1]; //aveIntegrals(0, 3) += d2Wdr2[i->m_id][0]; aveIntegrals(0, 4) += d2Wdr2[i->m_id][1]; aveIntegrals(0, 5) += d2Wdxdy[i->m_id];
					aveIntegrals(1, 0) += dr_W[i->m_id][0]; aveIntegrals(1, 1) += FPM_matrix[i->m_id](0, 0); aveIntegrals(1, 2) += FPM_matrix[i->m_id](1, 0); //aveIntegrals(0, 3) += d2Wdx2_dx[i->m_id]; aveIntegrals(0, 4) += d2Wdxy_dx[i->m_id]; aveIntegrals(0, 5) += d2Wdxdy_dx[i->m_id];
					aveIntegrals(2, 0) += dr_W[i->m_id][1]; aveIntegrals(2, 1) += FPM_matrix[i->m_id](0, 1); aveIntegrals(2, 2) += FPM_matrix[i->m_id](1, 1); //aveIntegrals(0, 3) += d2Wdx2_dy[i->m_id]; aveIntegrals(0, 4) += d2Wdxy_dy[i->m_id]; aveIntegrals(0, 5) += d2Wdxdy_dy[i->m_id];
					//aveIntegrals(0, 0) += dr2_W[i->m_id][0]; aveIntegrals(0, 1) += dWdx_dx2[i->m_id]; aveIntegrals(0, 2) += dWdy_dx2[i->m_id]; aveIntegrals(0, 3) += d2Wdx2_dx2[i->m_id]; aveIntegrals(0, 4) += d2Wdxy_dx2[i->m_id]; aveIntegrals(0, 5) += d2Wdxdy_dx2[i->m_id];
					//aveIntegrals(0, 0) += dr2_W[i->m_id][1]; aveIntegrals(0, 1) += dWdx_dy2[i->m_id]; aveIntegrals(0, 2) += dWdy_dy2[i->m_id]; aveIntegrals(0, 3) += d2Wdx2_dy2[i->m_id]; aveIntegrals(0, 4) += d2Wdxy_dy2[i->m_id]; aveIntegrals(0, 5) += d2Wdxdy_dy2[i->m_id];
					//aveIntegrals(0, 0) += W_dxdy[i->m_id]; aveIntegrals(0, 1) += dWdx_dxdy[i->m_id]; aveIntegrals(0, 2) += dWdy_dxdy[i->m_id]; aveIntegrals(0, 3) += d2Wdx2_dxdy[i->m_id]; aveIntegrals(0, 4) += d2Wdxy_dxdy[i->m_id]; aveIntegrals(0, 5) += d2Wdxdy_dxdy[i->m_id];


					//std::cout << FPM_matrix[i->m_id] << "\n";
					FPM[i->m_id].addMatVec(FPM_matrix[i->m_id], FPM_vecs[i->m_id]);
					FPM[i->m_id].Calculate();

					//std::cout << dr_dWdr[i->m_id] << "\n";

					// sumI_Wij += W_denominator[i->m_id];
					// if (std::abs(W_denominator[i->m_id] - 1.0) > maxI_Wij) maxI_Wij = std::abs(W_denominator[i->m_id] - 1.0);
					// sumI_dWij += dWdr[i->m_id];
					// if (std::abs(dWdr[i->m_id]) > maxI_dWij) maxI_dWij = std::abs(dWdr[i->m_id]);
					// sumI_d2Wij += d2Wdr2[i->m_id];
					// if (std::abs(d2Wdr2[i->m_id]) > maxI_d2Wij) maxI_d2Wij = std::abs(d2Wdr2[i->m_id]);
					// 
					// sumI_Wij_dr += dr_W[i->m_id];
					// if (std::abs(dr_W[i->m_id]) > maxI_Wij_dr) maxI_Wij_dr = std::abs(dr_W[i->m_id]);
					// sumI_dWij_dr += dr_dWdr[i->m_id];
					// if (std::abs(dr_dWdr[i->m_id] - 1.0) > maxI_dWij_dr) maxI_dWij_dr = std::abs(dr_dWdr[i->m_id] - 1.0);
					// sumI_d2Wij_dr += dr_d2Wdr2[i->m_id];
					// if (std::abs(dr_d2Wdr2[i->m_id]) > maxI_d2Wij_dr) maxI_d2Wij_dr = std::abs(dr_d2Wdr2[i->m_id]);
					// 
					// sumI_Wij_dr2 += dr2_W[i->m_id];
					// if (std::abs(dr2_W[i->m_id]) > maxI_Wij_dr2) maxI_Wij_dr2 = std::abs(dr2_W[i->m_id]);
					// sumI_dWij_dr2 += std::abs(dr2_dWdr[i->m_id]);
					// if (std::abs(dr2_dWdr[i->m_id]) > maxI_dWij_dr2) maxI_dWij_dr2 = std::abs(dr2_dWdr[i->m_id]);
					// sumI_d2Wij_dr2 += dr2_d2Wdr2[i->m_id];
					// if (std::abs(dr2_d2Wdr2[i->m_id]-2.0) > maxI_d2Wij_dr2) maxI_d2Wij_dr2 = std::abs(dr2_d2Wdr2[i->m_id]-2.0);

				}

				drho[i->m_id] = drho[i->m_id] / W_denominator[i->m_id];
				for (int d = 0;d < nrOfDim;d++) {
					dv[i->m_id][d] = dv[i->m_id][d] / W_denominator[i->m_id];
				}

				if (i->m_type == REAL) {
					//std::cout << i->m_density.val - drho[i->m_id] << "   ";
					i->m_density.val = drho[i->m_id];
					for (int d = 0;d < nrOfDim;d++) {
						i->m_velocity.val[d] = dv[i->m_id][d];
					}
					//std::cout << dh[i->m_id] << "/" << W_denominator[i->m_id] << "  ";
					dh[i->m_id] = dh[i->m_id] * (1 / W_denominator[i->m_id] - 1);
					//dh[i->m_id] -= hj_i_dWdr[i->m_id]/ dr_dWdr[i->m_id]* dr_W[i->m_id] / W_denominator[i->m_id];
					//std::cout << hj_i_dWdr[i->m_id] / dr_dWdr[i->m_id] * dr_W[i->m_id] / W_denominator[i->m_id] << "\n";
					for (int d = 0;d < nrOfDim;d++) {
						//std::cout << FPM[i->m_id].getd(static_cast<DIMENSIONS>(d)) << "  ";
						dh[i->m_id] -= FPM[i->m_id].getd(static_cast<DIMENSIONS>(d)) * dr_W[i->m_id][d] / W_denominator[i->m_id];
					}
					//std::cout << i->m_id << "  " << i->nrOfNeighbours << "  i->m_SmR =" << i->m_SmR << " + dh = " << dh[i->m_id] << "\n";
					i->m_SmR += dh[i->m_id];
				}
			}
		}
		if (m_options.SPH_iteration_algorithm == FINITE_PARTICLE_METHOD) {
			std::vector<ParticlePair::dsmthFunc> FunckeydX{ &Particle::dx ,&Particle::dy ,&Particle::dz,&Particle::distance };

			std::vector<Matrix> FPM_matrix;
			FPM_matrix.resize(SPH.Particles.size(), Matrix()); // на удаление
			for (int i = 0;i < SPH.Particles.size();i++) { FPM_matrix[i].resize(m_size, m_size); }

			std::vector<std::vector<part_prec>> FPM_vecs;
			FPM_vecs.resize(SPH.Particles.size(), std::vector<part_prec>()); // на удаление
			for (auto &jj : FPM_vecs) { jj.resize(m_size, 0.f); }

			std::vector<std::vector<part_prec>> FPM_results;
			FPM_results.resize(SPH.Particles.size(), std::vector<part_prec>()); // на удаление
			for (auto &jj : FPM_results) { jj.resize(m_size, 0.f); }

			std::vector<FiniteParticleMethod> FPM;
			FPM.resize(SPH.Particles.size(), FiniteParticleMethod(nrOfDim));

			//std::vector<std::vector<FiniteParticleMethod>> FPM;
			//FPM.resize(nrOfDim); //nrOfParemetersToCalculate
			//for (auto i : FPM) {
			//	i.resize(SPH.Particles.size(), FiniteParticleMethod(nrOfDim));
			//}

			//dh = hi - hi' = Sum(hj*Wij*mj/rhoj)*(1/Sum(Wij*mj/rhoj)-1) - {Sum((hj-hi)*Wij,x*mj/rhoj)/Sum((xj-xi)*Wij,x*mj/rhoj)*Sum((xj-xi)*Wij*mj/rhoj)/Sum(Wij*mj/rhoj) +
			//                                                            + Sum((hj-hi)*Wij,y*mj/rhoj)/Sum((yj-yi)*Wij,y*mj/rhoj)*Sum((yj-yi)*Wij*mj/rhoj)/Sum(Wij*mj/rhoj) + 
			//                                                            + Sum((hj-hi)*Wij,z*mj/rhoj)/Sum((zj-zi)*Wij,z*mj/rhoj)*Sum((zj-zi)*Wij*mj/rhoj)/Sum(Wij*mj/rhoj) }
			std::vector<part_prec> drho; drho.resize(SPH.Particles.size(), 0.f);
			std::vector<part_prec> dh; dh.resize(SPH.Particles.size(), 0.f);

			std::vector<part_prec> ddh; ddh.resize(SPH.Particles.size(), 0.f);
			//std::vector<part_prec_3> hj_i;
			//hj_i.resize(SPH.Particles.size(), { 0.0, 0.0, 0.0 });
			std::vector<part_prec> W_denominator; W_denominator.resize(SPH.Particles.size(), 0.0);
			std::vector<part_prec_3> dr_W; dr_W.resize(SPH.Particles.size(), { 0.0, 0.0, 0.0 });
			std::vector<part_prec_3> dr2_W; dr2_W.resize(SPH.Particles.size(), { 0.0, 0.0, 0.0 });

			std::vector<part_prec_3> dWdr; dWdr.resize(SPH.Particles.size(), { 0.0, 0.0, 0.0 });
			std::vector<part_prec_3> d2Wdr2; d2Wdr2.resize(SPH.Particles.size(), { 0.0, 0.0, 0.0 });


			//std::vector<part_prec> dr_W; dr_W.resize(SPH.Particles.size(), 0.f);
			std::vector<part_prec> dr_dWdr; dr_dWdr.resize(SPH.Particles.size(), 0.f);  //   Delta_h = dh(1/W_denominator - 1) - (dr_W)/(W_denominator)*(hj_i_dWdr)/(dr_dWdr)
			std::vector<part_prec> hj_i_dWdr; hj_i_dWdr.resize(SPH.Particles.size(), 0.f);
			std::vector<part_prec> dr_d2Wdr2; dr_d2Wdr2.resize(SPH.Particles.size(), 0.f);
			//std::vector<part_prec_3> dx_Wx_denominator;
			//dx_Wx_denominator.resize(SPH.Particles.size(), { 0.0, 0.0, 0.0 });
			std::vector<part_prec> dr2_dWdr; dr2_dWdr.resize(SPH.Particles.size(), 0.f);
			std::vector<part_prec> dr2_d2Wdr2; dr2_d2Wdr2.resize(SPH.Particles.size(), 0.f);


			std::vector<part_prec> hj_i_dWdx; hj_i_dWdx.resize(SPH.Particles.size(), 0.f);
			std::vector<part_prec> hj_i_dWdy; hj_i_dWdy.resize(SPH.Particles.size(), 0.f);
			std::vector<part_prec> hj_i_d2Wdx2; hj_i_d2Wdx2.resize(SPH.Particles.size(), 0.f);
			std::vector<part_prec> hj_i_d2Wdy2; hj_i_d2Wdy2.resize(SPH.Particles.size(), 0.f);
			std::vector<part_prec> hj_i_d2Wdxdy; hj_i_d2Wdxdy.resize(SPH.Particles.size(), 0.f);


			std::vector<part_prec_3> dv;
			dv.resize(SPH.Particles.size(), { 0.0, 0.0, 0.0 });
			for (auto*& i : SPH.ParticlePairs) {


				if (PPL_BOTH_PARTICLES_ARE(REAL) or PPL_ONE_OF_PARTICLE_IS_1_OTHER_IS_2(REAL, VIRTUAL) or PPL_ONE_OF_PARTICLE_IS_1_OTHER_IS_2(REAL, BOUNDARY)) {
					PPC_FPMStandart_General(i, &FPM_matrix);
					PPC_FPMStandart_vec(i, &FPM_vecs, "h");


					drho[i->Mpart_ptr->m_id] += (i->NBpart_ptr->m_mass) * i->W(i->Mpart_ptr);
					drho[i->NBpart_ptr->m_id] += (i->Mpart_ptr->m_mass) * i->W(i->NBpart_ptr);

					dh[i->Mpart_ptr->m_id] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val) * (i->NBpart_ptr->m_SmR) * i->W(i->Mpart_ptr);
					dh[i->NBpart_ptr->m_id] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val) * (i->Mpart_ptr->m_SmR) * i->W(i->NBpart_ptr);

					hj_i_dWdx[i->Mpart_ptr->m_id] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val) *(i->NBpart_ptr->m_SmR - i->Mpart_ptr->m_SmR)* (i->dWd(X, i->Mpart_ptr));
					hj_i_dWdx[i->NBpart_ptr->m_id] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val) *(i->Mpart_ptr->m_SmR - i->NBpart_ptr->m_SmR)* (i->dWd(X, i->NBpart_ptr));

					hj_i_dWdy[i->Mpart_ptr->m_id] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val) *(i->NBpart_ptr->m_SmR - i->Mpart_ptr->m_SmR)* (i->dWd(Y, i->Mpart_ptr));
					hj_i_dWdy[i->NBpart_ptr->m_id] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val) *(i->Mpart_ptr->m_SmR - i->NBpart_ptr->m_SmR)* (i->dWd(Y, i->NBpart_ptr));

					hj_i_d2Wdx2[i->Mpart_ptr->m_id] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val) *(i->NBpart_ptr->m_SmR - i->Mpart_ptr->m_SmR)* (i->d2Wd(X, X, i->Mpart_ptr));
					hj_i_d2Wdx2[i->NBpart_ptr->m_id] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val) *(i->Mpart_ptr->m_SmR - i->NBpart_ptr->m_SmR)*(i->d2Wd(X, X, i->Mpart_ptr));

					hj_i_d2Wdy2[i->Mpart_ptr->m_id] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val) *(i->NBpart_ptr->m_SmR - i->Mpart_ptr->m_SmR)* (i->d2Wd(Y, Y, i->Mpart_ptr));
					hj_i_d2Wdy2[i->NBpart_ptr->m_id] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val) *(i->Mpart_ptr->m_SmR - i->NBpart_ptr->m_SmR)* (i->d2Wd(Y, Y, i->Mpart_ptr));

					hj_i_d2Wdxdy[i->Mpart_ptr->m_id] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val) *(i->NBpart_ptr->m_SmR - i->Mpart_ptr->m_SmR)* (i->d2Wd(X, Y, i->Mpart_ptr));
					hj_i_d2Wdxdy[i->NBpart_ptr->m_id] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val) *(i->Mpart_ptr->m_SmR - i->NBpart_ptr->m_SmR)* (i->d2Wd(X, Y, i->Mpart_ptr));

					
					drho[i->Mpart_ptr->m_id] += (i->NBpart_ptr->m_mass) * i->W(i->Mpart_ptr);
					drho[i->NBpart_ptr->m_id] += (i->Mpart_ptr->m_mass) * i->W(i->NBpart_ptr);
					dh[i->Mpart_ptr->m_id] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val) * (i->NBpart_ptr->m_SmR) * i->W(i->Mpart_ptr);
					dh[i->NBpart_ptr->m_id] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val) * (i->Mpart_ptr->m_SmR) * i->W(i->NBpart_ptr);

					ddh[i->Mpart_ptr->m_id] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val) * (i->NBpart_ptr->m_SmR) * (i->dWd(R, i->Mpart_ptr));
					ddh[i->NBpart_ptr->m_id] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val) * (i->Mpart_ptr->m_SmR) * (i->dWd(R, i->NBpart_ptr));

					W_denominator[i->Mpart_ptr->m_id] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val) * i->W(i->Mpart_ptr);
					W_denominator[i->NBpart_ptr->m_id] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val) * i->W(i->NBpart_ptr);


					//std::cout << i->dWd(X, i->Mpart_ptr) << " , " << i->dWd(Y, i->Mpart_ptr) << " , " << i->dWd(R, i->Mpart_ptr) << "\n";

					//std::cout << i->dWd(R, i->Mpart_ptr)<<"\n";
					dr_dWdr[i->Mpart_ptr->m_id] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val) *(i->dvfunc(FunckeydX[R]))* (i->dWd(R, i->Mpart_ptr));
					dr_dWdr[i->NBpart_ptr->m_id] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val) *i->dvfunc(FunckeydX[R])* (i->dWd(R, i->NBpart_ptr));

					hj_i_dWdr[i->Mpart_ptr->m_id] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val) *(i->NBpart_ptr->m_SmR - i->Mpart_ptr->m_SmR)* (i->dWd(R, i->Mpart_ptr));
					hj_i_dWdr[i->NBpart_ptr->m_id] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val) *(i->Mpart_ptr->m_SmR - i->NBpart_ptr->m_SmR)* (i->dWd(R, i->NBpart_ptr));


					dr_d2Wdr2[i->Mpart_ptr->m_id] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val)*(i->dvfunc(FunckeydX[R])) * (i->d2Wd(R, R, i->Mpart_ptr));
					dr_d2Wdr2[i->NBpart_ptr->m_id] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val)*i->dvfunc(FunckeydX[R]) * (i->d2Wd(R, R, i->NBpart_ptr));




					dr2_dWdr[i->Mpart_ptr->m_id] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val) *pow(i->dvfunc(FunckeydX[R]), 2)* (i->dWd(R, i->Mpart_ptr));
					dr2_dWdr[i->NBpart_ptr->m_id] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val) *pow(i->dvfunc(FunckeydX[R]), 2)* (i->dWd(R, i->NBpart_ptr));

					dr2_d2Wdr2[i->Mpart_ptr->m_id] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val)*pow(i->dvfunc(FunckeydX[R]), 2) * (i->d2Wd(R, R, i->Mpart_ptr));
					dr2_d2Wdr2[i->NBpart_ptr->m_id] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val)*pow(i->dvfunc(FunckeydX[R]), 2) * (i->d2Wd(R, R, i->NBpart_ptr));


					for (int d = 0;d < nrOfDim;d++) {
						dv[i->Mpart_ptr->m_id][d] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val)*(i->NBpart_ptr->m_velocity.val[d]) * i->W(i->Mpart_ptr);
						dv[i->NBpart_ptr->m_id][d] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val)*(i->Mpart_ptr->m_velocity.val[d]) * i->W(i->NBpart_ptr);

						//   dx*W   dy*W   dz*W
						dr_W[i->Mpart_ptr->m_id][d] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val)*(-(i->dvfunc(FunckeydX[d]))) * i->W(i->Mpart_ptr);
						dr_W[i->NBpart_ptr->m_id][d] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val)*((i->dvfunc(FunckeydX[d]))) * i->W(i->NBpart_ptr);

						//   dx2*W   dy2*W   dz2*W
						dr2_W[i->Mpart_ptr->m_id] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val) *pow(i->dvfunc(FunckeydX[d]), 2)* i->W(i->Mpart_ptr);
						dr2_W[i->NBpart_ptr->m_id] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val) *pow(i->dvfunc(FunckeydX[d]), 2)* i->W(i->NBpart_ptr);

						//   dWdx   dWdy   dWdz
						dWdr[i->Mpart_ptr->m_id][d] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val) * (i->dWd(static_cast<DiffAxis>(d), i->Mpart_ptr));
						dWdr[i->NBpart_ptr->m_id][d] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val) * (i->dWd(static_cast<DiffAxis>(d), i->NBpart_ptr));

						//   d2Wdx2   d2Wdy2   d2Wdz2
						d2Wdr2[i->Mpart_ptr->m_id] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val) * (i->d2Wd(static_cast<DiffAxis>(d), static_cast<DiffAxis>(d), i->Mpart_ptr));
						d2Wdr2[i->NBpart_ptr->m_id] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val) * (i->d2Wd(static_cast<DiffAxis>(d), static_cast<DiffAxis>(d), i->NBpart_ptr));


						//dx_Wx_denominator[i->Mpart_ptr->m_id][d] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val)*(-(i->dvfunc(FunckeydX[d]))) * i->dWd(static_cast<DiffAxis>(d), i->Mpart_ptr);
						//dx_Wx_denominator[i->NBpart_ptr->m_id][d] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val)*((i->dvfunc(FunckeydX[d]))) * i->dWd(static_cast<DiffAxis>(d), i->NBpart_ptr);
						//hj_i[i->Mpart_ptr->m_id][d] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val)*(i->NBpart_ptr->m_SmR -i->Mpart_ptr->m_SmR) * i->dWd(static_cast<DiffAxis>(d), i->Mpart_ptr);
						//hj_i[i->NBpart_ptr->m_id][d] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val)*(i->Mpart_ptr->m_SmR - i->NBpart_ptr->m_SmR) * i->dWd(static_cast<DiffAxis>(d), i->NBpart_ptr);
					}
					
				}
			}


			for (auto*& i : SPH.Particles) {
				if (i->m_type == REAL) {
					//FPM[i->m_id].addMatVec(FPM_matrix[i->m_id], FPM_vecs[i->m_id]);
					//FPM[i->m_id].Calculate();
					//h[i->m_id] = dh[i->m_id] / FPM_matrix[i->m_id](0, 0) - FPM[i->m_id].getd(static_cast<DIMENSIONS>(1)) * FPM_matrix[i->m_id](0, 1) / FPM_matrix[i->m_id](0, 0) - FPM[i->m_id].getd(static_cast<DIMENSIONS>(2)) * FPM_matrix[i->m_id](0, 2) / FPM_matrix[i->m_id](0, 0) - FPM[i->m_id].getd(static_cast<DIMENSIONS>(1), static_cast<DIMENSIONS>(1)) * FPM_matrix[i->m_id](0, 3) / FPM_matrix[i->m_id](0, 0) - FPM[i->m_id].getd(static_cast<DIMENSIONS>(2), static_cast<DIMENSIONS>(2)) * FPM_matrix[i->m_id](0, 4) / FPM_matrix[i->m_id](0, 0) - FPM[i->m_id].getd(static_cast<DIMENSIONS>(1), static_cast<DIMENSIONS>(2)) * FPM_matrix[i->m_id](0, 5) / FPM_matrix[i->m_id](0, 0);
					//hx[i->m_id] = hj_i_dWdx[i->m_id] / FPM_matrix[i->m_id](1, 1) - hy[i->m_id] * FPM_matrix[i->m_id](1, 2) / FPM_matrix[i->m_id](1, 1) - FPM[i->m_id].getd(static_cast<DIMENSIONS>(1), static_cast<DIMENSIONS>(1)) * FPM_matrix[i->m_id](1, 3) / FPM_matrix[i->m_id](1, 1) - FPM[i->m_id].getd(static_cast<DIMENSIONS>(2), static_cast<DIMENSIONS>(2)) * FPM_matrix[i->m_id](1, 4) / FPM_matrix[i->m_id](1, 1) - FPM[i->m_id].getd(static_cast<DIMENSIONS>(1), static_cast<DIMENSIONS>(2)) * FPM_matrix[i->m_id](1, 5) / FPM_matrix[i->m_id](1, 1);
					//hy[i->m_id] = hj_i_dWdy[i->m_id] / FPM_matrix[i->m_id](2, 2) - FPM[i->m_id].getd(static_cast<DIMENSIONS>(1)) * FPM_matrix[i->m_id](2, 1) / FPM_matrix[i->m_id](2, 2) - FPM[i->m_id].getd(static_cast<DIMENSIONS>(1), static_cast<DIMENSIONS>(1)) * FPM_matrix[i->m_id](2, 3) / FPM_matrix[i->m_id](2, 2) - FPM[i->m_id].getd(static_cast<DIMENSIONS>(2), static_cast<DIMENSIONS>(2)) * FPM_matrix[i->m_id](2, 4) / FPM_matrix[i->m_id](2, 2) - FPM[i->m_id].getd(static_cast<DIMENSIONS>(1), static_cast<DIMENSIONS>(2)) * FPM_matrix[i->m_id](2, 5) / FPM_matrix[i->m_id](2, 2);
					//hxx[i->m_id] = hj_i_d2Wdx2[i->m_id] / FPM_matrix[i->m_id](3, 3) - FPM[i->m_id].getd(static_cast<DIMENSIONS>(1)) * FPM_matrix[i->m_id](3, 1) / FPM_matrix[i->m_id](3, 3) - FPM[i->m_id].getd(static_cast<DIMENSIONS>(2)) * FPM_matrix[i->m_id](3, 2) / FPM_matrix[i->m_id](3, 3) - FPM[i->m_id].getd(static_cast<DIMENSIONS>(2), static_cast<DIMENSIONS>(2)) * FPM_matrix[i->m_id](3, 4) / FPM_matrix[i->m_id](3, 3) - FPM[i->m_id].getd(static_cast<DIMENSIONS>(1), static_cast<DIMENSIONS>(2)) * FPM_matrix[i->m_id](3, 5) / FPM_matrix[i->m_id](3, 3);
					//hyy[i->m_id] = hj_i_d2Wdy2[i->m_id] / FPM_matrix[i->m_id](4, 4) - FPM[i->m_id].getd(static_cast<DIMENSIONS>(1)) * FPM_matrix[i->m_id](4, 1) / FPM_matrix[i->m_id](4, 4) - FPM[i->m_id].getd(static_cast<DIMENSIONS>(2)) * FPM_matrix[i->m_id](4, 2) / FPM_matrix[i->m_id](4, 4) - FPM[i->m_id].getd(static_cast<DIMENSIONS>(1), static_cast<DIMENSIONS>(1)) * FPM_matrix[i->m_id](4, 3) / FPM_matrix[i->m_id](4, 4) - FPM[i->m_id].getd(static_cast<DIMENSIONS>(1), static_cast<DIMENSIONS>(2)) * FPM_matrix[i->m_id](4, 5) / FPM_matrix[i->m_id](4, 4);
					//hxy[i->m_id] = hj_i_d2Wdxdy[i->m_id] / FPM_matrix[i->m_id](5, 5) - FPM[i->m_id].getd(static_cast<DIMENSIONS>(1)) * FPM_matrix[i->m_id](5, 1) / FPM_matrix[i->m_id](5, 5) - FPM[i->m_id].getd(static_cast<DIMENSIONS>(2)) * FPM_matrix[i->m_id](5, 2) / FPM_matrix[i->m_id](5, 5) - FPM[i->m_id].getd(static_cast<DIMENSIONS>(1), static_cast<DIMENSIONS>(1)) * FPM_matrix[i->m_id](5, 3) / FPM_matrix[i->m_id](5, 5) - FPM[i->m_id].getd(static_cast<DIMENSIONS>(2), static_cast<DIMENSIONS>(2)) * FPM_matrix[i->m_id](5, 4) / FPM_matrix[i->m_id](5, 5);



					

					h[i->m_id] = dh[i->m_id] *(1.0 / FPM_matrix[i->m_id](0, 0) - 1.0) - hx[i->m_id] * FPM_matrix[i->m_id](0, 1) / FPM_matrix[i->m_id](0, 0) - hy[i->m_id] * FPM_matrix[i->m_id](0, 2) / FPM_matrix[i->m_id](0, 0) - hxx[i->m_id] * FPM_matrix[i->m_id](0, 3) / FPM_matrix[i->m_id](0, 0) - hyy[i->m_id] * FPM_matrix[i->m_id](0, 4) / FPM_matrix[i->m_id](0, 0) - hxy[i->m_id] * FPM_matrix[i->m_id](0, 5) / FPM_matrix[i->m_id](0, 0);
					hx[i->m_id] = hj_i_dWdx[i->m_id] * (1.0 / FPM_matrix[i->m_id](1, 1) - 1.0) - hy[i->m_id] * FPM_matrix[i->m_id](1, 2) / FPM_matrix[i->m_id](1, 1) - hxx[i->m_id] * FPM_matrix[i->m_id](1, 3) / FPM_matrix[i->m_id](1, 1) - hyy[i->m_id] * FPM_matrix[i->m_id](1, 4) / FPM_matrix[i->m_id](1, 1) - hxy[i->m_id] * FPM_matrix[i->m_id](1, 5) / FPM_matrix[i->m_id](1, 1);
					hy[i->m_id] = hj_i_dWdy[i->m_id] * (1.0 / FPM_matrix[i->m_id](2, 2) - 1.0) - hx[i->m_id] * FPM_matrix[i->m_id](2, 1) / FPM_matrix[i->m_id](2, 2) - hxx[i->m_id] * FPM_matrix[i->m_id](2, 3) / FPM_matrix[i->m_id](2, 2) - hyy[i->m_id] * FPM_matrix[i->m_id](2, 4) / FPM_matrix[i->m_id](2, 2) - hxy[i->m_id] * FPM_matrix[i->m_id](2, 5) / FPM_matrix[i->m_id](2, 2);
					hxx[i->m_id] = hj_i_d2Wdx2[i->m_id] * (1.0 / FPM_matrix[i->m_id](3, 3) - 1.0) - hx[i->m_id] * FPM_matrix[i->m_id](3, 1) / FPM_matrix[i->m_id](3, 3) - hy[i->m_id] * FPM_matrix[i->m_id](3, 2) / FPM_matrix[i->m_id](3, 3) - hyy[i->m_id] * FPM_matrix[i->m_id](3, 4) / FPM_matrix[i->m_id](3, 3) - hxy[i->m_id] * FPM_matrix[i->m_id](3, 5) / FPM_matrix[i->m_id](3, 3);
					hyy[i->m_id] = hj_i_d2Wdy2[i->m_id] * (1.0 / FPM_matrix[i->m_id](4, 4) - 1.0) - hx[i->m_id] * FPM_matrix[i->m_id](4, 1) / FPM_matrix[i->m_id](4, 4) - hy[i->m_id] * FPM_matrix[i->m_id](4, 2) / FPM_matrix[i->m_id](4, 4) - hxx[i->m_id] * FPM_matrix[i->m_id](4, 3) / FPM_matrix[i->m_id](4, 4) - hxy[i->m_id] * FPM_matrix[i->m_id](4, 5) / FPM_matrix[i->m_id](4, 4);
					hxy[i->m_id] = hj_i_d2Wdxdy[i->m_id] * (1.0 / FPM_matrix[i->m_id](5, 5) - 1.0) - hx[i->m_id] * FPM_matrix[i->m_id](5, 1) / FPM_matrix[i->m_id](5, 5) - hy[i->m_id] * FPM_matrix[i->m_id](5, 2) / FPM_matrix[i->m_id](5, 5) - hxx[i->m_id] * FPM_matrix[i->m_id](5, 3) / FPM_matrix[i->m_id](5, 5) - hyy[i->m_id] * FPM_matrix[i->m_id](5, 4) / FPM_matrix[i->m_id](5, 5);





					//std::cout << h[i->m_id] << " " << hx[i->m_id] << " " << hy[i->m_id] << " " << hxx[i->m_id] << " " << hyy[i->m_id] << " " << hxy[i->m_id] << "\n";


					aveIntegrals += FPM_matrix[i->m_id];

					aveh += h[i->m_id];
					avehx += hx[i->m_id];
					avehy += hy[i->m_id];
					avehxx += hxx[i->m_id];
					avehyy += hyy[i->m_id];
					avehxy += hxy[i->m_id];



				}

				drho[i->m_id] = drho[i->m_id] / FPM_matrix[i->m_id](0, 0);

				if (i->m_type == REAL) {
					//std::cout << i->m_density.val - drho[i->m_id] << "   ";
					i->m_density.val = drho[i->m_id];
					for (int d = 0;d < nrOfDim;d++) {
						i->m_velocity.val[d] = dv[i->m_id][d];
					}

					avedrho += drho[i->m_id];

					//std::cout << dh[i->m_id] << "/" << W_denominator[i->m_id] << "  ";
					dh[i->m_id] = dh[i->m_id] * (1 / W_denominator[i->m_id] - 1);
					//dh[i->m_id] -= hj_i_dWdr[i->m_id]/ dr_dWdr[i->m_id]* dr_W[i->m_id] / W_denominator[i->m_id];
					//std::cout << hj_i_dWdr[i->m_id] / dr_dWdr[i->m_id] * dr_W[i->m_id] / W_denominator[i->m_id] << "\n";
					for (int d = 0;d < nrOfDim;d++) {
						//std::cout << FPM[i->m_id].getd(static_cast<DIMENSIONS>(d)) << "  ";
						dh[i->m_id] -= FPM[i->m_id].getd(static_cast<DIMENSIONS>(d)) * dr_W[i->m_id][d] / W_denominator[i->m_id];
					}
					//std::cout << i->m_id << "  " << i->nrOfNeighbours << "  i->m_SmR =" << i->m_SmR << " + dh = " << dh[i->m_id] << "\n";
					//i->m_SmR += dh[i->m_id];
					i->m_SmR += h[i->m_id];

				}
			}
		}

		if (computationalDomain_mode == MODE_CD::CD_DEBUG) {

			std::cout << aveh / m_options.nrOfParticles[REAL] << " " << avehx / m_options.nrOfParticles[REAL] << " " << avehy / m_options.nrOfParticles[REAL] << " " << avehxx / m_options.nrOfParticles[REAL] << " " << avehyy / m_options.nrOfParticles[REAL] << " " << avehxy / m_options.nrOfParticles[REAL] << "\n";


			std::cout << aveIntegrals * (1.0 / m_options.nrOfParticles[REAL]) << "\n";
			std::cout << avedrho / m_options.nrOfParticles[REAL] << "\n";


			//std::cout << " max / sum/nr  " << "| " << "         Kernel         " << "| " << "    Kernel Derivetive   " << "  | " << "   Kernel Second Derivetive " <<  "\n";
			//std::cout << "   *pow(dr,0)  " << "| " << maxI_Wij << " / " << sumI_Wij / m_options.nrOfParticles[REAL] << " | " <<  maxI_dWij << " / " << sumI_dWij / m_options.nrOfParticles[REAL] << " | " << maxI_d2Wij << " / " << sumI_d2Wij / m_options.nrOfParticles[REAL] << "\n";
			//std::cout << "   *pow(dr,1)  " << "| " << maxI_Wij_dr << " / " << sumI_Wij_dr / m_options.nrOfParticles[REAL] << " | " << maxI_dWij_dr << " / " << sumI_dWij_dr / m_options.nrOfParticles[REAL] << " | " << maxI_d2Wij_dr << " / " << sumI_d2Wij_dr / m_options.nrOfParticles[REAL] << "\n";
			//std::cout << "   *pow(dr,2)  " << "| " << maxI_Wij_dr2 << " / " << sumI_Wij_dr2 / m_options.nrOfParticles[REAL] << " | " << maxI_dWij_dr2 << " / " << sumI_dWij_dr2 / m_options.nrOfParticles[REAL] << " | " << maxI_d2Wij_dr2 << " / " << sumI_d2Wij_dr2 / m_options.nrOfParticles[REAL] << "\n";
			//std::cout << "(Kernel Integral - 1) : max " << maxI_Wij << " , sum/nr " << sumI_Wij / part_pairs_count << "   Kernel Derivetive Integral : max " << maxI_dWij << " , sum/nr " << sumI_dWij / part_pairs_count << "   Kernel Second Derivetive Integral : max " << maxI_d2Wij << " , sum/nr " << sumI_d2Wij / part_pairs_count << "\n";
			//std::cout << "Kernel*dr Integral  : max " << maxI_Wij_dr << " , sum/nr " << sumI_Wij_dr / part_pairs_count << "   (Kernel Derivetive*dr Integral - 1) : max " << maxI_dWij_dr << " , sum/nr " << sumI_dWij_dr / part_pairs_count << "   Kernel Second Derivetive*dr Integral : max " << maxI_d2Wij_dr << " , sum/nr " << sumI_d2Wij_dr / part_pairs_count << "\n";
			
			//if(nrOfDim == D1){
			//	std::cout << " average  " << "| " << "         W         " << "| " << "    W',x   " << "  | " << "   W'',xx " << "\n";
			//	std::cout << "   *  1        " << "| " << sumI_Wij / m_options.nrOfParticles[REAL] << " | " << sumI_dWdx / m_options.nrOfParticles[REAL] << " | " << sumI_d2Wdx2 / m_options.nrOfParticles[REAL] << "\n";
			//	std::cout << "   *pow(dx,1)  " << "| " << sumI_Wij_dx / m_options.nrOfParticles[REAL] << " | " << sumI_dWdx_dx / m_options.nrOfParticles[REAL] << " | " << sumI_d2Wdx2_dx / m_options.nrOfParticles[REAL] << "\n";
			//	std::cout << "   *pow(dx,2)  " << "| " << sumI_Wij_dx2 / m_options.nrOfParticles[REAL] << " | " << sumI_dWdx_dx2 / m_options.nrOfParticles[REAL] << " | " << sumI_d2Wdx2_dx2 / m_options.nrOfParticles[REAL] << "\n";
			//}
			//else if (nrOfDim == D2) {
			//	std::cout << " average  " << "| " << "         W         " << "| " << "    W',x   " << "  | " << "   W',y " << "| " << "    W'',xx   " << "  | " << "   W'',yy " << "  | " << "   W'',xy " << "\n";
			//	std::cout << "   *1          " << "| " << sumI_Wij / m_options.nrOfParticles[REAL] << " | " << sumI_dWdx / m_options.nrOfParticles[REAL] << " | " << sumI_dWdy / m_options.nrOfParticles[REAL] << " | " << sumI_d2Wdx2 / m_options.nrOfParticles[REAL] << " | " << sumI_d2Wdy2 / m_options.nrOfParticles[REAL] << " | " << sumI_d2Wdxdy / m_options.nrOfParticles[REAL] << "\n";
			//	std::cout << "   *pow(dx,1)  " << "| " << sumI_Wij_dx / m_options.nrOfParticles[REAL] << " | " << sumI_dWdx_dx / m_options.nrOfParticles[REAL] << " | " << sumI_dWdy_dx / m_options.nrOfParticles[REAL] << " | " << sumI_d2Wdx2_dx / m_options.nrOfParticles[REAL] << " | " << sumI_d2Wdy2_dx / m_options.nrOfParticles[REAL] << " | " << sumI_d2Wdxdy_dx / m_options.nrOfParticles[REAL] << "\n";
			//	std::cout << "   *pow(dy,1)  " << "| " << sumI_Wij_dy / m_options.nrOfParticles[REAL] << " | " << sumI_dWdx_dy / m_options.nrOfParticles[REAL] << " | " << sumI_dWdy_dy / m_options.nrOfParticles[REAL] << " | " << sumI_d2Wdx2_dy / m_options.nrOfParticles[REAL] << " | " << sumI_d2Wdy2_dy / m_options.nrOfParticles[REAL] << " | " << sumI_d2Wdxdy_dy / m_options.nrOfParticles[REAL] << "\n";
			//	std::cout << "   *pow(dx,2)  " << "| " << sumI_Wij_dx2 / m_options.nrOfParticles[REAL] << " | " << sumI_dWdx_dx2 / m_options.nrOfParticles[REAL] << " | " << sumI_dWdy_dx2 / m_options.nrOfParticles[REAL] << " | " << sumI_d2Wdx2_dx2 / m_options.nrOfParticles[REAL] << " | " << sumI_d2Wdy2_dx2 / m_options.nrOfParticles[REAL] << " | " << sumI_d2Wdxdy_dx2 / m_options.nrOfParticles[REAL] << "\n";
			//	std::cout << "   *pow(dy,2)  " << "| " << sumI_Wij_dy2 / m_options.nrOfParticles[REAL] << " | " << sumI_dWdx_dy2 / m_options.nrOfParticles[REAL] << " | " << sumI_dWdy_dy2 / m_options.nrOfParticles[REAL] << " | " << sumI_d2Wdx2_dy2 / m_options.nrOfParticles[REAL] << " | " << sumI_d2Wdy2_dy2 / m_options.nrOfParticles[REAL] << " | " << sumI_d2Wdxdy_dy2 / m_options.nrOfParticles[REAL] << "\n";
			//	std::cout << "   *dx*dy      " << "| " << sumI_Wij_dxdy / m_options.nrOfParticles[REAL] << " | " << sumI_dWdx_dxdy / m_options.nrOfParticles[REAL] << " | " << sumI_dWdy_dxdy / m_options.nrOfParticles[REAL] << " | " << sumI_d2Wdx2_dxdy / m_options.nrOfParticles[REAL] << " | " << sumI_d2Wdy2_dxdy / m_options.nrOfParticles[REAL] << " | " << sumI_d2Wdxdy_dxdy / m_options.nrOfParticles[REAL] << "\n";
			//}
			
		}

	//}

	
}
*/
/*
void SPH_CD::TimeDerivative() {
	if (computationalDomain_mode == MODE_CD::CD_DEBUG) {
		std::cout << "SPH_CD::TimeDerivative\n";
	}

	if (m_options.SPH_algorithm == STANDART_SPH) {
		std::vector<part_prec> drhodt;
		drhodt.resize(SPH.Particles.size(), 0.f);
		std::vector<part_prec_3> dvdt;
		dvdt.resize(SPH.Particles.size());
		std::vector<ParticlePair::dsmthFunc> FunckeydV{ &Particle::dVx ,&Particle::dVy ,&Particle::dVz };
		std::vector<ParticlePair::dsmthFunc> FunckeydX{ &Particle::dx ,&Particle::dy ,&Particle::dz,&Particle::distance };
		part_prec presspartM;
		part_prec presspartNB;

		switch (m_options.densityChangeUsing) {
		case(VAL):
			for (auto*& i : SPH.ParticlePairs) {
				//DENSITY CALCULATION
				if ((PPL_BOTH_PARTICLES_ARE(REAL)) or (PPL_MAIN_PARTICLE_IS(VIRTUAL) or PPL_NEIGHBOUR_PARTICLE_IS(VIRTUAL))) {
					drhodt[i->Mpart_ptr->m_id] += (i->NBpart_ptr->m_mass) * i->W(i->Mpart_ptr);
					drhodt[i->NBpart_ptr->m_id] += (i->Mpart_ptr->m_mass) * i->W(i->NBpart_ptr);
				}
				for (int d = 0;d < nrOfDim;d++) {
					if ((PPL_BOTH_PARTICLES_ARE(REAL)) or (PPL_MAIN_PARTICLE_IS(VIRTUAL) or PPL_NEIGHBOUR_PARTICLE_IS(VIRTUAL))) {
						//VELOCITY CALCULATION PRESSURE PART
						switch (m_options.velocityChangeAlgorithm_pressurePart) {
						case(0): // 0
							presspartM = 0;
							presspartNB = 0;
							break;
						case(1): //(-mj*(pi/rhoi^2+pj/rhoj^2))
							presspartM = -(i->NBpart_ptr->m_mass) / ((i->Mpart_ptr->m_density.val)*(i->NBpart_ptr->m_density.val))*((i->Mpart_ptr->m_pressure.val) - (i->NBpart_ptr->m_pressure.val));
							presspartNB = -(i->Mpart_ptr->m_mass) / ((i->NBpart_ptr->m_density.val)*(i->Mpart_ptr->m_density.val))*((i->NBpart_ptr->m_pressure.val) - (i->Mpart_ptr->m_pressure.val));
							break;
						case(2): //(mj*(pi-pj)/(rhoi*rhoj))
							presspartM = (i->NBpart_ptr->m_mass) *((i->Mpart_ptr->m_pressure.val / pow(i->Mpart_ptr->m_density.val, 2)) + (i->NBpart_ptr->m_pressure.val / pow(i->NBpart_ptr->m_density.val, 2)));
							presspartNB = (i->Mpart_ptr->m_mass) *((i->Mpart_ptr->m_pressure.val / pow(i->Mpart_ptr->m_density.val, 2)) + (i->NBpart_ptr->m_pressure.val / pow(i->NBpart_ptr->m_density.val, 2)));
							break;
						case(3): //(-mj*(pi+pj)/(rhoi*rhoj))
							presspartM = -(i->NBpart_ptr->m_mass) / ((i->Mpart_ptr->m_density.val)*(i->NBpart_ptr->m_density.val))*((i->Mpart_ptr->m_pressure.val) + (i->NBpart_ptr->m_pressure.val));
							presspartNB = -(i->Mpart_ptr->m_mass) / ((i->NBpart_ptr->m_density.val)*(i->Mpart_ptr->m_density.val))*((i->NBpart_ptr->m_pressure.val) + (i->Mpart_ptr->m_pressure.val));
							break;
						default:
							assert("SPH_CD::TimeDerivative::STANDART_SPH::m_options.velocityChangeAlgorithm_pressurePart has undefined value" && 0);
							break;
						}
					}
					dvdt[i->Mpart_ptr->m_id][d] += presspartM * (i->dWd(static_cast<DiffAxis>(d), i->Mpart_ptr));
					dvdt[i->NBpart_ptr->m_id][d] += presspartNB * (i->dWd(static_cast<DiffAxis>(d), i->NBpart_ptr));

					if ((PPL_BOTH_PARTICLES_ARE(REAL)) or (PPL_MAIN_PARTICLE_IS(VIRTUAL) or PPL_NEIGHBOUR_PARTICLE_IS(VIRTUAL))) {
						//VELOCITY CALCULATION VISCOCITY PART
						for (int dd = 0;dd < nrOfDim;dd++) {
							switch (m_options.velocityChangeAlgorithm_viscosityPart) {
							case(0): // 0
								break;
							case(1): //мой вывод (как я понимаю)
								//  как-то записать просто производные второго порядка
								break;
							case(2):
								(dvdt)[i->Mpart_ptr->m_id][d] += i->Mpart_ptr->m_DVisc *(i->NBpart_ptr->m_mass) / ((i->Mpart_ptr->m_density.val)*(i->NBpart_ptr->m_density.val)) * 2.f * static_cast<float>(nrOfDim + 2) * (-i->dvfunc(FunckeydV[dd])) * (-i->dvfunc(FunckeydX[dd])) / i->Mpart_ptr->distance(i->NBpart_ptr) * (i->dWd(static_cast<DiffAxis>(dd), i->Mpart_ptr));
								(dvdt)[i->NBpart_ptr->m_id][d] += i->NBpart_ptr->m_DVisc *(i->Mpart_ptr->m_mass) / ((i->NBpart_ptr->m_density.val)*(i->Mpart_ptr->m_density.val)) * 2.f * static_cast<float>(nrOfDim + 2) * (i->dvfunc(FunckeydV[dd])) * (i->dvfunc(FunckeydX[dd])) / i->NBpart_ptr->distance(i->Mpart_ptr) * (i->dWd(static_cast<DiffAxis>(dd), i->NBpart_ptr));
								break;
							case(3):
								(dvdt)[i->Mpart_ptr->m_id][d] += i->Mpart_ptr->m_DVisc *(i->NBpart_ptr->m_mass) / ((i->Mpart_ptr->m_density.val)*(i->NBpart_ptr->m_density.val)) * (i->dvfunc(FunckeydV[dd])) * (i->d2Wd(static_cast<DiffAxis>(dd), static_cast<DiffAxis>(dd), i->Mpart_ptr));
								(dvdt)[i->NBpart_ptr->m_id][d] += i->NBpart_ptr->m_DVisc *(i->Mpart_ptr->m_mass) / ((i->NBpart_ptr->m_density.val)*(i->Mpart_ptr->m_density.val)) * (-i->dvfunc(FunckeydV[dd])) * (i->d2Wd(static_cast<DiffAxis>(dd), static_cast<DiffAxis>(dd), i->NBpart_ptr));
								break;
							default:
								assert("SPH_CD::TimeDerivative::STANDART_SPH::m_options.velocityChangeAlgorithm_viscosityPart has undefined value" && 0);
								break;
							}
						}
					}
				}
			}
			break;
		case(DVAL_DT):
			for (auto*& i : SPH.ParticlePairs) {
				for (int d = 0;d < nrOfDim;d++) {
					if ((PPL_BOTH_PARTICLES_ARE(REAL)) or (PPL_MAIN_PARTICLE_IS(VIRTUAL) or PPL_NEIGHBOUR_PARTICLE_IS(VIRTUAL))) {
						//DENSITY CALCULATION
						switch (m_options.densityChangeAlgorithm) {
						case(0): //Monaghan                                  IS '-' correct?
							drhodt[i->Mpart_ptr->m_id] += (i->NBpart_ptr->m_mass)*(i->dvfunc(FunckeydV[d]))*(i->dWd(static_cast<DiffAxis>(d), i->Mpart_ptr));
							drhodt[i->NBpart_ptr->m_id] += (i->Mpart_ptr->m_mass)*(-i->dvfunc(FunckeydV[d]))*(i->dWd(static_cast<DiffAxis>(d), i->NBpart_ptr));
							break;
						case(1)://                                                                                                      IS '-' correct?
							drhodt[i->Mpart_ptr->m_id] += (i->Mpart_ptr->m_density.val)*(i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val)*(i->dvfunc(FunckeydV[d]))*(i->dWd(static_cast<DiffAxis>(d), i->Mpart_ptr));
							drhodt[i->NBpart_ptr->m_id] += (i->NBpart_ptr->m_density.val)*(i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val)*(-i->dvfunc(FunckeydV[d]))*(i->dWd(static_cast<DiffAxis>(d), i->NBpart_ptr));
							break;
						default:
							assert("SPH_CD::TimeDerivative::STANDART_SPH::DVAL_DT::m_options.densityChangeAlgorithm has undefined value" && 0);
							break;
						}
					}
					if ((PPL_BOTH_PARTICLES_ARE(REAL)) or (PPL_MAIN_PARTICLE_IS(VIRTUAL) or PPL_NEIGHBOUR_PARTICLE_IS(VIRTUAL))) {
						//VELOCITY CALCULATION PRESSURE PART
						switch (m_options.velocityChangeAlgorithm_pressurePart) {
						case(0): // 0
							presspartM = 0;
							presspartNB = 0;
							break;
						case(1): //(mj*(pi-pj)/(rhoi*rhoj))
							presspartM = -(i->NBpart_ptr->m_mass) / ((i->Mpart_ptr->m_density.val)*(i->NBpart_ptr->m_density.val))*((i->Mpart_ptr->m_pressure.val) - (i->NBpart_ptr->m_pressure.val));
							presspartNB = -(i->Mpart_ptr->m_mass) / ((i->NBpart_ptr->m_density.val)*(i->Mpart_ptr->m_density.val))*((i->NBpart_ptr->m_pressure.val) - (i->Mpart_ptr->m_pressure.val));
							break;
						case(2): //(-mj*(pi/rhoi^2+pj/rhoj^2))
							presspartM = (i->NBpart_ptr->m_mass) *((i->Mpart_ptr->m_pressure.val / pow(i->Mpart_ptr->m_density.val, 2)) + (i->NBpart_ptr->m_pressure.val / pow(i->NBpart_ptr->m_density.val, 2)));
							presspartNB = (i->Mpart_ptr->m_mass) *((i->Mpart_ptr->m_pressure.val / pow(i->Mpart_ptr->m_density.val, 2)) + (i->NBpart_ptr->m_pressure.val / pow(i->NBpart_ptr->m_density.val, 2)));
							break;
						case(3): //(-mj*(pi+pj)/(rhoi*rhoj))
							presspartM = -(i->NBpart_ptr->m_mass) / ((i->Mpart_ptr->m_density.val)*(i->NBpart_ptr->m_density.val))*((i->Mpart_ptr->m_pressure.val) + (i->NBpart_ptr->m_pressure.val));
							presspartNB = -(i->Mpart_ptr->m_mass) / ((i->NBpart_ptr->m_density.val)*(i->Mpart_ptr->m_density.val))*((i->NBpart_ptr->m_pressure.val) + (i->Mpart_ptr->m_pressure.val));
							break;
						default:
							assert("SPH_CD::TimeDerivative::STANDART_SPH::m_options.velocityChangeAlgorithm_pressurePart has undefined value" && 0);
							break;
						}
					}
					dvdt[i->Mpart_ptr->m_id][d] += presspartM * (i->dWd(static_cast<DiffAxis>(d), i->Mpart_ptr));
					dvdt[i->NBpart_ptr->m_id][d] += presspartNB * (i->dWd(static_cast<DiffAxis>(d), i->NBpart_ptr));

					if ((PPL_BOTH_PARTICLES_ARE(REAL)) or (PPL_MAIN_PARTICLE_IS(VIRTUAL) or PPL_NEIGHBOUR_PARTICLE_IS(VIRTUAL))) {
						//VELOCITY CALCULATION VISCOCITY PART
						for (int dd = 0;dd < nrOfDim;dd++) {
							switch (m_options.velocityChangeAlgorithm_viscosityPart) {
							case(0): // 0
								break;
							case(1): //мой вывод (как я понимаю)
								//  как-то записать просто производные второго порядка
								break;
							case(2):
								(dvdt)[i->Mpart_ptr->m_id][d] += i->Mpart_ptr->m_DVisc *(i->NBpart_ptr->m_mass) / ((i->Mpart_ptr->m_density.val)*(i->NBpart_ptr->m_density.val)) * 2.f * static_cast<float>(nrOfDim + 2) * (-i->dvfunc(FunckeydV[dd])) * (-i->dvfunc(FunckeydX[dd])) / i->Mpart_ptr->distance(i->NBpart_ptr) * (i->dWd(static_cast<DiffAxis>(dd), i->Mpart_ptr));
								(dvdt)[i->NBpart_ptr->m_id][d] += i->NBpart_ptr->m_DVisc *(i->Mpart_ptr->m_mass) / ((i->NBpart_ptr->m_density.val)*(i->Mpart_ptr->m_density.val)) * 2.f * static_cast<float>(nrOfDim + 2) * (i->dvfunc(FunckeydV[dd])) * (i->dvfunc(FunckeydX[dd])) / i->NBpart_ptr->distance(i->Mpart_ptr) * (i->dWd(static_cast<DiffAxis>(dd), i->NBpart_ptr));
								break;
							case(3):
								(dvdt)[i->Mpart_ptr->m_id][d] += i->Mpart_ptr->m_DVisc *(i->NBpart_ptr->m_mass) / ((i->Mpart_ptr->m_density.val)*(i->NBpart_ptr->m_density.val)) * (i->dvfunc(FunckeydV[dd])) * (i->d2Wd(static_cast<DiffAxis>(dd), static_cast<DiffAxis>(dd), i->Mpart_ptr));
								(dvdt)[i->NBpart_ptr->m_id][d] += i->NBpart_ptr->m_DVisc *(i->Mpart_ptr->m_mass) / ((i->NBpart_ptr->m_density.val)*(i->Mpart_ptr->m_density.val)) * (-i->dvfunc(FunckeydV[dd])) * (i->d2Wd(static_cast<DiffAxis>(dd), static_cast<DiffAxis>(dd), i->NBpart_ptr));
								break;
							default:
								assert("SPH_CD::TimeDerivative::STANDART_SPH::m_options.velocityChangeAlgorithm_viscosityPart has undefined value" && 0);
								break;
							}

						}
					}
				}
			}
			break;
		default:
			assert("SPH_CD::TimeDerivative::STANDART_SPH::m_options.densityChangeUsing has undefined value" && 0);
			break;
		}

		for (auto*& i : SPH.Particles) {
			i->m_density.dval = drhodt[i->m_id];
			for (int d = 0;d < nrOfDim;d++) {
				i->m_velocity.dval[d] = (dvdt)[i->m_id][d];
			}
		}



		// id dv drho

	}
	if (m_options.SPH_algorithm == CORRECTIVE_SPH) {
		std::vector<part_prec> drhodt;
		std::vector<part_prec> W_denominator;

		std::vector<part_prec_3> drhodt_numerator;
		std::vector<part_prec_3> dW_denominator;


		std::vector<part_prec_3> dp_numerator;


		std::vector<part_prec> VelDiv;

		std::vector<part_prec_3> dvdt;
		dvdt.resize(SPH.Particles.size());

		std::vector<ParticlePair::dsmthFunc> FunckeydV{ &Particle::dVx ,&Particle::dVy ,&Particle::dVz };
		std::vector<ParticlePair::dsmthFunc> FunckeydX{ &Particle::dx ,&Particle::dy ,&Particle::dz ,&Particle::distance };




		switch (m_options.densityChangeUsing) {
		case(VAL):
			drhodt.resize(SPH.Particles.size());
			W_denominator.resize(SPH.Particles.size());
			if (m_options.velocityChangeAlgorithm_pressurePart != 0) {
				dp_numerator.resize(SPH.Particles.size());
				dW_denominator.resize(SPH.Particles.size());
			}

			for (auto*& i : SPH.ParticlePairs) {
				if ((PPL_BOTH_PARTICLES_ARE(REAL)) or PPL_ONE_OF_PARTICLE_IS_1_OTHER_IS_2(REAL, VIRTUAL)) {
					//DENSITY CALCULATION
					drhodt[i->Mpart_ptr->m_id] += (i->NBpart_ptr->m_mass) * i->W(i->Mpart_ptr);
					drhodt[i->NBpart_ptr->m_id] += (i->Mpart_ptr->m_mass) * i->W(i->NBpart_ptr);
					W_denominator[i->Mpart_ptr->m_id] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val) * i->W(i->Mpart_ptr);
					W_denominator[i->NBpart_ptr->m_id] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val) * i->W(i->NBpart_ptr);
				}

				for (int d = 0;d < nrOfDim;d++) {
					if ((PPL_BOTH_PARTICLES_ARE(REAL)) or PPL_ONE_OF_PARTICLE_IS_1_OTHER_IS_2(REAL, VIRTUAL)) {
						//VELOCITY CALCULATION PRESSURE PART
						switch (m_options.velocityChangeAlgorithm_pressurePart) {
						case(0): // 0
							break;
						case(1): //(-mj*(pi/rhoi^2+pj/rhoj^2))
							dp_numerator[i->Mpart_ptr->m_id][d] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val)*((i->NBpart_ptr->m_pressure.val) - (i->Mpart_ptr->m_pressure.val))*(i->dWd(static_cast<DiffAxis>(d), i->Mpart_ptr));
							dp_numerator[i->Mpart_ptr->m_id][d] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val)*((i->Mpart_ptr->m_pressure.val) - (i->NBpart_ptr->m_pressure.val))* (i->dWd(static_cast<DiffAxis>(d), i->NBpart_ptr));
							dW_denominator[i->Mpart_ptr->m_id][d] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val)*(-(i->dvfunc(FunckeydX[d])))*(i->dWd(static_cast<DiffAxis>(d), i->Mpart_ptr));
							dW_denominator[i->NBpart_ptr->m_id][d] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val)*((i->dvfunc(FunckeydX[d])))*(i->dWd(static_cast<DiffAxis>(d), i->NBpart_ptr));
							//dW_denominator[i->Mpart_ptr->m_id][d] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val)*(-(i->dvfunc(FunckeydX[3])))*(i->dWd(static_cast<DiffAxis>(d), i->Mpart_ptr));
							//dW_denominator[i->NBpart_ptr->m_id][d] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val)*((i->dvfunc(FunckeydX[3])))*(i->dWd(static_cast<DiffAxis>(d), i->NBpart_ptr));
							break;
						default:
							assert("SPH_CD::TimeDerivative::STANDART_SPH::m_options.velocityChangeAlgorithm_pressurePart has undefined value" && 0);
							break;
						}
					}
				}
			}


			break;
		case(DVAL_DT): {
			dW_denominator.resize(SPH.Particles.size());
			drhodt_numerator.resize(SPH.Particles.size());

			for (auto*& i : SPH.ParticlePairs) {
				for (int d = 0;d < nrOfDim;d++) {
					//DENSITY CALCULATION
					if ((PPL_BOTH_PARTICLES_ARE(REAL)) or PPL_ONE_OF_PARTICLE_IS_1_OTHER_IS_2(REAL, VIRTUAL)) {
						drhodt_numerator[i->Mpart_ptr->m_id][d] += (i->Mpart_ptr->m_density.val)*(i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val)*(-(i->dvfunc(FunckeydV[d])))*(i->dWd(static_cast<DiffAxis>(d), i->Mpart_ptr));
						drhodt_numerator[i->NBpart_ptr->m_id][d] += (i->NBpart_ptr->m_density.val)*(i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val)*((i->dvfunc(FunckeydV[d])))*(i->dWd(static_cast<DiffAxis>(d), i->NBpart_ptr));
						dW_denominator[i->Mpart_ptr->m_id][d] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val)*(-(i->dvfunc(FunckeydX[d])))*(i->dWd(static_cast<DiffAxis>(d), i->Mpart_ptr));
						dW_denominator[i->NBpart_ptr->m_id][d] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val)*((i->dvfunc(FunckeydX[d])))*(i->dWd(static_cast<DiffAxis>(d), i->NBpart_ptr));
					}
				}

				for (int d = 0;d < nrOfDim;d++) {
					if ((PPL_BOTH_PARTICLES_ARE(REAL)) or PPL_ONE_OF_PARTICLE_IS_1_OTHER_IS_2(REAL, VIRTUAL)) {
						//VELOCITY CALCULATION PRESSURE PART
						switch (m_options.velocityChangeAlgorithm_pressurePart) {
						case(0): // 0
							break;
						case(1): //(-mj*(pi/rhoi^2+pj/rhoj^2))
							dp_numerator[i->Mpart_ptr->m_id][d] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val)*((i->NBpart_ptr->m_pressure.val) - (i->Mpart_ptr->m_pressure.val))*(i->dWd(static_cast<DiffAxis>(d), i->Mpart_ptr));
							dp_numerator[i->Mpart_ptr->m_id][d] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val)*((i->Mpart_ptr->m_pressure.val) - (i->NBpart_ptr->m_pressure.val))* (i->dWd(static_cast<DiffAxis>(d), i->NBpart_ptr));
							dW_denominator[i->Mpart_ptr->m_id][d] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val)*(-(i->dvfunc(FunckeydX[d])))*(i->dWd(static_cast<DiffAxis>(d), i->Mpart_ptr));
							dW_denominator[i->NBpart_ptr->m_id][d] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val)*((i->dvfunc(FunckeydX[d])))*(i->dWd(static_cast<DiffAxis>(d), i->NBpart_ptr));
							break;
						default:
							assert("SPH_CD::TimeDerivative::STANDART_SPH::m_options.velocityChangeAlgorithm_pressurePart has undefined value" && 0);
							break;
						}
					}
				}


			}
			//VELOCITY CALCULATION ?
			break;
		}
		default:
			assert("SPH_CD::TimeDerivative::CORRECTIVE_SPH::m_options.densityChangeUsing has undefined value" && 0);
			break;
		}


		for (auto*& i : SPH.Particles) {
			if (m_options.velocityChangeAlgorithm_pressurePart != 0) {
				for (int d = 0;d < nrOfDim;d++) {
					dvdt[i->m_id][d] = -dp_numerator[i->m_id][d] / dW_denominator[i->m_id][d] / (i->m_density.val);
				}
			}
		}

		for (auto*& i : SPH.ParticlePairs) {
			for (int d = 0;d < nrOfDim;d++) {
				if ((PPL_BOTH_PARTICLES_ARE(REAL)) or PPL_ONE_OF_PARTICLE_IS_1_OTHER_IS_2(REAL, VIRTUAL)) {
					//VELOCITY CALCULATION VISCOCITY PART
					//for (int dd = 0;dd < nrOfDim;dd++) {
					switch (m_options.velocityChangeAlgorithm_viscosityPart) {
					case(0): // 0
						break;
					case(1): //мой вывод (как я понимаю)
						//  как-то записать просто производные второго порядка
						break;
					case(2):
						if (i->dvfunc(FunckeydX[d]) != 0) {
							(dvdt)[i->Mpart_ptr->m_id][d] += i->Mpart_ptr->m_DVisc *(i->NBpart_ptr->m_mass) / ((i->Mpart_ptr->m_density.val)*(i->NBpart_ptr->m_density.val)) * 2.f * static_cast<float>(nrOfDim + 2) * (-i->dvfunc(FunckeydV[d])) * (-i->dvfunc(FunckeydX[d])) / i->Mpart_ptr->distance(i->NBpart_ptr)  * (i->dWd(static_cast<DiffAxis>(d), i->Mpart_ptr));
							(dvdt)[i->NBpart_ptr->m_id][d] += i->NBpart_ptr->m_DVisc *(i->Mpart_ptr->m_mass) / ((i->NBpart_ptr->m_density.val)*(i->Mpart_ptr->m_density.val)) * 2.f * static_cast<float>(nrOfDim + 2) * (i->dvfunc(FunckeydV[d])) * (i->dvfunc(FunckeydX[d])) / i->NBpart_ptr->distance(i->Mpart_ptr)  * (i->dWd(static_cast<DiffAxis>(d), i->NBpart_ptr));
						}
						//std::cout << (-i->dvfunc(FunckeydV[dd])) * (-i->dvfunc(FunckeydX[dd])) << "\n";
					case(3):
						(dvdt)[i->Mpart_ptr->m_id][d] += i->Mpart_ptr->m_DVisc *(i->NBpart_ptr->m_mass) / ((i->Mpart_ptr->m_density.val)*(i->NBpart_ptr->m_density.val)) * (i->dvfunc(FunckeydV[d])) * (i->d2Wd(static_cast<DiffAxis>(d), static_cast<DiffAxis>(d), i->Mpart_ptr));
						(dvdt)[i->NBpart_ptr->m_id][d] += i->NBpart_ptr->m_DVisc *(i->Mpart_ptr->m_mass) / ((i->NBpart_ptr->m_density.val)*(i->Mpart_ptr->m_density.val)) * (-i->dvfunc(FunckeydV[d])) * (i->d2Wd(static_cast<DiffAxis>(d), static_cast<DiffAxis>(d), i->NBpart_ptr));
						break;

						break;
					default:
						assert("SPH_CD::TimeDerivative::STANDART_SPH::m_options.velocityChangeAlgorithm_viscosityPart has undefined value" && 0);
						break;
					}
					//}
				}
			}
		}

		for (auto*& i : SPH.Particles) {
			switch (m_options.densityChangeUsing) {
			case(DVAL_DT):
				VelDiv.resize(SPH.Particles.size(), 0.f);
				for (int d = 0;d < nrOfDim;d++) {
					VelDiv[i->m_id] += drhodt_numerator[i->m_id][d] / dW_denominator[i->m_id][d];
				}
				i->m_density.dval = -i->m_density.val*VelDiv[i->m_id];

				for (int d = 0;d < nrOfDim;d++) {
					i->m_velocity.dval[d] = (dvdt)[i->m_id][d];
				}

				break;
			case(VAL):
				i->m_density.dval = drhodt[i->m_id] / W_denominator[i->m_id];
				for (int d = 0;d < nrOfDim;d++) {
					i->m_velocity.dval[d] = (dvdt)[i->m_id][d];
				}
				break;
			default:
				break;
			}
		}


		// id dv drho
	}
	if (m_options.SPH_algorithm == FINITE_PARTICLE_METHOD) {
		//INITIALIZING
		std::vector<part_prec> drhodt;
		std::vector<part_prec_3> dvdt;
		drhodt.resize(SPH.Particles.size(), 0.f);
		dvdt.resize(SPH.Particles.size(), part_prec_3(0.f));
		std::vector<std::string> Velocities;
		Velocities.push_back("Vx"); Velocities.push_back("Vy"); Velocities.push_back("Vz");
		std::vector<std::string> FPM_names;
		for (int d = 0;d < nrOfDim;d++) { FPM_names.push_back(Velocities[d]); }
		FPM_names.push_back("Pressure");
		std::vector<part_prec> drhodt_vec;
		drhodt_vec.resize(SPH.Particles.size(), 0.f);
		std::vector<part_prec> W_denominator;
		W_denominator.resize(SPH.Particles.size(), 0.f);
		int m_size = 0;
		for (int i = 1;i <= nrOfDim + 1;i++) { m_size += i; }

		std::vector<Matrix> FPM_matrix;
		FPM_matrix.resize(SPH.Particles.size(), Matrix()); // на удаление
		for (int i = 0;i < SPH.Particles.size();i++) { FPM_matrix[i].resize(m_size, m_size); }

		std::vector<std::vector<std::vector<part_prec>>> FPM_vecs;
		FPM_vecs.resize(nrOfDim + 1); //nrOfParemetersToCalculate
		for (auto &vec : FPM_vecs) {
			vec.resize(SPH.Particles.size(), std::vector<part_prec>()); // на удаление
			for (auto &jj : vec) {
				jj.resize(m_size, 0.f); // на удаление
			}
		}
		std::vector<std::vector<std::vector<part_prec>>> FPM_results;
		FPM_results.resize(nrOfDim + 1); //nrOfParemetersToCalculate
		for (auto &res : FPM_results) {
			res.resize(SPH.Particles.size(), std::vector<part_prec>()); // на удаление
			for (auto &jj : res) {
				jj.resize(m_size, 0.f);
			}
		}
		std::vector<std::vector<FiniteParticleMethod>> FPM;
		FPM.resize(nrOfDim + 1); //nrOfParemetersToCalculate
		//for (auto i : FPM) {
		//	i.resize(SPH.Particles.size(), FiniteParticleMethod(nrOfDim));
		//}

		//int counter=0;

		std::cout << "SPH_CD::TimeDerivative::staring calculation\n";
		for (ParticlePair* i : SPH.ParticlePairs) {

			if (PPL_BOTH_PARTICLES_ARE(REAL) or PPL_ONE_OF_PARTICLE_IS_1_OTHER_IS_2(REAL, VIRTUAL) or PPL_ONE_OF_PARTICLE_IS_1_OTHER_IS_2(REAL, BOUNDARY)) {
				//if ((i->Mpart_ptr->m_type == REAL) and (i->NBpart_ptr->m_type == REAL)) {
					//if((i->Mpart_ptr->m_id == 430)or(i->NBpart_ptr->m_id == 430)){
					//	std::cout << FPM_matrix[430](3,0) << " ->" ;
					//}
					//  АЛГОРИТМЫ
				PPC_FPMStandart_General(i, &FPM_matrix);
				//PPC_FPMDifference_General(i, &FPM_matrix);
				//std::cout << counter++ << " ";
				for (int vec = 0;vec < FPM_vecs.size();vec++) {
					PPC_FPMStandart_vec(i, &FPM_vecs[vec], FPM_names[vec]);
					//PPC_FPMDifference_vec(i, &FPM_vecs[vec], FPM_names[vec]);
				}

				switch (m_options.densityChangeUsing) {
				case(VAL):
					drhodt_vec[i->Mpart_ptr->m_id] += (i->NBpart_ptr->m_mass) * i->W(i->Mpart_ptr);
					drhodt_vec[i->NBpart_ptr->m_id] += (i->Mpart_ptr->m_mass) * i->W(i->NBpart_ptr);
					//std::cout << i->Mpart_ptr->m_id << "  mass: " << i->Mpart_ptr->m_mass << "\n";
					W_denominator[i->Mpart_ptr->m_id] += (i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val) * i->W(i->Mpart_ptr);
					W_denominator[i->NBpart_ptr->m_id] += (i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val) * i->W(i->NBpart_ptr);


					break;
				default:
					break;
				}
			}
		}
		//std::cout << "\n";
		//if (computationalDomain_mode == MODE_CD::CD_DEBUG) {
		//	std::cout << "	id: " << SPH.Particles[0]->m_id << " NrOfNeighbours: " << SPH.Particles[0]->nrOfNeighbours << "\n";
		//	std::cout << FPM_matrix[SPH.Particles[0]->m_id] << "\n";
		//	std::cout << "	id: " << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] / 2]->m_id << " NrOfNeighbours: " << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] / 2]->nrOfNeighbours << "\n";
		//	std::cout << FPM_matrix[SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] / 2]->m_id] << "\n";
		//	std::cout << "	id: " << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1]->m_id << " NrOfNeighbours: " << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1]->nrOfNeighbours << "\n";
		//	std::cout << FPM_matrix[SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1]->m_id] << "\n";
		//}
		for (auto*& i : SPH.Particles) {
			if ((i->m_type == PARTICLETYPE::REAL) or (i->m_type == PARTICLETYPE::BOUNDARY)) {

				//dvx_di[i->m_id] = FPM[i->m_id].Results(FPM_matrix[i->m_id], FPM_Vx_vec[i->m_id]);
				//dvy_di[i->m_id] = FPM[i->m_id].Results(FPM_matrix[i->m_id], FPM_Vy_vec[i->m_id]);
				//dvz_di[i->m_id] = FPM[i->m_id].Results(FPM_matrix[i->m_id], FPM_Vz_vec[i->m_id]);
				//dP_di[i->m_id] = FPM[i->m_id].Results(FPM_matrix[i->m_id], FPM_pressure_vec[i->m_id]);

				//std::cout << i->m_id << " MASS: " << i->m_mass  << " NrOfNeighbours: "<<i->nrOfNeighbours  << "\n";
				//std::cout << FPM_matrix[i->m_id] << "\n";
				
				std::cout << i->m_id << "\n";
				std::cout << FPM_matrix[i->m_id];
				for (int name = 0;name < FPM_names.size();name++) {
					std::cout << FPM_names[name] << "       ";

				}
				std::cout << "\n";
				for (int vecs = 0; vecs < FPM_vecs[0][0].size();vecs++) {
					for (int name = 0;name < FPM_names.size();name++) {
						std::cout << FPM_vecs[name][i->m_id][vecs] << "   ";
					}
					std::cout << "\n";
				}
				std::cout << determinant(FPM_matrix[i->m_id]) << "\n";
				



				for (int fpm = 0; fpm < FPM.size();fpm++) {
					FPM[fpm].push_back(FiniteParticleMethod(nrOfDim));
					FPM[fpm][i->m_id].addMatVec(FPM_matrix[i->m_id], FPM_vecs[fpm][i->m_id]);
					FPM[fpm][i->m_id].Calculate();
				}




				part_prec divV = 0.f;
				for (int d = 0;d < nrOfDim;d++) {
					divV += FPM[d][i->m_id].getd(static_cast<DIMENSIONS>(d + 1));
				}
				switch (m_options.densityChangeUsing) {
				case(DVAL_DT):
					i->m_density.dval = -i->m_density.val*(divV);
					break;
				case(VAL):
					i->m_density.dval = drhodt_vec[i->m_id] / W_denominator[i->m_id];
					break;
				default:
					break;
				}



				for (int d = 0;d < nrOfDim;d++) {

					//std::cout << "dv[" << d << "] = -1/rho*(dp/d(" << d << "))";

					switch (m_options.velocityChangeAlgorithm_pressurePart) {
					case(0): // 0
						i->m_velocity.dval[d] = 0;
						break;
					case(1):
						i->m_velocity.dval[d] = -1.f / (i->m_density.val)*FPM[nrOfDim][i->m_id].getd(static_cast<DIMENSIONS>(d + 1));
						break;
					default:
						assert("SPH_CD::TimeDerivative::STANDART_SPH::m_options.velocityChangeAlgorithm_pressurePart has undefined value" && 0);
						break;
					}

					for (int dd = 0;dd < nrOfDim;dd++) {
						switch (m_options.velocityChangeAlgorithm_viscosityPart) {
						case(0): // 0
							i->m_velocity.dval[d] += 0;
							break;
						case(1):
							i->m_velocity.dval[d] += 1.f / (i->m_density.val)*i->m_DVisc*(FPM[dd][i->m_id].getd(static_cast<DIMENSIONS>(d + 1), static_cast<DIMENSIONS>(dd + 1)) + FPM[d][i->m_id].getd(static_cast<DIMENSIONS>(dd + 1), static_cast<DIMENSIONS>(dd + 1))) - 1.f / (i->m_density.val)*i->m_DVisc*2.f / 3.f*FPM[dd][i->m_id].getd(static_cast<DIMENSIONS>(d + 1), static_cast<DIMENSIONS>(dd + 1));
							break;
						default:
							assert("SPH_CD::TimeDerivative::STANDART_SPH::m_options.velocityChangeAlgorithm_viscosityPart has undefined value" && 0);
							break;
						}//std::cout << "+ 1/rho*Visc*(d2V[" << dd << "]/(d" << d << "d" << dd << ")+ d2V[" << d << "]/(d" << dd << "d" << dd << ")" << "- 1/rho*Visc*2/3*(d2V[" << dd << "]/(d" << d << "d" << dd << ")";
					}
					//std::cout << "\n";




				}
			}


			//Частицы на краях надо как-то по-особому просчитывать
			
			if (FPM_matrix[i->m_id](0,0) < 1.f) {
				i->m_density.dval = 0.f;
				for (int d = 0;d < nrOfDim;d++) {
					i->m_velocity.dval[d] = 0.f;
				}
			}
			



		}
		//std::cout << i->m_id << " drho " << i->m_density.dval << " dvx " << i->m_velocity.dval[0] << " dvy " << i->m_velocity.dval[1] << "\n";

		// id dv drho
			

	}

	if (m_options.SPH_algorithm == THE_ONLY_RIGHT) {
		std::vector<part_prec> drhodt;
		drhodt.resize(SPH.Particles.size(), 0.f);
		std::vector<part_prec_3> dvdt;
		dvdt.resize(SPH.Particles.size());
		std::vector<ParticlePair::dsmthFunc> FunckeydV{ &Particle::dVx ,&Particle::dVy ,&Particle::dVz };
		std::vector<ParticlePair::dsmthFunc> FunckeydX{ &Particle::dx ,&Particle::dy ,&Particle::dz ,&Particle::distance };

		for (auto* i : SPH.ParticlePairs) {

			//if ((PPL_BOTH_PARTICLES_ARE(REAL)) or (PPL_MAIN_PARTICLE_IS(VIRTUAL) or PPL_NEIGHBOUR_PARTICLE_IS(VIRTUAL))) {
			//if(PPL_BOTH_PARTICLES_ARE_NOT(BOUNDARY)){
				DensityCalculation(&drhodt, i);
				PressurePartCalculation(&dvdt, i);
				ViscosityPartCalculation(&dvdt, i);
			
			
			//}


		}
		for (auto*& i : SPH.Particles) {
			//if (drhodt[i->m_id] != 0)
			//	std::cout << "	drho:		{" << drhodt[i->m_id] << "}\n";
			i->m_density.dval = drhodt[i->m_id];
			for (int d = 0;d < nrOfDim;d++) {
				i->m_velocity.dval[d] = (dvdt)[i->m_id][d];
			}
		}
	}

	if (computationalDomain_mode == MODE_CD::CD_DEBUG) {
		std::cout << "	id:		first:" << SPH.Particles[0]->m_id << "	";
		std::cout << "	middle:" << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] / 2]->m_id << "	";
		std::cout << "last:" << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1]->m_id << "\n";


		std::cout << "	dv:		{" << SPH.Particles[0]->m_velocity.dval.x << "," << SPH.Particles[0]->m_velocity.dval.y << "," << SPH.Particles[0]->m_velocity.dval.z << "} ";
		std::cout << "	{" << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] / 2]->m_velocity.dval.x << "," << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] / 2]->m_velocity.dval.y << "," << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] / 2]->m_velocity.dval.z << "} ";
		std::cout << "	{" << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1]->m_velocity.dval.x << "," << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1]->m_velocity.dval.y << "," << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1]->m_velocity.dval.z << "}\n";
		switch (m_options.densityChangeUsing) {
		case(DVAL_DT):
			std::cout << "	drho:	";
			break;
		case(VAL):
			std::cout << "	new_density:     ";
			break;
		}
		std::cout << SPH.Particles[0]->m_density.dval;
		std::cout << "	" << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] / 2]->m_density.dval << "	";
		std::cout << "	" << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1]->m_density.dval << "\n";

	}


}
*/

//  LOGIC CUALITY OF LIFE (SPH.ParticlePairs-loop)
// MAIN PARTICLE IS 
#define PPL_MAIN_PARTICLE_IS(type,pp) (pp->Mpart_ptr->m_type == PARTICLETYPE::type)
#define PPL_MAIN_PARTICLE_IS_NOT(type,pp) (pp->Mpart_ptr->m_type != PARTICLETYPE::type)
	// NEIGBOUR PARTICLE IS 
#define PPL_NEIGHBOUR_PARTICLE_IS(type,pp) (pp->NBpart_ptr->m_type == PARTICLETYPE::type)
#define PPL_NEIGHBOUR_PARTICLE_IS_NOT(type,pp) (pp->NBpart_ptr->m_type != PARTICLETYPE::type)


#define PPL_BOTH_PARTICLES_ARE(type,pp) ((pp->NBpart_ptr->m_type == PARTICLETYPE::type) and (pp->Mpart_ptr->m_type == PARTICLETYPE::type))
#define PPL_BOTH_PARTICLES_ARE_NOT(type,pp) ((pp->NBpart_ptr->m_type != PARTICLETYPE::type) and (pp->Mpart_ptr->m_type != PARTICLETYPE::type))


void SPH_CD::PPC_Density(ParticlePair* pp, std::vector<part_prec>* drhodt) {
	//std::cout << pp << "\n";
	part_prec_3 dV;
	std::vector<ParticlePair::dsmthFunc> Funckey{ &Particle::dVx ,&Particle::dVy ,&Particle::dVz };
	for (int d = 0;d < nrOfDim;d++) {
		dV[d] = pp->dvfunc(Funckey[d]);    //   (Vm - Vnb)
	}

	for (int d = 0;d < nrOfDim;d++) {
		//drhodt[i->Mpart_ptr->m_id] += (i->Mpart_ptr->m_density.val)*((i->NBpart_ptr->m_mass) / (i->NBpart_ptr->m_density.val) * dvxdWx);
		//drhodt[i->NBpart_ptr->m_id] += (i->NBpart_ptr->m_density.val)*((i->Mpart_ptr->m_mass) / (i->Mpart_ptr->m_density.val) * dvxdWx);
		//if ((PPL_BOTH_PARTICLES_ARE(REAL,pp)) or
		//	(PPL_MAIN_PARTICLE_IS(VIRTUAL,pp) or PPL_NEIGHBOUR_PARTICLE_IS(VIRTUAL,pp))) {


			switch (m_options.densityChangeAlgorithm) {
			case(0): //мой вывод (как я понимаю)
				(*drhodt)[pp->Mpart_ptr->m_id] += (pp->NBpart_ptr->m_mass) * (-dV[d]) * (pp->dWd(static_cast<DiffAxis>(d), pp->Mpart_ptr));
				(*drhodt)[pp->NBpart_ptr->m_id] += (pp->Mpart_ptr->m_mass) * (dV[d]) * (pp->dWd(static_cast<DiffAxis>(d), pp->NBpart_ptr));
				break;
			case(1):
				(*drhodt)[pp->Mpart_ptr->m_id] += -(pp->Mpart_ptr->m_density.val)*(pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val) *(-dV[d])*(pp->dWd(static_cast<DiffAxis>(d), pp->Mpart_ptr));
				(*drhodt)[pp->NBpart_ptr->m_id] += -(pp->NBpart_ptr->m_density.val)*(pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val) *(dV[d])*(pp->dWd(static_cast<DiffAxis>(d), pp->NBpart_ptr));
				break;
			case(2):
				(*drhodt)[pp->Mpart_ptr->m_id] += (pp->Mpart_ptr->m_mass) * pp->W(pp->Mpart_ptr);
				(*drhodt)[pp->NBpart_ptr->m_id] += (pp->NBpart_ptr->m_mass) * pp->W(pp->NBpart_ptr);
				break;
			default:
				break;
			}
		//}
	}


}
void SPH_CD::PPC_InternalForces(ParticlePair* pp, std::vector<part_prec_3>* dvdt) {



	part_prec_3 dV;
	part_prec_3 dX;

	dX = { pp->Mpart_ptr->dx(pp->NBpart_ptr), pp->Mpart_ptr->dy(pp->NBpart_ptr), pp->Mpart_ptr->dz(pp->NBpart_ptr) };
	dV = { pp->Mpart_ptr->dVx(pp->NBpart_ptr), pp->Mpart_ptr->dVy(pp->NBpart_ptr), pp->Mpart_ptr->dVz(pp->NBpart_ptr) };

	float presspartM;
	float presspartNB;

	switch (m_options.velocityChangeAlgorithm_pressurePart) {
	case(0): //мой вывод (как я понимаю)
		presspartM = (pp->NBpart_ptr->m_mass) *((pp->Mpart_ptr->m_pressure.val / pow(pp->Mpart_ptr->m_density.val, 2)) + (pp->NBpart_ptr->m_pressure.val / pow(pp->NBpart_ptr->m_density.val, 2)));
		presspartNB = (pp->Mpart_ptr->m_mass) *((pp->Mpart_ptr->m_pressure.val / pow(pp->Mpart_ptr->m_density.val, 2)) + (pp->NBpart_ptr->m_pressure.val / pow(pp->NBpart_ptr->m_density.val, 2)));
		break;
	case(1):
		presspartM = (pp->NBpart_ptr->m_mass) / ((pp->Mpart_ptr->m_density.val)*(pp->NBpart_ptr->m_density.val))*((pp->Mpart_ptr->m_pressure.val) + (pp->NBpart_ptr->m_pressure.val));
		presspartNB = (pp->Mpart_ptr->m_mass) / ((pp->NBpart_ptr->m_density.val)*(pp->Mpart_ptr->m_density.val))*((pp->NBpart_ptr->m_pressure.val) + (pp->Mpart_ptr->m_pressure.val));
		//std::cout << presspartM << "  =  " << (pp->NBpart_ptr->m_mass) <<"/(" << (pp->Mpart_ptr->m_density.val)  << "*" << (pp->NBpart_ptr->m_density.val) <<  ")*(" << (pp->Mpart_ptr->m_pressure.val)  << " + " << (pp->NBpart_ptr->m_pressure.val) << ")\n";
		break;
	default:
		break;
	}

	for (int d = 0;d < nrOfDim;d++) {
		if ((PPL_BOTH_PARTICLES_ARE(REAL, pp)) or
			(PPL_MAIN_PARTICLE_IS(VIRTUAL, pp) or PPL_NEIGHBOUR_PARTICLE_IS(VIRTUAL, pp)) or (PPL_MAIN_PARTICLE_IS(BOUNDARY, pp) or PPL_NEIGHBOUR_PARTICLE_IS(BOUNDARY, pp))) {

			switch (m_options.velocityChangeAlgorithm_pressurePart) {
			case(0): //мой вывод (как я понимаю)
				(*dvdt)[pp->Mpart_ptr->m_id][d] += presspartM * (pp->dWd(static_cast<DiffAxis>(d), pp->Mpart_ptr));
				(*dvdt)[pp->NBpart_ptr->m_id][d] += presspartNB * (pp->dWd(static_cast<DiffAxis>(d), pp->NBpart_ptr));
				break;
			case(1):
				(*dvdt)[pp->Mpart_ptr->m_id][d] += -presspartM * (pp->dWd(static_cast<DiffAxis>(d), pp->Mpart_ptr));
				(*dvdt)[pp->NBpart_ptr->m_id][d] += -presspartNB * (pp->dWd(static_cast<DiffAxis>(d), pp->NBpart_ptr));
				break;
			default:
				break;
			}

			//if((pp->dWd(static_cast<DiffAxis>(d), pp->Mpart_ptr))and(pp->dWd(static_cast<DiffAxis>(d), pp->NBpart_ptr))){
			//	std::cout << pp->Mpart_ptr->m_id << " dWdx" << d << " = " << pp->dWd(static_cast<DiffAxis>(d), pp->Mpart_ptr) << " PRESSPART = "<< presspartM << "\n";
			//	std::cout << "pairs   " << pp->dWd(static_cast<DiffAxis>(d), pp->Mpart_ptr) << "   " << pp->dWd(static_cast<DiffAxis>(d), pp->NBpart_ptr) << "\n";
			//}

			for (int dd = 0;dd < nrOfDim;dd++) {

				switch (m_options.velocityChangeAlgorithm_viscosityPart) {
				case(0): //мой вывод (как я понимаю)
					//(*dvdt)[pp->Mpart_ptr->m_id][d] += -pp->NBpart_ptr->m_mass * ((Tau[pp->Mpart_ptr->m_id][d][dd] / pow(pp->Mpart_ptr->m_density.val, 2)) + (Tau[pp->NBpart_ptr->m_id][d][dd] / pow(pp->NBpart_ptr->m_density.val, 2))) * (pp->dWd(static_cast<DiffAxis>(dd), pp->Mpart_ptr));
					//(*dvdt)[pp->NBpart_ptr->m_id][d] += -pp->Mpart_ptr->m_mass * ((Tau[pp->Mpart_ptr->m_id][d][dd] / pow(pp->Mpart_ptr->m_density.val, 2)) + (Tau[pp->NBpart_ptr->m_id][d][dd] / pow(pp->NBpart_ptr->m_density.val, 2))) * (pp->dWd(static_cast<DiffAxis>(dd), pp->NBpart_ptr));
					break;
				case(1):
					(*dvdt)[pp->Mpart_ptr->m_id][d] += pp->Mpart_ptr->m_DVisc *(pp->NBpart_ptr->m_mass) / ((pp->Mpart_ptr->m_density.val)*(pp->NBpart_ptr->m_density.val)) * 2.f * static_cast<float>(nrOfDim + 2) * (-dV[dd]) * (-dX[dd]) / pp->Mpart_ptr->distance(pp->NBpart_ptr) * (pp->dWd(static_cast<DiffAxis>(dd), pp->Mpart_ptr));
					(*dvdt)[pp->NBpart_ptr->m_id][d] += pp->NBpart_ptr->m_DVisc *(pp->Mpart_ptr->m_mass) / ((pp->NBpart_ptr->m_density.val)*(pp->Mpart_ptr->m_density.val)) * 2.f * static_cast<float>(nrOfDim + 2) * (dV[dd]) * (dX[dd]) / pp->NBpart_ptr->distance(pp->Mpart_ptr) * (pp->dWd(static_cast<DiffAxis>(dd), pp->NBpart_ptr));
					break;
				default:
					break;
				}

			}

		}


	}
}
/*
void SPH_CD::PPC_DVx_FPM(ParticlePair* pp, std::vector<FiniteParticleMethod>* FPM_matrix) {

	glm::vec4 vecM = glm::vec4(0.f);     glm::vec4 vecNB = glm::vec4(0.f);
	glm::mat4x4 matM = glm::mat4x4(0.f);   glm::mat4x4 matNB = glm::mat4x4(0.f);
	vecM[0] = pp->NBpart_ptr->m_velocity.val.x *pp->W(pp->Mpart_ptr)*(pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
	vecNB[0] = pp->Mpart_ptr->m_velocity.val.x *pp->W(pp->NBpart_ptr)*(pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);

	glm::vec4 tempVecM;  glm::vec4 tempVecNB;

	glm::vec3 Xji = { pp->NBpart_ptr->dx(pp->Mpart_ptr),pp->NBpart_ptr->dy(pp->Mpart_ptr),pp->NBpart_ptr->dz(pp->Mpart_ptr) };
	glm::vec3 Xij = { pp->Mpart_ptr->dx(pp->NBpart_ptr),pp->Mpart_ptr->dy(pp->NBpart_ptr),pp->Mpart_ptr->dz(pp->NBpart_ptr) };

	tempVecM[0] = pp->W(pp->Mpart_ptr)*(pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
	tempVecNB[0] = pp->W(pp->Mpart_ptr)*(pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);

	for (int dd = 0;dd < nrOfDim;dd++) {
		tempVecM[dd + 1] = Xji[dd] * pp->W(pp->Mpart_ptr)*(pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
		tempVecNB[dd + 1] = Xij[dd] * pp->W(pp->NBpart_ptr)*(pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
	}

	matM[0] = tempVecM;
	matNB[0] = tempVecNB;

	for (int d = 0;d < nrOfDim;d++) {


		vecM[d + 1] = pp->NBpart_ptr->m_velocity.val.x *(pp->dWd(static_cast<DiffAxis>(d), pp->Mpart_ptr))*(pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
		vecNB[d + 1] = pp->Mpart_ptr->m_velocity.val.x *(pp->dWd(static_cast<DiffAxis>(d), pp->NBpart_ptr))*(pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);


		tempVecM[0] = (pp->dWd(static_cast<DiffAxis>(d), pp->Mpart_ptr))*(pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
		tempVecNB[0] = (pp->dWd(static_cast<DiffAxis>(d), pp->NBpart_ptr))*(pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);

		for (int dd = 0;dd < nrOfDim;dd++) {
			tempVecM[dd + 1] = Xji[dd] * (pp->dWd(static_cast<DiffAxis>(d), pp->Mpart_ptr))*(pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
			tempVecNB[dd + 1] = Xij[dd] * (pp->dWd(static_cast<DiffAxis>(d), pp->NBpart_ptr))*(pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
		}

		matM[d + 1] = tempVecM;
		matNB[d + 1] = tempVecNB;
	}

}
*/
void SPH_CD::PPC_FPMStandart_General(ParticlePair* pp, std::vector<Matrix>* FPM_matrix) {


	int m_size = (*FPM_matrix)[0].getRows();
	

	Matrix matM;                    Matrix matNB;
	matM.resize(m_size, m_size);    matNB.resize(m_size, m_size);

	std::vector<float> tempVecM;     std::vector<float> tempVecNB;
	tempVecM.resize(m_size,0.f); tempVecNB.resize(m_size,0.f);

	glm::vec3 Xji = { pp->NBpart_ptr->dx(pp->Mpart_ptr),pp->NBpart_ptr->dy(pp->Mpart_ptr),pp->NBpart_ptr->dz(pp->Mpart_ptr) };
	glm::vec3 Xij = { pp->Mpart_ptr->dx(pp->NBpart_ptr),pp->Mpart_ptr->dy(pp->NBpart_ptr),pp->Mpart_ptr->dz(pp->NBpart_ptr) };

	tempVecM[0] = pp->W(pp->Mpart_ptr)*(pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
	tempVecNB[0] = pp->W(pp->NBpart_ptr)*(pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
	

	for (int dd = 0;dd < nrOfDim;dd++) {
		tempVecM[dd + 1] = Xji[dd] * pp->W(pp->Mpart_ptr)*(pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
		tempVecNB[dd + 1] = Xij[dd] * pp->W(pp->NBpart_ptr)*(pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);

		tempVecM[dd + nrOfDim + 1] = pow(Xji[dd],2)/2.f * pp->W(pp->Mpart_ptr)*(pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
		tempVecNB[dd + nrOfDim + 1] = pow(Xij[dd],2)/2.f * pp->W(pp->NBpart_ptr)*(pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
	}

	if (nrOfDim == 2) {
		tempVecM[2 * nrOfDim + 1] = Xji[X]* Xji[Y] * pp->W(pp->Mpart_ptr)*(pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
		tempVecNB[2 * nrOfDim + 1] = Xij[X] * Xij[Y] * pp->W(pp->NBpart_ptr)*(pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
	}
	else if (nrOfDim == 3) {
		tempVecM[2 * nrOfDim + 1] = Xji[X] * Xji[Y] * pp->W(pp->Mpart_ptr)*(pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
		tempVecNB[2 * nrOfDim + 1] = Xij[X] * Xij[Y] * pp->W(pp->NBpart_ptr)*(pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);

		tempVecM[2 * nrOfDim + 2] = Xji[X] * Xji[Z] * pp->W(pp->Mpart_ptr)*(pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
		tempVecNB[2 * nrOfDim + 2] = Xij[X] * Xij[Z] * pp->W(pp->NBpart_ptr)*(pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);

		tempVecM[2 * nrOfDim + 3] = Xji[Y] * Xji[Z] * pp->W(pp->Mpart_ptr)*(pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
		tempVecNB[2 * nrOfDim + 3] = Xij[Y] * Xij[Z] * pp->W(pp->NBpart_ptr)*(pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
	}

	
	for (int c = 0;c < m_size;c++) {


		matM(0,c) = tempVecM[c];
		matNB(0,c) = tempVecNB[c];
	}




	float derivativeM;
	float derivativeNB;

	for (int r = 1;r < m_size;r++) {
		//Определяем порядок производной
		switch (nrOfDim) {
		case(1):
			switch (r) {
			case(1):
				derivativeM = (pp->dWd(static_cast<DiffAxis>(r - 1), pp->Mpart_ptr));
				derivativeNB = (pp->dWd(static_cast<DiffAxis>(r - 1), pp->NBpart_ptr));
				break;
			case(2):
				derivativeM = (pp->d2Wd(static_cast<DiffAxis>(r - nrOfDim - 1), static_cast<DiffAxis>(r - nrOfDim - 1), pp->Mpart_ptr));
				derivativeNB = (pp->d2Wd(static_cast<DiffAxis>(r - nrOfDim - 1), static_cast<DiffAxis>(r - nrOfDim - 1), pp->NBpart_ptr));
				break;
			default:
				assert("SPH_CD::PPC_FPM_General" && 0);
				break;
			}
			break;
		case(2):
			switch (r) {
			case(1):
			case(2):
				derivativeM = (pp->dWd(static_cast<DiffAxis>(r - 1), pp->Mpart_ptr));
				derivativeNB = (pp->dWd(static_cast<DiffAxis>(r - 1), pp->NBpart_ptr));
				break;
			case(3):
			case(4):
				derivativeM = (pp->d2Wd(static_cast<DiffAxis>(r - nrOfDim - 1), static_cast<DiffAxis>(r - nrOfDim - 1), pp->Mpart_ptr));
				derivativeNB = (pp->d2Wd(static_cast<DiffAxis>(r - nrOfDim - 1), static_cast<DiffAxis>(r - nrOfDim - 1), pp->NBpart_ptr));
				break;
			case(5):
				derivativeM = (pp->d2Wd(X, Y, pp->Mpart_ptr));
				derivativeNB = (pp->d2Wd(X, Y, pp->NBpart_ptr));
				break;
			default:
				assert("SPH_CD::PPC_FPM_General" && 0);
				break;
			}
			break;
		case(3):
			switch (r) {
			case(1):
			case(2):
			case(3):
				derivativeM = (pp->dWd(static_cast<DiffAxis>(r - 1), pp->Mpart_ptr));
				derivativeNB = (pp->dWd(static_cast<DiffAxis>(r - 1), pp->NBpart_ptr));
				break;
			case(4):
			case(5):
			case(6):
				derivativeM = (pp->d2Wd(static_cast<DiffAxis>(r - nrOfDim - 1), static_cast<DiffAxis>(r - nrOfDim - 1), pp->Mpart_ptr));
				derivativeNB = (pp->d2Wd(static_cast<DiffAxis>(r - nrOfDim - 1), static_cast<DiffAxis>(r - nrOfDim - 1), pp->NBpart_ptr));
				break;
			case(7):
				derivativeM = (pp->d2Wd(X, Y, pp->Mpart_ptr));
				derivativeNB = (pp->d2Wd(X, Y, pp->NBpart_ptr));
				break;
			case(8):
				derivativeM = (pp->d2Wd(X, Z, pp->Mpart_ptr));
				derivativeNB = (pp->d2Wd(X, Z, pp->NBpart_ptr));
				break;
			case(9):
				derivativeM = (pp->d2Wd(Y, Z, pp->Mpart_ptr));
				derivativeNB = (pp->d2Wd(Y, Z, pp->NBpart_ptr));
				break;
			default:
				assert("SPH_CD::PPC_FPM_General" && 0);
				break;
			}
			break;
		default:
			break;
		}
		for (int c = 0;c < m_size;c++) {
			switch (nrOfDim) {
			case(1):
				switch (c) {
				case(0):
					tempVecM[c] = derivativeM * (pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
					tempVecNB[c] = derivativeNB * (pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
					break;
				case(1):
					tempVecM[c] = Xji[c - 1] * derivativeM * (pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
					tempVecNB[c] = Xij[c - 1] * derivativeNB * (pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
					break;
				case(2):
					tempVecM[c] = pow(Xji[c - nrOfDim - 1], 2) / 2.f * derivativeM * (pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
					tempVecNB[c] = pow(Xij[c - nrOfDim - 1], 2) / 2.f * derivativeNB * (pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
					break;
				default:
					break;
				}
				break;
			case(2):
				switch (c) {
				case(0):
					tempVecM[c] = derivativeM * (pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
					tempVecNB[c] = derivativeNB * (pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
					break;
				case(1):
				case(2):
					tempVecM[c] = Xji[c - 1] * derivativeM * (pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
					tempVecNB[c] = Xij[c - 1] * derivativeNB * (pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
					break;
				case(3):
				case(4):
					tempVecM[c] = pow(Xji[c - nrOfDim - 1], 2) / 2.0 * derivativeM * (pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
					tempVecNB[c] = pow(Xij[c - nrOfDim - 1], 2) / 2.0 * derivativeNB * (pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
					break;
				case(5):
					tempVecM[c] = Xji[X] * Xji[Y] * derivativeM * (pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
					tempVecNB[c] = Xij[X] * Xij[Y] * derivativeNB * (pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
					break;
				default:
					break;
				}
				break;
			case(3):
				switch (c) {
				case(0):
					tempVecM[c] = derivativeM * (pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
					tempVecNB[c] = derivativeNB * (pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
					break;
				case(1):
				case(2):
				case(3):
					tempVecM[c] = Xji[c - 1] * derivativeM * (pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
					tempVecNB[c] = Xij[c - 1] * derivativeNB * (pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
					break;
				case(4):
				case(5):
				case(6):
					tempVecM[c] = pow(Xji[c - nrOfDim - 1], 2)/2.f * derivativeM * (pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
					tempVecNB[c] = pow(Xij[c - nrOfDim - 1], 2)/2.f * derivativeNB * (pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
					break;
				case(7):
					tempVecM[c] = Xji[X] * Xji[Y] * derivativeM * (pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
					tempVecNB[c] = Xij[X] * Xij[Y] * derivativeNB * (pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
					break;
				case(8):
					tempVecM[c] = Xji[X] * Xji[Z] * derivativeM * (pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
					tempVecNB[c] = Xij[X] * Xij[Z] * derivativeNB * (pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
					break;
				case(9):
					tempVecM[c] = Xji[Y] * Xji[Z] * derivativeM * (pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
					tempVecNB[c] = Xij[Y] * Xij[Z] * derivativeNB * (pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
					break;
				default:
					break;
				}
				break;
			default:
				break;
			}
		}

		//matM[r] = tempVecM;
		//matNB[r] = tempVecNB;
		for (int c = 0;c < m_size;c++) {
			matM(r, c) = tempVecM[c];
			matNB(r, c) = tempVecNB[c];
		}

	}


	/*
	for (int r = 0;r < m_size;r++) {
		for (int c = 0;c < m_size;c++) {
			(*FPM_matrix)[pp->Mpart_ptr->m_id](r,c) += matM(r,c);
			(*FPM_matrix)[pp->NBpart_ptr->m_id](r,c) += matNB(r,c);
		}
	}
	*/

	//std::cout << matM << "\n";
	//std::cout << matNB << "\n";

	(*FPM_matrix)[pp->Mpart_ptr->m_id] += matM;
	(*FPM_matrix)[pp->NBpart_ptr->m_id] += matNB;


}
void SPH_CD::PPC_FPMStandart_vec(ParticlePair* pp, std::vector<std::vector<part_prec>>* FPM_vec,std::string param) {


	float MpParam;
	float NBpParam;

	if (param == "Vx") { 
		MpParam = pp->Mpart_ptr->m_velocity.val.x;
		NBpParam = pp->NBpart_ptr->m_velocity.val.x;
	}else if (param == "Vy") {
		MpParam = pp->Mpart_ptr->m_velocity.val.y;
		NBpParam = pp->NBpart_ptr->m_velocity.val.y;
	}else if (param == "Vz") {
		MpParam = pp->Mpart_ptr->m_velocity.val.z;
		NBpParam = pp->NBpart_ptr->m_velocity.val.z;
	}else if (param == "Pressure") {
		MpParam = pp->Mpart_ptr->m_pressure.val;
		NBpParam = pp->NBpart_ptr->m_pressure.val;
	}else if (param == "h") {
		MpParam = pp->Mpart_ptr->m_SmR;
		NBpParam = pp->NBpart_ptr->m_SmR;
	}



	int m_size = (*FPM_vec)[0].size();

	std::vector<part_prec> vecM;         std::vector<part_prec> vecNB;
	vecM.resize(m_size,0.f);         vecNB.resize(m_size,0.f);

	vecM[0] = NBpParam *pp->W(pp->Mpart_ptr)*(pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
	vecNB[0] = MpParam *pp->W(pp->NBpart_ptr)*(pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);

	for (int d = 0;d < nrOfDim;d++) {
		vecM[d + 1] = NBpParam *(pp->dWd(static_cast<DiffAxis>(d), pp->Mpart_ptr))*(pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
		vecNB[d + 1] = MpParam *(pp->dWd(static_cast<DiffAxis>(d), pp->NBpart_ptr))*(pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);


		vecM[d + nrOfDim + 1] = NBpParam *(pp->d2Wd(static_cast<DiffAxis>(d), static_cast<DiffAxis>(d), pp->Mpart_ptr))*(pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
		vecNB[d + nrOfDim + 1] = MpParam *(pp->d2Wd(static_cast<DiffAxis>(d), static_cast<DiffAxis>(d), pp->NBpart_ptr))*(pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
	}

	if (nrOfDim == 2) {
		vecM[2*nrOfDim + 1] = NBpParam *(pp->d2Wd(X, Y, pp->Mpart_ptr))*(pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
		vecNB[2*nrOfDim + 1] = MpParam *(pp->d2Wd(X, Y, pp->NBpart_ptr))*(pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
	}
	else if(nrOfDim == 3)	{
		vecM[2 * nrOfDim + 1] = NBpParam *(pp->d2Wd(X, Y, pp->Mpart_ptr))*(pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
		vecNB[2 * nrOfDim + 1] = MpParam *(pp->d2Wd(X, Y, pp->NBpart_ptr))*(pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);

		vecM[2 * nrOfDim + 2] = NBpParam *(pp->d2Wd(X, Z, pp->Mpart_ptr))*(pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
		vecNB[2 * nrOfDim + 2] = MpParam *(pp->d2Wd(X, Z, pp->NBpart_ptr))*(pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);

		vecM[2 * nrOfDim + 3] = NBpParam *(pp->d2Wd(Y, Z, pp->Mpart_ptr))*(pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
		vecNB[2 * nrOfDim + 3] = MpParam *(pp->d2Wd(Y, Z, pp->NBpart_ptr))*(pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
	}

	//std::cout << &((*FPM_vec)[pp->Mpart_ptr->m_id]) << "   " << &(vecM) << "\n";


	for (int i = 0;i < vecM.size();i++) {
		//std::cout << ((*FPM_vec)[pp->Mpart_ptr->m_id][i]);
		//std::cout << "   " << (vecM[i]) << "\n";
		(*FPM_vec)[pp->Mpart_ptr->m_id][i] += vecM[i];
		(*FPM_vec)[pp->NBpart_ptr->m_id][i] += vecNB[i];
	}
}
void SPH_CD::PPC_FPMDifference_General(ParticlePair* pp, std::vector<Matrix>* FPM_matrix) {


	int m_size = (*FPM_matrix)[0].getRows();

	Matrix matM;                    Matrix matNB;
	matM.resize(m_size, m_size);    matNB.resize(m_size, m_size);

	std::vector<float> tempVecM;     std::vector<float> tempVecNB;
	tempVecM.resize(m_size, 0.f); tempVecNB.resize(m_size, 0.f);

	part_prec_3 Xji = { pp->NBpart_ptr->dx(pp->Mpart_ptr),pp->NBpart_ptr->dy(pp->Mpart_ptr),pp->NBpart_ptr->dz(pp->Mpart_ptr) };
	part_prec_3 Xij = { pp->Mpart_ptr->dx(pp->NBpart_ptr),pp->Mpart_ptr->dy(pp->NBpart_ptr),pp->Mpart_ptr->dz(pp->NBpart_ptr) };

	tempVecM[0] = pp->W(pp->Mpart_ptr)*(pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
	tempVecNB[0] = pp->W(pp->NBpart_ptr)*(pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);


	for (int dd = 0;dd < nrOfDim;dd++) {
		tempVecM[dd + 1] = Xji[dd] * pp->W(pp->Mpart_ptr)*(pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
		tempVecNB[dd + 1] = Xij[dd] * pp->W(pp->NBpart_ptr)*(pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);

		tempVecM[dd + nrOfDim + 1] = pow(Xji[dd], 2) / 2.f * pp->W(pp->Mpart_ptr)*(pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
		tempVecNB[dd + nrOfDim + 1] = pow(Xij[dd], 2) / 2.f * pp->W(pp->NBpart_ptr)*(pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
	}

	if (nrOfDim == 2) {
		tempVecM[2 * nrOfDim + 1] = Xji[X] * Xji[Y] * pp->W(pp->Mpart_ptr)*(pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
		tempVecNB[2 * nrOfDim + 1] = Xij[X] * Xij[Y] * pp->W(pp->NBpart_ptr)*(pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
	}
	else if (nrOfDim == 3) {
		tempVecM[2 * nrOfDim + 1] = Xji[X] * Xji[Y] * pp->W(pp->Mpart_ptr)*(pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
		tempVecNB[2 * nrOfDim + 1] = Xij[X] * Xij[Y] * pp->W(pp->NBpart_ptr)*(pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);

		tempVecM[2 * nrOfDim + 2] = Xji[X] * Xji[Z] * pp->W(pp->Mpart_ptr)*(pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
		tempVecNB[2 * nrOfDim + 2] = Xij[X] * Xij[Z] * pp->W(pp->NBpart_ptr)*(pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);

		tempVecM[2 * nrOfDim + 3] = Xji[Y] * Xji[Z] * pp->W(pp->Mpart_ptr)*(pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
		tempVecNB[2 * nrOfDim + 3] = Xij[Y] * Xij[Z] * pp->W(pp->NBpart_ptr)*(pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
	}


	for (int c = 0;c < m_size;c++) {


		matM(0, c) = tempVecM[c];
		matNB(0, c) = tempVecNB[c];
	}



	float derivativeM;
	float derivativeNB;

	for (int r = 1;r < m_size;r++) {
		//Определяем порядок производной
		switch (nrOfDim) {
		case(1):
			switch (r) {
			case(1):
				derivativeM = (pp->dWd(static_cast<DiffAxis>(r - 1), pp->Mpart_ptr));
				derivativeNB = (pp->dWd(static_cast<DiffAxis>(r - 1), pp->NBpart_ptr));
				break;
			case(2):
				derivativeM = (pp->d2Wd(static_cast<DiffAxis>(r - nrOfDim - 1), static_cast<DiffAxis>(r - nrOfDim - 1), pp->Mpart_ptr));
				derivativeNB = (pp->d2Wd(static_cast<DiffAxis>(r - nrOfDim - 1), static_cast<DiffAxis>(r - nrOfDim - 1), pp->NBpart_ptr));
				break;
			default:
				assert("SPH_CD::PPC_FPM_General" && 0);
				break;
			}
			break;
		case(2):
			switch (r) {
			case(1):
			case(2):
				derivativeM = (pp->dWd(static_cast<DiffAxis>(r - 1), pp->Mpart_ptr));
				derivativeNB = (pp->dWd(static_cast<DiffAxis>(r - 1), pp->NBpart_ptr));
				break;
			case(3):
			case(4):
				derivativeM = (pp->d2Wd(static_cast<DiffAxis>(r - nrOfDim - 1), static_cast<DiffAxis>(r - nrOfDim - 1), pp->Mpart_ptr));
				derivativeNB = (pp->d2Wd(static_cast<DiffAxis>(r - nrOfDim - 1), static_cast<DiffAxis>(r - nrOfDim - 1), pp->NBpart_ptr));
				break;
			case(5):
				derivativeM = (pp->d2Wd(X, Y, pp->Mpart_ptr));
				derivativeNB = (pp->d2Wd(X, Y, pp->NBpart_ptr));
				break;
			default:
				assert("SPH_CD::PPC_FPM_General" && 0);
				break;
			}
			break;
		case(3):
			switch (r) {
			case(1):
			case(2):
			case(3):
				derivativeM = (pp->dWd(static_cast<DiffAxis>(r - 1), pp->Mpart_ptr));
				derivativeNB = (pp->dWd(static_cast<DiffAxis>(r - 1), pp->NBpart_ptr));
				break;
			case(4):
			case(5):
			case(6):
				derivativeM = (pp->d2Wd(static_cast<DiffAxis>(r - nrOfDim - 1), static_cast<DiffAxis>(r - nrOfDim - 1), pp->Mpart_ptr));
				derivativeNB = (pp->d2Wd(static_cast<DiffAxis>(r - nrOfDim - 1), static_cast<DiffAxis>(r - nrOfDim - 1), pp->NBpart_ptr));
				break;
			case(7):
				derivativeM = (pp->d2Wd(X, Y, pp->Mpart_ptr));
				derivativeNB = (pp->d2Wd(X, Y, pp->NBpart_ptr));
				break;
			case(8):
				derivativeM = (pp->d2Wd(X, Z, pp->Mpart_ptr));
				derivativeNB = (pp->d2Wd(X, Z, pp->NBpart_ptr));
				break;
			case(9):
				derivativeM = (pp->d2Wd(Y, Z, pp->Mpart_ptr));
				derivativeNB = (pp->d2Wd(Y, Z, pp->NBpart_ptr));
				break;
			default:
				assert("SPH_CD::PPC_FPM_General" && 0);
				break;
			}
			break;
		default:
			break;
		}
		for (int c = 0;c < m_size;c++) {
			switch (nrOfDim) {
			case(1):
				switch (c) {
				case(0):
					tempVecM[c] = derivativeM * (pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
					tempVecNB[c] = derivativeNB * (pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
					break;
				case(1):
					tempVecM[c] = Xji[c - 1] * derivativeM * (pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
					tempVecNB[c] = Xij[c - 1] * derivativeNB * (pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
					break;
				case(2):
					tempVecM[c] = pow(Xji[c - nrOfDim - 1], 2) / 2.f * derivativeM * (pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
					tempVecNB[c] = pow(Xij[c - nrOfDim - 1], 2) / 2.f * derivativeNB * (pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
					break;
				default:
					break;
				}
				break;
			case(2):
				switch (c) {
				case(0):
					tempVecM[c] = derivativeM * (pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
					tempVecNB[c] = derivativeNB * (pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
					break;
				case(1):
				case(2):
					tempVecM[c] = Xji[c - 1] * derivativeM * (pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
					tempVecNB[c] = Xij[c - 1] * derivativeNB * (pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
					break;
				case(3):
				case(4):
					tempVecM[c] = pow(Xji[c - nrOfDim - 1], 2) / 2.f * derivativeM * (pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
					tempVecNB[c] = pow(Xij[c - nrOfDim - 1], 2) / 2.f * derivativeNB * (pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
					break;
				case(5):
					tempVecM[c] = Xji[X] * Xji[Y] * derivativeM * (pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
					tempVecNB[c] = Xij[X] * Xij[Y] * derivativeNB * (pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
					break;
				default:
					break;
				}
				break;
			case(3):
				switch (c) {
				case(0):
					tempVecM[c] = derivativeM * (pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
					tempVecNB[c] = derivativeNB * (pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
					break;
				case(1):
				case(2):
				case(3):
					tempVecM[c] = Xji[c - 1] * derivativeM * (pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
					tempVecNB[c] = Xij[c - 1] * derivativeNB * (pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
					break;
				case(4):
				case(5):
				case(6):
					tempVecM[c] = pow(Xji[c - nrOfDim - 1], 2) / 2.f * derivativeM * (pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
					tempVecNB[c] = pow(Xij[c - nrOfDim - 1], 2) / 2.f * derivativeNB * (pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
					break;
				case(7):
					tempVecM[c] = Xji[X] * Xji[Y] * derivativeM * (pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
					tempVecNB[c] = Xij[X] * Xij[Y] * derivativeNB * (pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
					break;
				case(8):
					tempVecM[c] = Xji[X] * Xji[Z] * derivativeM * (pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
					tempVecNB[c] = Xij[X] * Xij[Z] * derivativeNB * (pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
					break;
				case(9):
					tempVecM[c] = Xji[Y] * Xji[Z] * derivativeM * (pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
					tempVecNB[c] = Xij[Y] * Xij[Z] * derivativeNB * (pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
					break;
				default:
					break;
				}
				break;
			default:
				break;
			}
		}

		//matM[r] = tempVecM;
		//matNB[r] = tempVecNB;
		for (int c = 0;c < m_size;c++) {
			matM(r, c) = tempVecM[c];
			matNB(r, c) = tempVecNB[c];
		}

	}


	/*
	for (int r = 0;r < m_size;r++) {
		for (int c = 0;c < m_size;c++) {
			(*FPM_matrix)[pp->Mpart_ptr->m_id](r,c) += matM(r,c);
			(*FPM_matrix)[pp->NBpart_ptr->m_id](r,c) += matNB(r,c);
		}
	}
	*/

	//std::cout << matM << "\n";
	//std::cout << matNB << "\n";

	(*FPM_matrix)[pp->Mpart_ptr->m_id] += matM;
	(*FPM_matrix)[pp->NBpart_ptr->m_id] += matNB;



	for (int r = 1;r < m_size;r++) {
		(*FPM_matrix)[pp->Mpart_ptr->m_id](r, 0) = 1.f;
	}


}
void SPH_CD::PPC_FPMDifference_vec(ParticlePair* pp, std::vector<std::vector<part_prec>>* FPM_vec, std::string param) {


	float MpParam;
	float NBpParam;

	if (param == "Vx") {
		MpParam = pp->Mpart_ptr->m_velocity.val.x;
		NBpParam = pp->NBpart_ptr->m_velocity.val.x;
	}
	if (param == "Vy") {
		MpParam = pp->Mpart_ptr->m_velocity.val.y;
		NBpParam = pp->NBpart_ptr->m_velocity.val.y;
	}
	if (param == "Vz") {
		MpParam = pp->Mpart_ptr->m_velocity.val.z;
		NBpParam = pp->NBpart_ptr->m_velocity.val.z;
	}
	if (param == "Pressure") {
		MpParam = pp->Mpart_ptr->m_pressure.val;
		NBpParam = pp->NBpart_ptr->m_pressure.val;
	}



	int m_size = (*FPM_vec)[0].size();

	std::vector<float> vecM;         std::vector<float> vecNB;
	vecM.resize(m_size, 0.f);         vecNB.resize(m_size, 0.f);

	vecM[0] = (NBpParam - MpParam) * pp->W(pp->Mpart_ptr)*(pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
	vecNB[0] = (MpParam - NBpParam) * pp->W(pp->NBpart_ptr)*(pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);

	for (int d = 0;d < nrOfDim;d++) {
		vecM[d + 1] = (NBpParam - MpParam) * (pp->dWd(static_cast<DiffAxis>(d), pp->Mpart_ptr))*(pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
		vecNB[d + 1] = (MpParam - NBpParam) * (pp->dWd(static_cast<DiffAxis>(d), pp->NBpart_ptr))*(pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);


		vecM[d + nrOfDim + 1] = (NBpParam - MpParam) * (pp->d2Wd(static_cast<DiffAxis>(d), static_cast<DiffAxis>(d), pp->Mpart_ptr))*(pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
		vecNB[d + nrOfDim + 1] = (MpParam - NBpParam) * (pp->d2Wd(static_cast<DiffAxis>(d), static_cast<DiffAxis>(d), pp->NBpart_ptr))*(pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
	}

	if (nrOfDim == 2) {
		vecM[2 * nrOfDim + 1] = (NBpParam - MpParam) * (pp->d2Wd(X, Y, pp->Mpart_ptr))*(pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
		vecNB[2 * nrOfDim + 1] = (MpParam - NBpParam) * (pp->d2Wd(X, Y, pp->NBpart_ptr))*(pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
	}
	else if (nrOfDim == 3) {
		vecM[2 * nrOfDim + 1] = (NBpParam - MpParam) * (pp->d2Wd(X, Y, pp->Mpart_ptr))*(pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
		vecNB[2 * nrOfDim + 1] = (MpParam - NBpParam) * (pp->d2Wd(X, Y, pp->NBpart_ptr))*(pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);

		vecM[2 * nrOfDim + 2] = (NBpParam - MpParam) * (pp->d2Wd(X, Z, pp->Mpart_ptr))*(pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
		vecNB[2 * nrOfDim + 2] = (MpParam - NBpParam) * (pp->d2Wd(X, Z, pp->NBpart_ptr))*(pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);

		vecM[2 * nrOfDim + 3] = (NBpParam - MpParam) * (pp->d2Wd(Y, Z, pp->Mpart_ptr))*(pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
		vecNB[2 * nrOfDim + 3] = (MpParam - NBpParam) * (pp->d2Wd(Y, Z, pp->NBpart_ptr))*(pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
	}

	//std::cout << &((*FPM_vec)[pp->Mpart_ptr->m_id]) << "   " << &(vecM) << "\n";


	for (int i = 0;i < vecM.size();i++) {
		//std::cout << ((*FPM_vec)[pp->Mpart_ptr->m_id][i]);
		//std::cout << "   " << (vecM[i]) << "\n";
		(*FPM_vec)[pp->Mpart_ptr->m_id][i] += vecM[i];
		(*FPM_vec)[pp->NBpart_ptr->m_id][i] += vecNB[i];
	}




}
void SPH_CD::PPC_CorSPH_General(ParticlePair* pp, std::vector<Matrix>* FPM_matrix) {
	int m_size = (*FPM_matrix)[0].getRows();
	Matrix matM;                    Matrix matNB;
	matM.resize(m_size, m_size);    matNB.resize(m_size, m_size);

	part_prec_3 Xji = { pp->NBpart_ptr->dx(pp->Mpart_ptr),pp->NBpart_ptr->dy(pp->Mpart_ptr),pp->NBpart_ptr->dz(pp->Mpart_ptr) };
	part_prec_3 Xij = { pp->Mpart_ptr->dx(pp->NBpart_ptr),pp->Mpart_ptr->dy(pp->NBpart_ptr),pp->Mpart_ptr->dz(pp->NBpart_ptr) };

	for (int d = 0;d < nrOfDim;d++) {
		for (int dd = 0;dd < nrOfDim;dd++) {
			matM(d, dd) = Xji[dd] * (pp->dWd(static_cast<DiffAxis>(d), pp->Mpart_ptr))*(pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
			matNB(d, dd) = Xij[dd] * (pp->dWd(static_cast<DiffAxis>(d), pp->NBpart_ptr))*(pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
	
		}

	}
	//if (pp->Mpart_ptr->m_id == 0) std::cout << matM;
	//if (pp->NBpart_ptr->m_id == 0) std::cout << matNB;
	(*FPM_matrix)[pp->Mpart_ptr->m_id] += matM;
	(*FPM_matrix)[pp->NBpart_ptr->m_id] += matNB;
}
void SPH_CD::PPC_CorSPH_vec(ParticlePair* pp, std::vector<std::vector<part_prec>>* FPM_vec, std::string param) {
	part_prec MpParam;
	part_prec NBpParam;
	if (param == "h") {
		MpParam = pp->Mpart_ptr->m_SmR;
		NBpParam = pp->NBpart_ptr->m_SmR;
	}
	int m_size = (*FPM_vec)[0].size();
	std::vector<part_prec> vecM;         std::vector<part_prec> vecNB;
	vecM.resize(m_size, 0.f);         vecNB.resize(m_size, 0.f);
	for (int d = 0;d < nrOfDim;d++) {
		vecM[d] = (NBpParam - MpParam) * (pp->dWd(static_cast<DiffAxis>(d), pp->Mpart_ptr))*(pp->NBpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val);
		vecNB[d] = (MpParam - NBpParam) * (pp->dWd(static_cast<DiffAxis>(d), pp->NBpart_ptr))*(pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
	}
	for (int i = 0;i < vecM.size();i++) {
		//std::cout << ((*FPM_vec)[pp->Mpart_ptr->m_id][i]);
		//std::cout << "   " << (vecM[i]) << "\n";
		(*FPM_vec)[pp->Mpart_ptr->m_id][i] += vecM[i];
		(*FPM_vec)[pp->NBpart_ptr->m_id][i] += vecNB[i];
	}
}


void SPH_CD::RenormFactorCalc(std::vector<part_prec>* gamma, ParticlePair * pp) {
	(*gamma)[pp->Mpart_ptr->m_id] += pp->NBpart_ptr->m_mass / pp->NBpart_ptr->m_density.val * pp->W(pp->Mpart_ptr);
	(*gamma)[pp->NBpart_ptr->m_id] += pp->Mpart_ptr->m_mass / pp->Mpart_ptr->m_density.val *pp->W(pp->NBpart_ptr);

	//if (pp->NBpart_ptr->m_mass / pp->NBpart_ptr->m_density.val * pp->W(pp->Mpart_ptr) != 0) {
	//	std::cout << pp->NBpart_ptr->m_mass / pp->NBpart_ptr->m_density.val * pp->W(pp->Mpart_ptr)<< " \n";
	//}
	//std::cout << drhodt[i->m_id] << "   " << gamma[i->m_id] << "   " << i->m_density.val << "   " << dot(gamma_deriv[i->m_id], i->m_velocity.val) << "\n";
}
void SPH_CD::RenormFactorDerivCalc(std::vector<part_prec_3>* gamma_deriv, ParticlePair * pp) {
	(*gamma_deriv)[pp->Mpart_ptr->m_id] += pp->NBpart_ptr->m_mass / pp->NBpart_ptr->m_density.val * pp->vec_e(pp->Mpart_ptr) * pp->dWd(R, pp->Mpart_ptr);
	(*gamma_deriv)[pp->NBpart_ptr->m_id] += pp->Mpart_ptr->m_mass / pp->Mpart_ptr->m_density.val *pp->vec_e(pp->NBpart_ptr) *pp->dWd(R, pp->NBpart_ptr);
}


void SPH_CD::DensityCalculation(std::vector<part_prec>* drhodt, ParticlePair * pp) {

	switch (m_options.densityChangeAlgorithm) {
	case(0):
		// k=0
		switch (m_options.densityChangeUsing) {
		case(VAL):
			(*drhodt)[pp->Mpart_ptr->m_id] += pp->NBpart_ptr->m_mass * pp->Mpart_ptr->m_density.val / pp->NBpart_ptr->m_density.val * pp->W(pp->Mpart_ptr);
			(*drhodt)[pp->NBpart_ptr->m_id] += pp->Mpart_ptr->m_mass * pp->NBpart_ptr->m_density.val / pp->Mpart_ptr->m_density.val *pp->W(pp->NBpart_ptr);
			break;
		case(DVAL_DT):
			(*drhodt)[pp->Mpart_ptr->m_id] += pp->Mpart_ptr->m_density.val*pp->NBpart_ptr->m_mass / pp->NBpart_ptr->m_density.val*glm::dot(pp->Mpart_ptr->dV(pp->NBpart_ptr), pp->dWd(R, pp->Mpart_ptr)*pp->vec_e(pp->Mpart_ptr));
			(*drhodt)[pp->NBpart_ptr->m_id] += pp->NBpart_ptr->m_density.val*pp->Mpart_ptr->m_mass / pp->Mpart_ptr->m_density.val*glm::dot(pp->NBpart_ptr->dV(pp->Mpart_ptr), pp->dWd(R, pp->NBpart_ptr)*pp->vec_e(pp->NBpart_ptr));
			break;
		default:
			break;
		}
		break;
	case(1):
		// k=1
		switch (m_options.densityChangeUsing) {
		case(VAL):
			(*drhodt)[pp->Mpart_ptr->m_id] += pp->NBpart_ptr->m_mass * pp->W(pp->Mpart_ptr);
			(*drhodt)[pp->NBpart_ptr->m_id] += pp->Mpart_ptr->m_mass * pp->W(pp->NBpart_ptr);
			break;
		case(DVAL_DT):
			(*drhodt)[pp->Mpart_ptr->m_id] += pp->NBpart_ptr->m_mass * glm::dot(pp->Mpart_ptr->dV(pp->NBpart_ptr), pp->dWd(R, pp->Mpart_ptr)* pp->vec_e(pp->Mpart_ptr));
			(*drhodt)[pp->NBpart_ptr->m_id] += pp->Mpart_ptr->m_mass * glm::dot(pp->NBpart_ptr->dV(pp->Mpart_ptr), pp->dWd(R, pp->NBpart_ptr)*pp->vec_e(pp->NBpart_ptr));
			break;
		default:
			break;
		}
		break;
	case(2):
		// k=2
		switch (m_options.densityChangeUsing) {
		case(VAL):
			(*drhodt)[pp->Mpart_ptr->m_id] += pp->NBpart_ptr->m_mass * pp->NBpart_ptr->m_density.val / pp->Mpart_ptr->m_density.val * pp->W(pp->Mpart_ptr);
			(*drhodt)[pp->NBpart_ptr->m_id] += pp->Mpart_ptr->m_mass * pp->Mpart_ptr->m_density.val / pp->NBpart_ptr->m_density.val *pp->W(pp->NBpart_ptr);
			break;
		case(DVAL_DT):
			(*drhodt)[pp->Mpart_ptr->m_id] += pp->Mpart_ptr->m_density.val*pp->NBpart_ptr->m_mass / pp->NBpart_ptr->m_density.val*glm::dot(pp->Mpart_ptr->dV(pp->NBpart_ptr), pp->dWd(R, pp->Mpart_ptr)*pp->vec_e(pp->Mpart_ptr));
			(*drhodt)[pp->NBpart_ptr->m_id] += pp->NBpart_ptr->m_density.val*pp->Mpart_ptr->m_mass / pp->Mpart_ptr->m_density.val*glm::dot(pp->NBpart_ptr->dV(pp->Mpart_ptr), pp->dWd(R, pp->NBpart_ptr)*pp->vec_e(pp->NBpart_ptr));
			break;
		default:
			break;
		}
		break;
	default:
		break;
	}

	//if ((pp->NBpart_ptr->m_id == m_options.nrOfParticles[PARTICLETYPE::REAL] - 1)) {
	//	std::cout << "Density (" << pp->Mpart_ptr->m_id << ") = (" << pp->NBpart_ptr->m_mass * glm::dot(pp->Mpart_ptr->dV(pp->NBpart_ptr), pp->dWd(R, pp->Mpart_ptr)* pp->vec_e(pp->Mpart_ptr)) << ")    -    (" << pp->NBpart_ptr->m_id << ") = (" << pp->Mpart_ptr->m_mass * glm::dot(pp->NBpart_ptr->dV(pp->Mpart_ptr), pp->dWd(R, pp->NBpart_ptr)*pp->vec_e(pp->NBpart_ptr)) << ")\n";
	//}
	//if ((pp->Mpart_ptr->m_id == m_options.nrOfParticles[PARTICLETYPE::REAL] - 1)) {
	//	std::cout << "Density (" << pp->NBpart_ptr->m_id << ") = (" << pp->Mpart_ptr->m_mass * glm::dot(pp->NBpart_ptr->dV(pp->Mpart_ptr), pp->dWd(R, pp->NBpart_ptr)*pp->vec_e(pp->NBpart_ptr)) << ")    -    (" << pp->Mpart_ptr->m_id << ") = (" << pp->NBpart_ptr->m_mass * glm::dot(pp->Mpart_ptr->dV(pp->NBpart_ptr), pp->dWd(R, pp->Mpart_ptr)* pp->vec_e(pp->Mpart_ptr)) << ")\n";
	//}




	//if ((pp->Mpart_ptr->m_id == m_options.nrOfParticles[PARTICLETYPE::REAL] - 1) and (pp->NBpart_ptr->m_mass * glm::dot(pp->Mpart_ptr->dV(pp->NBpart_ptr), pp->dWd(R, pp->Mpart_ptr)*pp->vec_e(pp->Mpart_ptr)) != 0)) {
	//	//std::cout << "(" << pp->NBpart_ptr->m_id << ")(" << pp->NBpart_ptr->m_mass << "*" << glm::dot(pp->Mpart_ptr->dV(pp->NBpart_ptr), pp->dWd(R, pp->Mpart_ptr)*pp->vec_e(pp->Mpart_ptr)) << ")+ ";
	//	std::cout << "(" << pp->NBpart_ptr->m_id << ")(" << pp->NBpart_ptr->m_mass * glm::dot(pp->Mpart_ptr->dV(pp->NBpart_ptr), pp->dWd(R, pp->Mpart_ptr)*pp->vec_e(pp->Mpart_ptr)) << ")\n";
	//}
	//if ((pp->NBpart_ptr->m_id == m_options.nrOfParticles[PARTICLETYPE::REAL] - 1) and (pp->Mpart_ptr->m_mass * glm::dot(pp->NBpart_ptr->dV(pp->Mpart_ptr), pp->dWd(R, pp->NBpart_ptr)* pp->vec_e(pp->NBpart_ptr)) != 0)) {
	//	//std::cout << "(" << pp->Mpart_ptr->m_id << ")(" << pp->Mpart_ptr->m_mass << "*" << glm::dot(pp->NBpart_ptr->dV(pp->Mpart_ptr), pp->dWd(R, pp->NBpart_ptr)*pp->vec_e(pp->NBpart_ptr)) << ")+ ";
	//	std::cout << "(" << pp->Mpart_ptr->m_id << ")(" << pp->Mpart_ptr->m_mass * glm::dot(pp->NBpart_ptr->dV(pp->Mpart_ptr), pp->dWd(R, pp->NBpart_ptr)*pp->vec_e(pp->NBpart_ptr)) << ")\n";
	//}
	//if (pp->NBpart_ptr->m_mass * glm::dot(pp->Mpart_ptr->dV(pp->NBpart_ptr), pp->dWd(R, pp->Mpart_ptr)*pp->vec_e(pp->Mpart_ptr)) != 0) {
	//	std::cout << pp->NBpart_ptr->m_mass * glm::dot(pp->Mpart_ptr->dV(pp->NBpart_ptr), pp->dWd(R, pp->Mpart_ptr)*pp->vec_e(pp->Mpart_ptr)) << " = "<< pp->Mpart_ptr->m_mass * glm::dot(pp->NBpart_ptr->dV(pp->Mpart_ptr), pp->dWd(R, pp->NBpart_ptr)*pp->vec_e(pp->NBpart_ptr)) <<"\n";
	//}

}

void SPH_CD::DensityDiffusiveTerm(std::vector<part_prec>* drhodt, ParticlePair * pp) {
	(*drhodt)[pp->Mpart_ptr->m_id] += dot((static_cast<part_prec>(2.0)*(pp->Mpart_ptr->m_density.val - pp->NBpart_ptr->m_density.val) / pp->Mpart_ptr->distance(pp->NBpart_ptr)*pp->vec_e(pp->Mpart_ptr)),(pp->dWd(R, pp->Mpart_ptr) * pp->vec_e(pp->Mpart_ptr)))*(pp->NBpart_ptr->m_mass)/(pp->NBpart_ptr->m_density.val);
	(*drhodt)[pp->NBpart_ptr->m_id] += dot((static_cast<part_prec>(2.0)*(pp->NBpart_ptr->m_density.val - pp->Mpart_ptr->m_density.val) / pp->NBpart_ptr->distance(pp->Mpart_ptr)*pp->vec_e(pp->NBpart_ptr)), (pp->dWd(R, pp->NBpart_ptr) * pp->vec_e(pp->NBpart_ptr)))*(pp->Mpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val);
}


void SPH_CD::PressurePartCalculation(std::vector<part_prec_3>* dvdt, ParticlePair * pp) {
	switch (m_options.velocityChangeAlgorithm_pressurePart){
	case(0):
		// k=0
		(*dvdt)[pp->Mpart_ptr->m_id] += (-pp->NBpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val*pp->NBpart_ptr->m_density.val)*(pp->Mpart_ptr->m_pressure.val + pp->NBpart_ptr->m_pressure.val)*pp->dWd(R, pp->Mpart_ptr) * pp->vec_e(pp->Mpart_ptr);
		(*dvdt)[pp->Mpart_ptr->m_id] += (-pp->Mpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val*pp->Mpart_ptr->m_density.val)*(pp->NBpart_ptr->m_pressure.val + pp->Mpart_ptr->m_pressure.val)*pp->dWd(R, pp->NBpart_ptr) * pp->vec_e(pp->NBpart_ptr);
		break;
	case(1):
		// k=1
		(*dvdt)[pp->Mpart_ptr->m_id] += (-pp->NBpart_ptr->m_mass)*(pp->Mpart_ptr->m_pressure.val / pow(pp->Mpart_ptr->m_density.val, 2) + pp->NBpart_ptr->m_pressure.val / pow(pp->NBpart_ptr->m_density.val, 2))*pp->dWd(R, pp->Mpart_ptr) * pp->vec_e(pp->Mpart_ptr);
		(*dvdt)[pp->NBpart_ptr->m_id] += (-pp->Mpart_ptr->m_mass)*(pp->NBpart_ptr->m_pressure.val / pow(pp->NBpart_ptr->m_density.val, 2) + pp->Mpart_ptr->m_pressure.val / pow(pp->Mpart_ptr->m_density.val, 2))*pp->dWd(R, pp->NBpart_ptr) * pp->vec_e(pp->NBpart_ptr);
		break;
	case(2):
		// k=2
		(*dvdt)[pp->Mpart_ptr->m_id] += (-pp->NBpart_ptr->m_mass)*(pp->Mpart_ptr->m_pressure.val*pp->NBpart_ptr->m_density.val / pow(pp->Mpart_ptr->m_density.val, 3) + pp->NBpart_ptr->m_pressure.val*pp->Mpart_ptr->m_density.val / pow(pp->NBpart_ptr->m_density.val, 3))*pp->dWd(R, pp->Mpart_ptr) * pp->vec_e(pp->Mpart_ptr);
		(*dvdt)[pp->NBpart_ptr->m_id] += (-pp->Mpart_ptr->m_mass)*(pp->NBpart_ptr->m_pressure.val*pp->Mpart_ptr->m_density.val / pow(pp->NBpart_ptr->m_density.val, 3) + pp->Mpart_ptr->m_pressure.val*pp->NBpart_ptr->m_density.val / pow(pp->Mpart_ptr->m_density.val, 3))*pp->dWd(R, pp->NBpart_ptr) * pp->vec_e(pp->NBpart_ptr);
		break;
	default:
		break;
	}

	//if ((pp->NBpart_ptr->m_id == m_options.nrOfParticles[PARTICLETYPE::REAL] - 31)) {
	//	std::cout << "Pressure part (" << pp->Mpart_ptr->m_id << ") = (" << (-pp->NBpart_ptr->m_mass) << "*(" <<pp->Mpart_ptr->m_pressure.val / pow(pp->Mpart_ptr->m_density.val, 2) << "+"<< pp->NBpart_ptr->m_pressure.val / pow(pp->NBpart_ptr->m_density.val, 2) << ")*" << pp->dWd(R, pp->Mpart_ptr) * pp->vec_e(pp->Mpart_ptr).x << ", " << (-pp->NBpart_ptr->m_mass)<<"*("<<pp->Mpart_ptr->m_pressure.val / pow(pp->Mpart_ptr->m_density.val, 2) <<"+"<< pp->NBpart_ptr->m_pressure.val / pow(pp->NBpart_ptr->m_density.val, 2)<<")*"<<pp->dWd(R, pp->Mpart_ptr) * pp->vec_e(pp->Mpart_ptr).y << ")    -    (" << pp->NBpart_ptr->m_id << ") = (" << (-pp->Mpart_ptr->m_mass)<<"*("<<pp->NBpart_ptr->m_pressure.val / pow(pp->NBpart_ptr->m_density.val, 2) <<"+"<< pp->Mpart_ptr->m_pressure.val / pow(pp->Mpart_ptr->m_density.val, 2)<<")*"<<pp->dWd(R, pp->NBpart_ptr) * pp->vec_e(pp->NBpart_ptr).x << ", " << (-pp->Mpart_ptr->m_mass) <<"*("<<pp->NBpart_ptr->m_pressure.val / pow(pp->NBpart_ptr->m_density.val, 2) <<"+"<< pp->Mpart_ptr->m_pressure.val / pow(pp->Mpart_ptr->m_density.val, 2)<<")*"<< pp->dWd(R, pp->NBpart_ptr) * pp->vec_e(pp->NBpart_ptr).y << ")\n";
	//}
	//if ((pp->Mpart_ptr->m_id == m_options.nrOfParticles[PARTICLETYPE::REAL] - 31)) {
	//	std::cout << "Pressure part (" << pp->NBpart_ptr->m_id << ") = (" << ((-pp->Mpart_ptr->m_mass)*(pp->NBpart_ptr->m_pressure.val / pow(pp->NBpart_ptr->m_density.val, 2) + pp->Mpart_ptr->m_pressure.val / pow(pp->Mpart_ptr->m_density.val, 2))*pp->dWd(R, pp->NBpart_ptr) * pp->vec_e(pp->NBpart_ptr)).x << ", " << ((-pp->Mpart_ptr->m_mass)*(pp->NBpart_ptr->m_pressure.val / pow(pp->NBpart_ptr->m_density.val, 2) + pp->Mpart_ptr->m_pressure.val / pow(pp->Mpart_ptr->m_density.val, 2))*pp->dWd(R, pp->NBpart_ptr) * pp->vec_e(pp->NBpart_ptr)).y << ")    -    (" << pp->Mpart_ptr->m_id << ") = (" << ((-pp->NBpart_ptr->m_mass)*(pp->Mpart_ptr->m_pressure.val / pow(pp->Mpart_ptr->m_density.val, 2) + pp->NBpart_ptr->m_pressure.val / pow(pp->NBpart_ptr->m_density.val, 2))*pp->dWd(R, pp->Mpart_ptr) * pp->vec_e(pp->Mpart_ptr)).x << ", " << ((-pp->NBpart_ptr->m_mass)*(pp->Mpart_ptr->m_pressure.val / pow(pp->Mpart_ptr->m_density.val, 2) + pp->NBpart_ptr->m_pressure.val / pow(pp->NBpart_ptr->m_density.val, 2))*pp->dWd(R, pp->Mpart_ptr) * pp->vec_e(pp->Mpart_ptr)).y << ")\n";
	//}


	//if (((-pp->NBpart_ptr->m_mass)*(pp->Mpart_ptr->m_pressure.val*pp->NBpart_ptr->m_density.val / pow(pp->Mpart_ptr->m_density.val, 3) + pp->NBpart_ptr->m_pressure.val*pp->Mpart_ptr->m_density.val / pow(pp->NBpart_ptr->m_density.val, 3))*pp->dWd(R, pp->Mpart_ptr) * pp->vec_e(pp->Mpart_ptr)).x != 0)
	//	std::cout << ((-pp->NBpart_ptr->m_mass)*(pp->Mpart_ptr->m_pressure.val*pp->NBpart_ptr->m_density.val / pow(pp->Mpart_ptr->m_density.val, 3) + pp->NBpart_ptr->m_pressure.val*pp->Mpart_ptr->m_density.val / pow(pp->NBpart_ptr->m_density.val, 3))*pp->dWd(R, pp->Mpart_ptr) * pp->vec_e(pp->Mpart_ptr)).x << "   " << ((-pp->Mpart_ptr->m_mass)*(pp->NBpart_ptr->m_pressure.val*pp->Mpart_ptr->m_density.val / pow(pp->NBpart_ptr->m_density.val, 3) + pp->Mpart_ptr->m_pressure.val*pp->NBpart_ptr->m_density.val / pow(pp->Mpart_ptr->m_density.val, 3))*pp->dWd(R, pp->NBpart_ptr) * pp->vec_e(pp->NBpart_ptr)).x << "\n";

	//if ((pp->Mpart_ptr->m_id == m_options.nrOfParticles[PARTICLETYPE::REAL] - 32)/* and ((-pp->NBpart_ptr->m_mass)*(pp->Mpart_ptr->m_pressure.val / pow(pp->Mpart_ptr->m_density.val, 2) + pp->NBpart_ptr->m_pressure.val / pow(pp->NBpart_ptr->m_density.val, 2))*pp->dWd(R, pp->Mpart_ptr) * pp->vec_e(pp->Mpart_ptr)).x != 0)*/) {
	//	std::cout << "(" << pp->NBpart_ptr->m_id << ")(" << (-pp->NBpart_ptr->m_mass) << "*("<< pp->Mpart_ptr->m_pressure.val<<"/" << pow(pp->Mpart_ptr->m_density.val, 2) << "+" << pp->NBpart_ptr->m_pressure.val << "/" << pow(pp->NBpart_ptr->m_density.val, 2) <<")*" << pp->dWd(R, pp->Mpart_ptr) <<"*"<< pp->vec_e(pp->Mpart_ptr).x << "="<<(-pp->NBpart_ptr->m_mass)*(pp->Mpart_ptr->m_pressure.val / pow(pp->Mpart_ptr->m_density.val, 2) + pp->NBpart_ptr->m_pressure.val / pow(pp->NBpart_ptr->m_density.val, 2))*pp->dWd(R, pp->Mpart_ptr) * pp->vec_e(pp->Mpart_ptr).x << ")+ ";
	//	//std::cout << "(" << pp->NBpart_ptr->m_id << ")(" << (*dvdt)[pp->Mpart_ptr->m_id].x << ")+";
	//}
	//if ((pp->NBpart_ptr->m_id == m_options.nrOfParticles[PARTICLETYPE::REAL] - 32)/* and ((-pp->Mpart_ptr->m_mass)*(pp->NBpart_ptr->m_pressure.val / pow(pp->NBpart_ptr->m_density.val, 2) + pp->Mpart_ptr->m_pressure.val / pow(pp->Mpart_ptr->m_density.val, 2))*pp->dWd(R, pp->NBpart_ptr) * pp->vec_e(pp->NBpart_ptr)).x != 0)*/) {
	//	std::cout << "(" << pp->Mpart_ptr->m_id << ")(" << (-pp->Mpart_ptr->m_mass) << "*(" << pp->NBpart_ptr->m_pressure.val << "/" << pow(pp->NBpart_ptr->m_density.val, 2) << "+" << pp->Mpart_ptr->m_pressure.val << "/" << pow(pp->Mpart_ptr->m_density.val, 2) << ")*" << pp->dWd(R, pp->NBpart_ptr) << "*" << pp->vec_e(pp->NBpart_ptr).x << "=" << ((-pp->Mpart_ptr->m_mass)*(pp->NBpart_ptr->m_pressure.val / pow(pp->NBpart_ptr->m_density.val, 2) + pp->Mpart_ptr->m_pressure.val / pow(pp->Mpart_ptr->m_density.val, 2))*pp->dWd(R, pp->NBpart_ptr) * pp->vec_e(pp->NBpart_ptr)).x << ")+ ";
	//	//std::cout << "(" << pp->Mpart_ptr->m_id << ")(" << (*dvdt)[pp->NBpart_ptr->m_id].x << ")+";
	//}
}

void SPH_CD::ViscosityPartCalculation(std::vector<part_prec_3>* dvdt, ParticlePair * pp) {
	switch (m_options.velocityChangeAlgorithm_viscosityPart) {
	case(0):
		(*dvdt)[pp->Mpart_ptr->m_id] += (static_cast<part_prec>(2.0)*pp->Mpart_ptr->m_DVisc)*pp->Mpart_ptr->dV(pp->NBpart_ptr) / (pp->getDist())*pp->dWd(R, pp->Mpart_ptr);
		(*dvdt)[pp->NBpart_ptr->m_id] += (static_cast<part_prec>(2.0)*pp->NBpart_ptr->m_DVisc)*pp->Mpart_ptr->dV(pp->Mpart_ptr) / (pp->getDist())*pp->dWd(R, pp->NBpart_ptr);
		break;
	case(1):
		(*dvdt)[pp->Mpart_ptr->m_id] += (static_cast<part_prec>(2.0)*pp->Mpart_ptr->m_DVisc)*static_cast<part_prec>(nrOfDim + 2)*glm::dot(pp->Mpart_ptr->dV(pp->NBpart_ptr), pp->vec_e(pp->Mpart_ptr)) / (pp->Mpart_ptr->m_density.val*pp->NBpart_ptr->m_density.val*pp->getDist())*pp->dWd(R, pp->Mpart_ptr) * pp->vec_e(pp->Mpart_ptr);
		(*dvdt)[pp->NBpart_ptr->m_id] += (static_cast<part_prec>(2.0)*pp->NBpart_ptr->m_DVisc)*static_cast<part_prec>(nrOfDim + 2)*glm::dot(pp->Mpart_ptr->dV(pp->Mpart_ptr), pp->vec_e(pp->NBpart_ptr)) / (pp->NBpart_ptr->m_density.val*pp->Mpart_ptr->m_density.val*pp->getDist())*pp->dWd(R, pp->NBpart_ptr) * pp->vec_e(pp->NBpart_ptr);
		break;
	case(2):
		(*dvdt)[pp->Mpart_ptr->m_id] += ((pp->Mpart_ptr->m_DVisc + pp->NBpart_ptr->m_DVisc)) / (static_cast<part_prec>(2.0)*pp->Mpart_ptr->m_density.val*pp->NBpart_ptr->m_density.val*pp->getDist())*(static_cast<part_prec>(nrOfDim + 2)*glm::dot(pp->Mpart_ptr->dV(pp->NBpart_ptr), pp->vec_e(pp->Mpart_ptr))*pp->vec_e(pp->Mpart_ptr) + pp->Mpart_ptr->dV(pp->NBpart_ptr))*pp->dWd(R, pp->Mpart_ptr);
		(*dvdt)[pp->NBpart_ptr->m_id] += ((pp->NBpart_ptr->m_DVisc + pp->Mpart_ptr->m_DVisc)) / (static_cast<part_prec>(2.0)*pp->NBpart_ptr->m_density.val*pp->Mpart_ptr->m_density.val*pp->getDist())*(static_cast<part_prec>(nrOfDim + 2)*glm::dot(pp->NBpart_ptr->dV(pp->Mpart_ptr), pp->vec_e(pp->NBpart_ptr))*pp->vec_e(pp->NBpart_ptr) + pp->NBpart_ptr->dV(pp->Mpart_ptr))*pp->dWd(R, pp->NBpart_ptr);
		break;			
	default:
		break;
	}

	//if ((pp->NBpart_ptr->m_id == m_options.nrOfParticles[PARTICLETYPE::REAL] - 1)) {
	//	std::cout << "Viscosity part (" << pp->Mpart_ptr->m_id << ") = (" << (((pp->Mpart_ptr->m_DVisc + pp->NBpart_ptr->m_DVisc)) / (2.0*pp->Mpart_ptr->m_density.val*pp->NBpart_ptr->m_density.val*pp->getDist())*(static_cast<part_prec>(nrOfDim + 2)*glm::dot(pp->Mpart_ptr->dV(pp->NBpart_ptr), pp->vec_e(pp->Mpart_ptr))*pp->vec_e(pp->Mpart_ptr) + pp->Mpart_ptr->dV(pp->NBpart_ptr))*pp->dWd(R, pp->Mpart_ptr)).x << ", " << (((pp->Mpart_ptr->m_DVisc + pp->NBpart_ptr->m_DVisc)) / (2.0*pp->Mpart_ptr->m_density.val*pp->NBpart_ptr->m_density.val*pp->getDist())*(static_cast<part_prec>(nrOfDim + 2)*glm::dot(pp->Mpart_ptr->dV(pp->NBpart_ptr), pp->vec_e(pp->Mpart_ptr))*pp->vec_e(pp->Mpart_ptr) + pp->Mpart_ptr->dV(pp->NBpart_ptr))*pp->dWd(R, pp->Mpart_ptr)).y << ")    -    (" << pp->NBpart_ptr->m_id << ") = (" << (((pp->NBpart_ptr->m_DVisc + pp->Mpart_ptr->m_DVisc)) / (2.0*pp->NBpart_ptr->m_density.val*pp->Mpart_ptr->m_density.val*pp->getDist())*(static_cast<part_prec>(nrOfDim + 2)*glm::dot(pp->NBpart_ptr->dV(pp->Mpart_ptr), pp->vec_e(pp->NBpart_ptr))*pp->vec_e(pp->NBpart_ptr) + pp->NBpart_ptr->dV(pp->Mpart_ptr))*pp->dWd(R, pp->NBpart_ptr)).x << ", " << (((pp->NBpart_ptr->m_DVisc + pp->Mpart_ptr->m_DVisc)) / (2.0*pp->NBpart_ptr->m_density.val*pp->Mpart_ptr->m_density.val*pp->getDist())*(static_cast<part_prec>(nrOfDim + 2)*glm::dot(pp->NBpart_ptr->dV(pp->Mpart_ptr), pp->vec_e(pp->NBpart_ptr))*pp->vec_e(pp->NBpart_ptr) + pp->NBpart_ptr->dV(pp->Mpart_ptr))*pp->dWd(R, pp->NBpart_ptr)).y << ")\n";
	//}
	//if ((pp->Mpart_ptr->m_id == m_options.nrOfParticles[PARTICLETYPE::REAL] - 1)) {
	//	std::cout << "Viscosity part (" << pp->NBpart_ptr->m_id << ") = (" << (((pp->NBpart_ptr->m_DVisc + pp->Mpart_ptr->m_DVisc)) / (2.0*pp->NBpart_ptr->m_density.val*pp->Mpart_ptr->m_density.val*pp->getDist())*(static_cast<part_prec>(nrOfDim + 2)*glm::dot(pp->NBpart_ptr->dV(pp->Mpart_ptr), pp->vec_e(pp->NBpart_ptr))*pp->vec_e(pp->NBpart_ptr) + pp->NBpart_ptr->dV(pp->Mpart_ptr))*pp->dWd(R, pp->NBpart_ptr)).x << ", " << (((pp->NBpart_ptr->m_DVisc + pp->Mpart_ptr->m_DVisc)) / (2.0*pp->NBpart_ptr->m_density.val*pp->Mpart_ptr->m_density.val*pp->getDist())*(static_cast<part_prec>(nrOfDim + 2)*glm::dot(pp->NBpart_ptr->dV(pp->Mpart_ptr), pp->vec_e(pp->NBpart_ptr))*pp->vec_e(pp->NBpart_ptr) + pp->NBpart_ptr->dV(pp->Mpart_ptr))*pp->dWd(R, pp->NBpart_ptr)).y << ")    -    (" << pp->Mpart_ptr->m_id << ") = (" << (((pp->Mpart_ptr->m_DVisc + pp->NBpart_ptr->m_DVisc)) / (2.0*pp->Mpart_ptr->m_density.val*pp->NBpart_ptr->m_density.val*pp->getDist())*(static_cast<part_prec>(nrOfDim + 2)*glm::dot(pp->Mpart_ptr->dV(pp->NBpart_ptr), pp->vec_e(pp->Mpart_ptr))*pp->vec_e(pp->Mpart_ptr) + pp->Mpart_ptr->dV(pp->NBpart_ptr))*pp->dWd(R, pp->Mpart_ptr)).x << ", " << (((pp->Mpart_ptr->m_DVisc + pp->NBpart_ptr->m_DVisc)) / (2.0*pp->Mpart_ptr->m_density.val*pp->NBpart_ptr->m_density.val*pp->getDist())*(static_cast<part_prec>(nrOfDim + 2)*glm::dot(pp->Mpart_ptr->dV(pp->NBpart_ptr), pp->vec_e(pp->Mpart_ptr))*pp->vec_e(pp->Mpart_ptr) + pp->Mpart_ptr->dV(pp->NBpart_ptr))*pp->dWd(R, pp->Mpart_ptr)).y << ")\n";
	//}


	//if (((-(pp->Mpart_ptr->m_DVisc + pp->NBpart_ptr->m_DVisc)) / (2.0*pp->Mpart_ptr->m_density.val*pp->NBpart_ptr->m_density.val*pp->getDist())*(static_cast<part_prec>(nrOfDim + 2)*glm::dot(pp->Mpart_ptr->dV(pp->NBpart_ptr), pp->vec_e(pp->Mpart_ptr))*pp->vec_e(pp->Mpart_ptr) + pp->Mpart_ptr->dV(pp->NBpart_ptr))*pp->dWd(R, pp->Mpart_ptr)).x != 0)
	//	std::cout << ((-(pp->Mpart_ptr->m_DVisc + pp->NBpart_ptr->m_DVisc)) / (2.0*pp->Mpart_ptr->m_density.val*pp->NBpart_ptr->m_density.val*pp->getDist())*(static_cast<part_prec>(nrOfDim + 2)*glm::dot(pp->Mpart_ptr->dV(pp->NBpart_ptr), pp->vec_e(pp->Mpart_ptr))*pp->vec_e(pp->Mpart_ptr) + pp->Mpart_ptr->dV(pp->NBpart_ptr))*pp->dWd(R, pp->Mpart_ptr)).x << "   " << ((-(pp->NBpart_ptr->m_DVisc + pp->Mpart_ptr->m_DVisc)) / (2.0*pp->NBpart_ptr->m_density.val*pp->Mpart_ptr->m_density.val*pp->getDist())*(static_cast<part_prec>(nrOfDim + 2)*glm::dot(pp->NBpart_ptr->dV(pp->Mpart_ptr), pp->vec_e(pp->NBpart_ptr))*pp->vec_e(pp->NBpart_ptr) + pp->NBpart_ptr->dV(pp->Mpart_ptr))*pp->dWd(R, pp->NBpart_ptr)).x << "\n\n";


	//if ((pp->Mpart_ptr->m_id == m_options.nrOfParticles[PARTICLETYPE::REAL] - 1) and ((((pp->Mpart_ptr->m_DVisc + pp->NBpart_ptr->m_DVisc)) / (2.0*pp->Mpart_ptr->m_density.val*pp->NBpart_ptr->m_density.val*pp->getDist())*(static_cast<part_prec>(nrOfDim + 2)*glm::dot(pp->Mpart_ptr->dV(pp->NBpart_ptr), pp->vec_e(pp->Mpart_ptr))*pp->vec_e(pp->Mpart_ptr) + pp->Mpart_ptr->dV(pp->NBpart_ptr))*pp->dWd(R, pp->Mpart_ptr)).x != 0)) {
	//	//std::cout << "(" << pp->NBpart_ptr->m_id << ")(" << (((pp->Mpart_ptr->m_DVisc + pp->NBpart_ptr->m_DVisc)) / (2.0*pp->Mpart_ptr->m_density.val*pp->NBpart_ptr->m_density.val*pp->getDist())*(static_cast<part_prec>(nrOfDim + 2)*glm::dot(pp->Mpart_ptr->dV(pp->NBpart_ptr), pp->vec_e(pp->Mpart_ptr))*pp->vec_e(pp->Mpart_ptr) + pp->Mpart_ptr->dV(pp->NBpart_ptr))*pp->dWd(R, pp->Mpart_ptr)).x << ")+ ";
	//	std::cout << "(" << pp->NBpart_ptr->m_id << ")(" << (*dvdt)[pp->Mpart_ptr->m_id].x << ")+";
	//}
	//if ((pp->NBpart_ptr->m_id == m_options.nrOfParticles[PARTICLETYPE::REAL] - 1) and ((((pp->NBpart_ptr->m_DVisc + pp->Mpart_ptr->m_DVisc)) / (2.0*pp->NBpart_ptr->m_density.val*pp->Mpart_ptr->m_density.val*pp->getDist())*(static_cast<part_prec>(nrOfDim + 2)*glm::dot(pp->NBpart_ptr->dV(pp->Mpart_ptr), pp->vec_e(pp->NBpart_ptr))*pp->vec_e(pp->NBpart_ptr) + pp->NBpart_ptr->dV(pp->Mpart_ptr))*pp->dWd(R, pp->NBpart_ptr)).x != 0)) {
	//	//std::cout << "(" << pp->Mpart_ptr->m_id << ")(" << (((pp->NBpart_ptr->m_DVisc + pp->Mpart_ptr->m_DVisc)) / (2.0*pp->NBpart_ptr->m_density.val*pp->Mpart_ptr->m_density.val*pp->getDist())*(static_cast<part_prec>(nrOfDim + 2)*glm::dot(pp->NBpart_ptr->dV(pp->Mpart_ptr), pp->vec_e(pp->NBpart_ptr))*pp->vec_e(pp->NBpart_ptr) + pp->NBpart_ptr->dV(pp->Mpart_ptr))*pp->dWd(R, pp->NBpart_ptr)).x << ")+ ";
	//	std::cout << "(" << pp->Mpart_ptr->m_id << ")(" << (*dvdt)[pp->NBpart_ptr->m_id].x << ")+";
	//}
}

void SPH_CD::RenormPressurePartCalculation(std::vector<part_prec_3>* dvdt, ParticlePair * pp, std::vector<part_prec>* gamma) {
	switch (m_options.velocityChangeAlgorithm_pressurePart) {
	case(0):
		// k=0
		(*dvdt)[pp->Mpart_ptr->m_id] += (-pp->NBpart_ptr->m_mass) / (pp->Mpart_ptr->m_density.val*pp->NBpart_ptr->m_density.val)*(pp->Mpart_ptr->m_pressure.val / (*gamma)[pp->Mpart_ptr->m_id] + pp->NBpart_ptr->m_pressure.val / (*gamma)[pp->NBpart_ptr->m_id])*pp->dWd(R, pp->Mpart_ptr) * pp->vec_e(pp->Mpart_ptr);
		(*dvdt)[pp->Mpart_ptr->m_id] += (-pp->Mpart_ptr->m_mass) / (pp->NBpart_ptr->m_density.val*pp->Mpart_ptr->m_density.val)*(pp->NBpart_ptr->m_pressure.val / (*gamma)[pp->NBpart_ptr->m_id] + pp->Mpart_ptr->m_pressure.val / (*gamma)[pp->Mpart_ptr->m_id])*pp->dWd(R, pp->NBpart_ptr) * pp->vec_e(pp->NBpart_ptr);
		break;
	case(1):
		// k=1
		(*dvdt)[pp->Mpart_ptr->m_id] += (-pp->NBpart_ptr->m_mass)*(pp->Mpart_ptr->m_pressure.val / pow(pp->Mpart_ptr->m_density.val, 2) / (*gamma)[pp->Mpart_ptr->m_id] + pp->NBpart_ptr->m_pressure.val / pow(pp->NBpart_ptr->m_density.val, 2) / (*gamma)[pp->NBpart_ptr->m_id])*pp->dWd(R, pp->Mpart_ptr) * pp->vec_e(pp->Mpart_ptr);
		(*dvdt)[pp->NBpart_ptr->m_id] += (-pp->Mpart_ptr->m_mass)*(pp->NBpart_ptr->m_pressure.val / pow(pp->NBpart_ptr->m_density.val, 2) / (*gamma)[pp->NBpart_ptr->m_id] + pp->Mpart_ptr->m_pressure.val / pow(pp->Mpart_ptr->m_density.val, 2) / (*gamma)[pp->Mpart_ptr->m_id])*pp->dWd(R, pp->NBpart_ptr) * pp->vec_e(pp->NBpart_ptr);
		break;
	case(2):
		// k=2
		(*dvdt)[pp->Mpart_ptr->m_id] += (-pp->NBpart_ptr->m_mass)*(pp->Mpart_ptr->m_pressure.val*pp->NBpart_ptr->m_density.val / pow(pp->Mpart_ptr->m_density.val, 3) / (*gamma)[pp->Mpart_ptr->m_id] + pp->NBpart_ptr->m_pressure.val*pp->Mpart_ptr->m_density.val / pow(pp->NBpart_ptr->m_density.val, 3) / (*gamma)[pp->NBpart_ptr->m_id])*pp->dWd(R, pp->Mpart_ptr) * pp->vec_e(pp->Mpart_ptr);
		(*dvdt)[pp->NBpart_ptr->m_id] += (-pp->Mpart_ptr->m_mass)*(pp->NBpart_ptr->m_pressure.val*pp->Mpart_ptr->m_density.val / pow(pp->NBpart_ptr->m_density.val, 3) / (*gamma)[pp->NBpart_ptr->m_id] + pp->Mpart_ptr->m_pressure.val*pp->NBpart_ptr->m_density.val / pow(pp->Mpart_ptr->m_density.val, 3) / (*gamma)[pp->Mpart_ptr->m_id])*pp->dWd(R, pp->NBpart_ptr) * pp->vec_e(pp->NBpart_ptr);
		break;
	default:
		break;
	}
	
}


void SPH_CD::DensityRecalculation() {
	std::vector<part_prec> drhodt;
	drhodt.resize(SPH.Particles.size(), 0.f);
	std::vector<part_prec> gamma;
	gamma.resize(SPH.Particles.size(), 0.f);
	std::vector<part_prec_3> gamma_deriv;
	gamma_deriv.resize(SPH.Particles.size(), { 0.0,0.0,0.0 });
	if (m_options.boundary_handling == RENORMALIZATION) {
		for (auto* i : SPH.ParticlePairs) {
			if (PPL_BOTH_PARTICLES_ARE(REAL, i)) {
				DensityCalculation(&drhodt, i);
				RenormFactorCalc(&gamma, i);
				RenormFactorDerivCalc(&gamma_deriv, i);
			}
		}
		for (auto*& i : SPH.Particles) {
			//std::cout << drhodt[i->m_id] << "   " << gamma[i->m_id] <<"   "<< i->m_density.val << "   " << dot(gamma_deriv[i->m_id], i->m_velocity.val)  << "\n";
			i->m_density.dval = drhodt[i->m_id] / gamma[i->m_id] - (i->m_density.val* dot(gamma_deriv[i->m_id], i->m_velocity.val)) / gamma[i->m_id];
		}

	}
	else {
		for (auto* i : SPH.ParticlePairs) {
			if ((PPL_BOTH_PARTICLES_ARE(REAL, i)) or (PPL_MAIN_PARTICLE_IS(VIRTUAL, i) or PPL_NEIGHBOUR_PARTICLE_IS(VIRTUAL, i))) {
				DensityCalculation(&drhodt, i);
			}
		}
		for (auto*& i : SPH.Particles) {
			i->m_density.dval = drhodt[i->m_id];
		}
	}

}

void SPH_CD::VelocityRecalculation() {
	std::vector<part_prec_3> dvdt_press;
	dvdt_press.resize(SPH.Particles.size(), { 0.0,0.0,0.0 });
	std::vector<part_prec_3> dvdt_dvisc;
	dvdt_dvisc.resize(SPH.Particles.size(), { 0.0,0.0,0.0 });
	if (m_options.boundary_handling == RENORMALIZATION) {
		std::vector<part_prec> gamma;
		gamma.resize(SPH.Particles.size(), 0.f);
		std::vector<part_prec_3> gamma_deriv;
		gamma_deriv.resize(SPH.Particles.size(), { 0.0,0.0,0.0 });
		for (auto* i : SPH.ParticlePairs) {
			if (PPL_BOTH_PARTICLES_ARE(REAL, i)) {
				RenormFactorCalc(&gamma, i);
				RenormFactorDerivCalc(&gamma_deriv, i);
			}
		}
		for (auto* i : SPH.ParticlePairs) {
			if (PPL_BOTH_PARTICLES_ARE(REAL, i)) {
				RenormPressurePartCalculation(&dvdt_press, i, &gamma);
				ViscosityPartCalculation(&dvdt_dvisc, i);
			}
		}
		for (auto*& i : SPH.Particles) {
			if(i->m_type == REAL){
				//std::cout << (dvdt)[i->m_id].x << ", " << (dvdt)[i->m_id].y <<"   "<< (gamma)[i->m_id] << "   " <<(gamma_deriv)[i->m_id].x << ", " << (gamma_deriv)[i->m_id].y << "\n";
				i->m_velocity.dval = (dvdt_dvisc)[i->m_id] + (dvdt_press)[i->m_id] + i->m_pressure.val/i->m_density.val*(gamma_deriv)[i->m_id] / (gamma)[i->m_id];
			}
		}
	}
	else{
		for (auto* i : SPH.ParticlePairs) {
			if ((PPL_BOTH_PARTICLES_ARE(REAL, i)) or (PPL_MAIN_PARTICLE_IS(VIRTUAL, i) or PPL_NEIGHBOUR_PARTICLE_IS(VIRTUAL, i))) {
				PressurePartCalculation(&dvdt_press, i);
			}
			if ((PPL_BOTH_PARTICLES_ARE(REAL, i)) or (PPL_MAIN_PARTICLE_IS(VIRTUAL, i) or PPL_NEIGHBOUR_PARTICLE_IS(VIRTUAL, i))) {
				ViscosityPartCalculation(&dvdt_dvisc, i);
			}
		}
		for (auto*& i : SPH.Particles) { 
			i->m_velocity.dval = (dvdt_dvisc)[i->m_id] + (dvdt_press)[i->m_id];
		}
	}
}

void SPH_CD::DensityAndVelocityRecalculation(){
	std::vector<part_prec> drhodt;
	drhodt.resize(SPH.Particles.size(), 0.f);
	std::vector<part_prec_3> dvdt_press;
	dvdt_press.resize(SPH.Particles.size(), { 0.0,0.0,0.0 });
	std::vector<part_prec_3> dvdt_dvisc;
	dvdt_dvisc.resize(SPH.Particles.size(), { 0.0,0.0,0.0 });
	std::vector<part_prec> gamma;
	std::vector<part_prec_3> gamma_deriv;
	std::vector<part_prec> dens_diff;
	//std::cout << " = " << drhodt[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1] <<" "<< SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1]->m_density.dval<< "\n";

	double ksi = 0.2;
	if (m_options.Density_Diffusion_term) {
		dens_diff.resize(SPH.Particles.size(), 0.f);
	}

	if (m_options.boundary_handling == RENORMALIZATION) {
		gamma.resize(SPH.Particles.size(), 0.f);
		gamma_deriv.resize(SPH.Particles.size(), { 0.0,0.0,0.0 });
		std::cout << "Renormalization\n";
		for (auto* i : SPH.ParticlePairs) {
			if (PPL_BOTH_PARTICLES_ARE(REAL, i)) {
				DensityCalculation(&drhodt, i);
				RenormFactorCalc(&gamma, i);
				RenormFactorDerivCalc(&gamma_deriv, i);
				if (m_options.Density_Diffusion_term) {
					DensityDiffusiveTerm(&dens_diff, i);
				}
			}
		}
		for (auto* i : SPH.ParticlePairs) {
			if (PPL_BOTH_PARTICLES_ARE(REAL, i)) {
				RenormPressurePartCalculation(&dvdt_press, i, &gamma);
				ViscosityPartCalculation(&dvdt_dvisc, i);
			}
		}
		for (auto*& i : SPH.Particles) {
			i->m_density.dval = drhodt[i->m_id] / gamma[i->m_id] - (i->m_density.val* dot(gamma_deriv[i->m_id], i->m_velocity.val)) / gamma[i->m_id];
			if (m_options.Density_Diffusion_term) {	i->m_density.dval += dens_diff[i->m_id]; }
			i->m_velocity.dval = (dvdt_dvisc)[i->m_id] + (dvdt_press)[i->m_id] + i->m_pressure.val / i->m_density.val*(gamma_deriv)[i->m_id] / (gamma)[i->m_id];
		}
	}
	else {
		//std::vector<part_prec> gamma;
		//gamma.resize(SPH.Particles.size(), 0.f);
		//std::cout << drhodt[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1] << " === " << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1]->m_density.dval << " -> ";
		//std::cout << dvdt[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1].x << " === " << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1]->m_velocity.dval.x << " -> ";

		//gamma.resize(SPH.Particles.size(), 0.f);
		for (auto* i : SPH.ParticlePairs) {
			if ((PPL_BOTH_PARTICLES_ARE(REAL, i))) {
			//if (PPL_BOTH_PARTICLES_ARE(REAL, i)) {
				DensityCalculation(&drhodt, i);
				//RenormFactorCalc(&gamma, i);
				if (m_options.Density_Diffusion_term) {
					DensityDiffusiveTerm(&dens_diff, i);
				}
			}
			if ((PPL_BOTH_PARTICLES_ARE(REAL, i)) or (PPL_MAIN_PARTICLE_IS(VIRTUAL, i) or PPL_NEIGHBOUR_PARTICLE_IS(VIRTUAL, i))) {
				PressurePartCalculation(&dvdt_press, i);
			}
			if ((PPL_BOTH_PARTICLES_ARE(REAL, i)) or (PPL_MAIN_PARTICLE_IS(VIRTUAL, i) or PPL_NEIGHBOUR_PARTICLE_IS(VIRTUAL, i))) {
				ViscosityPartCalculation(&dvdt_dvisc, i);
			}

		}
		for (auto*& i : SPH.Particles) {
			//if(i->m_type == REAL)   std::cout << gamma[i->m_id] << "\n";
			if (i->m_type == REAL) {
				i->m_density.dval = drhodt[i->m_id];
				if (m_options.Density_Diffusion_term) { i->m_density.dval += dens_diff[i->m_id]; }
				i->m_velocity.dval = (dvdt_dvisc)[i->m_id] + (dvdt_press)[i->m_id];
			}
		}
	}



	//std::cout << "id" << " " << "\t  drho" << " " << "\t  dv_visc.x" << ", " << "\t  dv_visc.y" << " " << "\t  dv_press.x" << ", " << "\t  dv_press.y" << "\n";
	//for (auto*& i : SPH.Particles) {
	//	if(i->m_type == REAL){
	//		if (i->m_id > (m_options.nrOfParticles[PARTICLETYPE::REAL] - 32)) {
	//			std::cout << i->m_id << " " << drhodt[i->m_id] << " " << (dvdt_dvisc)[i->m_id].x << ", " << (dvdt_dvisc)[i->m_id].y << " " << (dvdt_press)[i->m_id].x << ", " << (dvdt_press)[i->m_id].y << "\n";
	//		}
	//	}
	//}
	//std::cout << " = " << drhodt[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1] << " " << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1]->m_density.dval << "\n";
	//std::cout << " = " << drhodt[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1] << " == " << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1]->m_density.dval << "\n";
	//std::cout << " = " << dvdt[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1].x << " == " << SPH.Particles[m_options.nrOfParticles[PARTICLETYPE::REAL] - 1]->m_velocity.dval.x << "\n";
}


void SPH_CD::RecalcPrep() {
	if (m_options.boundary_handling == RENORMALIZATION) {
		neighbourSearch();
	}
	else {
		neighbourSearch();
		creatingVirtualParticles();
		addingVirtualTOPairs();
	}

}




int SPH_CD::getNrOfParticles(PARTICLETYPE type)
{
	return m_options.nrOfParticles[type];
}

float SPH_CD::getGlobalStats(int number)
{
	switch (number)	{
	case(0):
		return static_cast<float>(getNrOfParticles(REAL));
		break;
	case(1):
		return static_cast<float>(getNrOfParticles(BOUNDARY));
		break;
	case(2):
		return static_cast<float>(getNrOfParticles(VIRTUAL));
		break;
	default:
		return 0;
		break;
	}

}

std::string SPH_CD::getLocalStats(int number) { 
	std::stringstream stats;
	if((number < m_options.nrOfParticles[PARTICLETYPE::REAL] + m_options.nrOfParticles[PARTICLETYPE::BOUNDARY] + m_options.nrOfParticles[PARTICLETYPE::VIRTUAL])and(number>=0)){
		stats << "Paricle " << std::to_string(SPH.Particles[number]->m_id) << " is ";
		if (SPH.Particles[number]->m_type == REAL) { stats << "real "; }
		if (SPH.Particles[number]->m_type == BOUNDARY) { stats << "boundary "; }
		if (SPH.Particles[number]->m_type == VIRTUAL) { stats << "virtual "; }
		stats << "\n";
		stats << "Position:{";
		if (nrOfDim > 0) { stats << std::to_string(SPH.Particles[number]->m_position.val.x); }
		if (nrOfDim > 1) { stats << ", " << std::to_string(SPH.Particles[number]->m_position.val.y); }
		if (nrOfDim > 2) { stats << ", " << std::to_string(SPH.Particles[number]->m_position.val.z); }
		stats << "} , Velocity:{";
		if (nrOfDim > 0) { stats << std::to_string(SPH.Particles[number]->m_velocity.val.x); }
		if (nrOfDim > 1) { stats << ", " << std::to_string(SPH.Particles[number]->m_velocity.val.y); }
		if (nrOfDim > 2) { stats << ", " << std::to_string(SPH.Particles[number]->m_velocity.val.z); }
		stats << "}\n";
		stats << "Smoothing radius:" << std::to_string(SPH.Particles[number]->m_SmR) << ", Density:" << std::to_string(SPH.Particles[number]->m_density.val) << ", Mass:" << std::to_string(SPH.Particles[number]->m_mass) << "\n";
		stats << "Number of neighboirs: " << std::to_string(SPH.Particles[number]->nrOfNeighbours) << "\n";
	}
	else {
		stats << "There's no components with such number.\n";
	}
	return stats.str().c_str();
}
