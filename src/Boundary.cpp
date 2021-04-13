#include "Boundary.h"




bool operator==(const �D_Boundary& CDB1, const �D_Boundary& CDB2) {
	return(CDB1.B_mesh == CDB2.B_mesh);
}
bool operator!=(const �D_Boundary& CDB1, const �D_Boundary& CDB2) {
	return !(CDB1.B_mesh == CDB2.B_mesh);
}

void ScanForPeriodicBC(std::vector<�D_Boundary>* arrayOfBC) {

	Mesh* main_mesh;
	Mesh* periodic_mesh;
	bool pair_found;
	for (�D_Boundary i : *arrayOfBC) {
		pair_found = true;
		if (i.isPeriodic()) {
			main_mesh = i.ReturnMesh();
			periodic_mesh = i.ReturnMeshToLoop();
			pair_found = false;
		}
		for (�D_Boundary j : *arrayOfBC) {
			if (i != j) {
				if (j.isPeriodic()) {
					if ((main_mesh == j.ReturnMesh()) and (periodic_mesh == j.ReturnMeshToLoop())) {
						pair_found = true;
						break;
					}
				}
			}
		}
		if (pair_found == false) {
			assert("�D_Boundary::ScanForPeriodicBC boundary pair wasn't found" && 0);
		}
		else {
			if (&main_mesh == &periodic_mesh) {

			}
			else {
				assert("�D_Boundary::ScanForPeriodicBC meshes are not similar" && 0);
			}
		}
	}

}

