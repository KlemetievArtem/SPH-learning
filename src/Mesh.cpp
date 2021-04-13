#include <Mesh.h>

bool operator==(const Mesh& mesh1, const Mesh& mesh2) {
	return(mesh1.position == mesh2.position &&
		mesh1.origin == mesh2.origin &&
		mesh1.rotation == mesh2.rotation &&
		mesh1.scale == mesh2.scale &&
		mesh1.nrOfVertices == mesh2.nrOfVertices// &&
		//mesh1.vertexArray->position == mesh2.vertexArray->position &&
		//mesh1.vertexArray->normal == mesh2.vertexArray->normal &&
		//mesh1.nrOfIndices == mesh2.nrOfIndices &&
		//mesh1.indexArray == mesh2.indexArray
		);
}
