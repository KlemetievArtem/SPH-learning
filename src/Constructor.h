#pragma once
#include <sstream>
#include <iostream>


class ToolConstructor {
private:
	std::vector<glm::vec3*> vertices_memory;
public:
	ToolConstructor() {}
	~ToolConstructor() {
		for (size_t i = 0;i < this->vertices_memory.size();i++)
			delete this->vertices_memory[i];
	}
	void remember(float x, float y, float z) {
		std::cout << "point: " << x << "," << y << "," << z << "was saved. ";
		vertices_memory.push_back(new glm::vec3(x, y, z));
		std::cout << vertices_memory.size() << " vertices saved.\n";
	}

	void addTriangleTo(std::vector<Mesh*>* meshes, std::vector<Mesh*>* cursor_mesh) {
		int vsize = vertices_memory.size();
		if (vertices_memory.size() < 3)
			std::cout << "NOT_ENOUGH_VERTICES\n";
		meshes->push_back(new Mesh(&Triangle(*vertices_memory[vsize-3], *vertices_memory[vsize-2], *vertices_memory[vsize-1]), cursor_mesh->front()->getPosition(), glm::vec3(0.f), glm::vec3(0.f), glm::vec3(1.f)));
		std::cout << cursor_mesh->front()->getPosition().x << " " << cursor_mesh->front()->getPosition().y << " " << cursor_mesh->front()->getPosition().z<< "\n";
		/*
		for (size_t i = 0;i < this->vertices_memory.size();i++)
			delete this->vertices_memory[i];
		vertices_memory.resize(0);
		*/
	}


	inline int checkSize() {
		return vertices_memory.size();
	}

	void savePrimitive(const std::string& filePath, const std::string& className) {	


	std::ofstream in_file;
	in_file.open(filePath, std::ios::app);

	std::stringstream ss;
	ss << "#include \"Primitives.h\"\n";
	ss << "class "<< className<<" : public Primitive { \n";
	ss << "public:\n";
	ss <<	 className << "(): Primitive() { \n";
	ss <<		"Vertex vertices[]={ \n";
	int triCount=0;
	glm::vec3 p1;
	glm::vec3 p2;
	glm::vec3 p3;
	for (int i = 0;i < this->checkSize();i++) {
		if ((i)%3==0){
			p1 = *vertices_memory[triCount * 3];
			p2 = *vertices_memory[triCount * 3 + 1];
			p3 = *vertices_memory[triCount * 3 + 2];
			++triCount;
		}
		glm::vec3 normal;
		normal.x = (p2.y - p1.y)*(p3.z - p1.z) - (p2.z - p1.z)*(p3.y - p1.y);
		normal.y = (p2.z - p1.z)*(p3.x - p1.x) - (p2.x - p1.x)*(p3.z - p1.z);
		normal.z = (p2.x - p1.x)*(p3.y - p1.y) - (p2.y - p1.y)*(p3.x - p1.x);
		normal = glm::normalize(normal);

		ss <<		 "glm::vec3(" << vertices_memory[i]->x << "," << vertices_memory[i]->y << "," << vertices_memory[i]->z  << "), glm::vec3(1.f,1.f,1.f), glm::vec2(1.f,1.f), glm::vec3(" << normal.x << "," << normal.y << "," << normal.z <<")";
		if (i != this->checkSize() - 1)
			ss << ",";
		ss << "\n";

	}
	ss <<		"};\n";
	ss <<		"unsigned nrOfVertices = sizeof(vertices) / sizeof(Vertex);\n";
	ss <<		"GLuint indices[]={ \n";
	for (int i = 0;i < this->checkSize();i++) {
		ss		 << i;
		if (i != this->checkSize() - 1)
			ss << ",";
		ss << "\n";
	}
	ss <<		"};\n";
	ss <<		"unsigned nrOfIndices = sizeof(indices) / sizeof(GLuint);\n";
	ss <<		"this->set(vertices, nrOfVertices, indices, nrOfIndices);\n";
	ss <<	"}\n";
	ss << "};\n";
	std::cout << ss.str();

	in_file << ss.str();
	in_file.close(); // закрываем файл

	std::cout << "Primitive was saved\n";
	
	}

		

};