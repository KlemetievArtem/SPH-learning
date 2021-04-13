#pragma once

#include<iostream>
#include<vector>


#include "Vertex.h"
#include "Primitives.h"
#include "Shader.h"
#include "Texture.h"
#include "Material.h"


class Mesh {
private:
	Vertex* vertexArray;
	unsigned nrOfVertices;
	GLuint* indexArray;
	unsigned nrOfIndices;

	//std::vector<Vertex> vertices;
	//std::vector<GLuint> indices;

	GLuint VAO;
	GLuint VBO;
	GLuint EBO;

	glm::vec3 position;
	glm::vec3 origin;
	glm::vec3 rotation;
	glm::vec3 scale;
	glm::mat4 ModelMatrix;


	void initVAO() {


		//VAO,VBO,EBO
		//GEN VAO AND BIND
		glCreateVertexArrays(1, &this->VAO);
		glBindVertexArray(this->VAO);


		
			
		//GEN VBO AND BIND AND SEND DATA
		glGenBuffers(1, &this->VBO);
		glBindBuffer(GL_ARRAY_BUFFER, this->VBO);
		glBufferData(GL_ARRAY_BUFFER, this->nrOfVertices *sizeof(Vertex), this->vertexArray, GL_STREAM_DRAW);

		//GEN EBO AND BIND AND SEND DATA
		if (this->nrOfIndices > 0) {
			glGenBuffers(1, &this->EBO);
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->EBO);
			glBufferData(GL_ELEMENT_ARRAY_BUFFER, this->nrOfIndices * sizeof(GLuint), this->indexArray, GL_STREAM_DRAW);
		}



		//SET VERTEXATTRIBUTEPOINTERS AND ENABLE
		//GLuint attribLoc = glGetAttribLocation(core_program, "vertex_position");
		//Position
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (GLvoid*)offsetof(Vertex, position));
		glEnableVertexAttribArray(0);
		//Color
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (GLvoid*)offsetof(Vertex, color));
		glEnableVertexAttribArray(1);
		//Texcoord
		glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (GLvoid*)offsetof(Vertex, texcoord));
		glEnableVertexAttribArray(2);
		//Normal
		glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (GLvoid*)offsetof(Vertex, normal));
		glEnableVertexAttribArray(3);

		//BIND VAO 0
	}
	/*

	void initVAO(Vertex* VertexArray, const unsigned& nrOfVertices,
		GLuint* IndexArray, const unsigned& nrOfIndices) {

		//Set variables
		this->nrOfVertices = nrOfVertices;
		this->nrOfIndices = nrOfIndices;
		//VAO,VBO,EBO
		//GEN VAO AND BIND
		glCreateVertexArrays(1, &this->VAO);
		glBindVertexArray(this->VAO);

		//GEN VBO AND BIND AND SEND DATA
		glGenBuffers(1, &this->VBO);
		glBindBuffer(GL_ARRAY_BUFFER, this->VBO);
		glBufferData(GL_ARRAY_BUFFER, this->nrOfVertices * sizeof(Vertex), VertexArray, GL_STATIC_DRAW);

		//GEN EBO AND BIND AND SEND DATA

		if (this->nrOfIndices > 0) {
			glGenBuffers(1, &this->EBO);
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->EBO);
			glBufferData(GL_ELEMENT_ARRAY_BUFFER, this->nrOfIndices * sizeof(GLuint), IndexArray, GL_STATIC_DRAW);
		}

		//SET VERTEXATTRIBUTEPOINTERS AND ENABLE
		//GLuint attribLoc = glGetAttribLocation(core_program, "vertex_position");
		//Position
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (GLvoid*)offsetof(Vertex, position));
		glEnableVertexAttribArray(0);
		//Color
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (GLvoid*)offsetof(Vertex, color));
		glEnableVertexAttribArray(1);
		//Texcoord
		glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (GLvoid*)offsetof(Vertex, texcoord));
		glEnableVertexAttribArray(2);
		//Normal
		glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (GLvoid*)offsetof(Vertex, normal));
		glEnableVertexAttribArray(3);

		//BIND VAO 0
		glBindVertexArray(0);
	}

	*/

	void updateUniforms(Shader* shader) {
		shader->use();
		shader->setMat4fv(this->ModelMatrix, "ModelMatrix");

	}
	void updaeteModelMatrix() {
		this->ModelMatrix = glm::mat4(1.0f);
		this->ModelMatrix = glm::translate(this->ModelMatrix, this->origin);
		this->ModelMatrix = glm::rotate(this->ModelMatrix, glm::radians(this->rotation.x), glm::vec3(1.0f, 0.0f, 0.0f));
		this->ModelMatrix = glm::rotate(this->ModelMatrix, glm::radians(this->rotation.y), glm::vec3(0.0f, 1.0f, 0.0f));
		this->ModelMatrix = glm::rotate(this->ModelMatrix, glm::radians(this->rotation.z), glm::vec3(0.0f, 0.0f, 1.0f));
		this->ModelMatrix = glm::translate(this->ModelMatrix, this->position - this->origin);
		this->ModelMatrix = glm::scale(this->ModelMatrix, this->scale);
	}


public:
	Mesh(Vertex* VertexArray, const unsigned& nrOfVertices, GLuint* IndexArray, const unsigned& nrOfIndices,
		glm::vec3 position = glm::vec3(0.f), glm::vec3 origin = glm::vec3(0.f) , glm::vec3 rotation = glm::vec3(0.f), glm::vec3 scale = glm::vec3(1.f)) {
		this->position = position;
		this->origin = origin;
		this->rotation = rotation;
		this->scale = scale;

		this->nrOfVertices = nrOfVertices;
		this->nrOfIndices = nrOfIndices;

		this->vertexArray = new Vertex[this->nrOfVertices];
		for (size_t i = 0; i < nrOfVertices; i++) {
			this->vertexArray[i] = VertexArray[i];
		}
		this->indexArray = new GLuint[this->nrOfIndices];
		for (size_t i = 0; i < nrOfIndices; i++) {
			this->indexArray[i] = IndexArray[i];
		}

		this->initVAO();
		this->updaeteModelMatrix();
	}
	Mesh(Primitive* primitive,
		glm::vec3 position = glm::vec3(0.f), glm::vec3 origin = glm::vec3(0.f), glm::vec3 rotation = glm::vec3(0.f), glm::vec3 scale = glm::vec3(1.f)) {
		this->position = position;
		this->origin = origin;
		this->rotation = rotation;
		this->scale = scale;

		this->nrOfVertices = primitive->getNrOfVertices();
		this->nrOfIndices = primitive->getNrOfIndices();

		this->vertexArray = new Vertex[this->nrOfVertices];
		for (size_t i = 0; i < this->nrOfVertices; i++) {
			this->vertexArray[i] = primitive->getVertices()[i];
		}
		this->indexArray = new GLuint[this->nrOfIndices];
		for (size_t i = 0; i < this->nrOfIndices; i++) {
			this->indexArray[i] = primitive->getIndices()[i];
		}

		this->initVAO();
		this->updaeteModelMatrix();
	}
	Mesh(const Mesh& obj) {
		this->position = obj.position;
		this->origin = obj.origin;
		this->rotation = obj.rotation;
		this->scale = obj.scale;

		this->nrOfVertices = obj.nrOfVertices;
		this->nrOfIndices = obj.nrOfIndices;

		this->vertexArray = new Vertex[this->nrOfVertices];
		for (size_t i = 0; i < this->nrOfVertices; i++) {
			this->vertexArray[i] = obj.vertexArray[i];
		}
		this->indexArray = new GLuint[this->nrOfIndices];
		for (size_t i = 0; i < this->nrOfIndices; i++) {
			this->indexArray[i] = obj.indexArray[i];
		}

		this->initVAO();
		this->updaeteModelMatrix();
	}
	Mesh(Mesh* newMesh) {
		this->position = newMesh->position;
		this->origin = newMesh->origin;
		this->rotation = newMesh->rotation;
		this->scale = newMesh->scale;

		this->nrOfVertices = newMesh->nrOfVertices;
		this->nrOfIndices = newMesh->nrOfIndices;

		delete[] this->vertexArray;
		this->vertexArray = new Vertex[this->nrOfVertices];
		for (size_t i = 0; i < this->nrOfVertices; i++) {
			this->vertexArray[i] = newMesh->vertexArray[i];
		}
		delete[] this->indexArray;
		this->indexArray = new GLuint[this->nrOfIndices];
		for (size_t i = 0; i < this->nrOfIndices; i++) {
			this->indexArray[i] = newMesh->indexArray[i];
		}

		this->initVAO();
		this->updaeteModelMatrix();
	}


	~Mesh() {
		glDeleteVertexArrays(1, &this->VAO);
		glDeleteBuffers(1, &this->VBO);
		if (this->nrOfIndices > 0) {
			glDeleteBuffers(1, &this->EBO);
		}
		delete[] this->vertexArray;
		delete[] this->indexArray;
	}


	//Accessors
	glm::vec3  getPosition() const { return this->position; }
	//Modifires
	void setPosition(const glm::vec3 position) {
		this->position = position;
	}
	void setOrigin(const glm::vec3 origin) {
		this->origin = origin;
	}
	void setRotation(const glm::vec3 rotation) {
		this->rotation = rotation;
	}
	void setScale(const glm::vec3 scale) {
		this->scale = scale;
	}

	//Functions
	void move(const glm::vec3 position) {
		this->position += position;
	}
	void rotate(const glm::vec3 rotation) {
		this->rotation += rotation;
		if (this->rotation.x > 360)
			this->rotation.x -= 360;
		else if(this->rotation.x < 0)
			this->rotation.x += 360;
		if (this->rotation.y > 360)
			this->rotation.y -= 360;
		else if (this->rotation.y < 0)
			this->rotation.y += 360;
		if (this->rotation.z > 360)
			this->rotation.z -= 360;
		else if (this->rotation.z < 0)
			this->rotation.z += 360;
	}
	void scaleUp(const glm::vec3 scale) {
		this->scale += scale;
	}

	void update() {

	}
	void render(Shader* shader) {
		//Update uniforms
		this->updaeteModelMatrix();
		this->updateUniforms(shader);

		//Bind vertex array object
		glBindVertexArray(this->VAO);

		//Draw
		if(this->nrOfIndices==0)
			glDrawArrays(GL_TRIANGLES, 0 , this->nrOfVertices);
		else
		glDrawElements(GL_TRIANGLES, this->nrOfIndices, GL_UNSIGNED_INT, 0);


		//Reset
		glBindVertexArray(0);
		glUseProgram(0);
		glActiveTexture(0);
		glBindTexture(GL_TEXTURE_2D, 0);
	}





	void changeColorTo(glm::vec3 newcolor = glm::vec3(0.f)) {
		for (size_t i = 0; i < nrOfVertices; i++) {
			this->vertexArray[i].color = newcolor;
		}
		//initVAO();   //нвемэ ме ясоеп
	}

	void printAll() {
		std::cout << "\n";
		std::cout << this<<"\n";
		std::cout << "vertexArray["<< nrOfVertices <<"] = {\n";
		std::cout << "position     " << "color\n";
		for (size_t i = 0; i < nrOfVertices; i++) {
			std::cout <<"{"<< vertexArray[i].position.x << "," << vertexArray[i].position.y << "," << vertexArray[i].position.z << "},{" << vertexArray[i].color.x << "," << vertexArray[i].color.y << "," << vertexArray[i].color.z << "}\n";
		}
		std::cout << "}\n";
		std::cout << "IndexArray[" << nrOfIndices << "] = {\n";
		for (size_t i = 0; i < nrOfIndices; i++) {
			std::cout << indexArray[i] << ", ";
		}
		std::cout << "}\n";

		std::cout << "Mesh.position{" << position.x << ", " << position.y << ", " << position.z << "}\n";
		std::cout << "Mesh.origin{" << origin.x << ", " << origin.y << ", " << origin.z << "}\n";
		std::cout << "Mesh.rotation{" << rotation.x << ", " << rotation.y << ", " << rotation.z << "}\n";
		std::cout << "Mesh.scale{" << scale.x << ", " << scale.y << ", " << scale.z << "}\n";
	}



	int getNumberOfIndice() { return nrOfIndices; }
	int getNrOfVertices() { return nrOfVertices; }

	GLuint* getIndexArray() { return indexArray; }
	Vertex* getVertexArray() { return vertexArray; }



	friend bool operator==(const Mesh& mesh1, const Mesh& mesh2);


};