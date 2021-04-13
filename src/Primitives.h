#pragma once

#include <GL\glew.h>
#include <GLFW\glfw3.h>

#include <vector>
#include "Vertex.h"

class Primitive {
private:
	std::vector<Vertex> vertices;
	std::vector<GLuint> indices;

public:
	Primitive() {}
	virtual ~Primitive() {}

	//Functions
	void set(const Vertex* vertices, const unsigned nrOfVertices, const GLuint* indices, const unsigned nrOfIndices) {
		for (size_t i = 0;i < nrOfVertices;i++) {
			this->vertices.push_back(vertices[i]);
		}
		for (size_t i = 0;i < nrOfIndices;i++) {
			this->indices.push_back(indices[i]);

		}

	}

	inline Vertex* getVertices() { return this->vertices.data(); }
	inline GLuint* getIndices() { return this->indices.data(); }
	inline const unsigned getNrOfVertices() { return this->vertices.size(); }
	inline const unsigned getNrOfIndices() { return this->indices.size(); }

};
class Point :public Primitive {
public:
	Point() :Primitive() {
		Vertex vertices[] =
		{
			//Position                         //Color                           //Texcoords                      //Normals
			glm::vec3(0.f,  0.f, 0.f),   glm::vec3(1.0f, 0.0f, 0.0f),      glm::vec2(0.0f, 1.0f),          glm::vec3(0.0f, 0.0f, 1.0f),
		};
		unsigned nrOfVertices = sizeof(vertices) / sizeof(Vertex);

		GLuint indices[] = { 0 };
		unsigned nrOfIndices = sizeof(indices) / sizeof(GLuint);

		this->set(vertices, nrOfVertices, indices, nrOfIndices);
	}

};


class Quad : public Primitive {
public:
	enum NORMALAXIS { X, Y, Z };
	enum NORMALDIRECTION {PLUS, MINUS};
	struct QuadNormal {
		NORMALAXIS normalAxis;
		NORMALDIRECTION normalDirection;
		QuadNormal(NORMALAXIS axis, NORMALDIRECTION direction)
			: normalAxis(axis), normalDirection(direction){}
	};

	Quad() :Primitive() {
		Vertex vertices[] =
		{
			//Position                         //Color                           //Texcoords                      //Normals
			glm::vec3(-0.5f,  0.5f, 0.0f),   glm::vec3(1.0f, 0.0f, 0.0f),      glm::vec2(0.0f, 1.0f),          glm::vec3(0.0f, 0.0f, 1.0f),
			glm::vec3(-0.5f, -0.5f, 0.0f),   glm::vec3(0.0f, 1.0f, 0.0f),      glm::vec2(0.0f, 0.0f),		   glm::vec3(0.0f, 0.0f, 1.0f),
			glm::vec3(0.5f,  -0.5f, 0.0f),   glm::vec3(0.0f, 0.0f, 1.0f),      glm::vec2(1.0f, 0.0f),		   glm::vec3(0.0f, 0.0f, 1.0f),
			glm::vec3(0.5f,   0.5f, 0.0f),   glm::vec3(1.0f, 1.0f, 0.0f),      glm::vec2(1.0f, 1.0f),		   glm::vec3(0.0f, 0.0f, 1.0f)
		};
		unsigned nrOfVertices = sizeof(vertices) / sizeof(Vertex);

		GLuint indices[] =
		{
			0, 1, 2,
			0, 2, 3
		};
		unsigned nrOfIndices = sizeof(indices) / sizeof(GLuint);

		this->set(vertices, nrOfVertices, indices, nrOfIndices);
	}

	Quad(glm::vec3 start, QuadNormal normal = QuadNormal(Z, PLUS), float v1 = 1.f, float v2 = 1.f, glm::vec3 color = glm::vec3(1.0f, 1.0f, 1.0f)) :Primitive() {
		glm::vec3 normalvec;
		switch (normal.normalAxis)
		{
		case Quad::X:
			normalvec = glm::vec3(1.0f, 0.0f, 0.0f);
			break;
		case Quad::Y:
			normalvec = glm::vec3(0.0f, 1.0f, 0.0f);
			break;
		case Quad::Z:
			normalvec = glm::vec3(0.0f, 0.0f, 1.0f);
			break;
		default:
			normalvec = glm::vec3(0.0f, 0.0f, 0.0f);
			break;
		}
		if (normal.normalDirection == Quad::MINUS)
			normalvec = -normalvec;

		Vertex vertices[4];
		switch (normal.normalAxis) {
		case Quad::X:
			                  //Position                                     //Color       //Texcoords                      //Normals
			vertices[0] = { glm::vec3(start.x, start.y, start.z + v2),         color,      glm::vec2(0.0f, 1.0f),           normalvec };
			vertices[1] = { glm::vec3(start.x, start.y, start.z),              color,      glm::vec2(0.0f, 0.0f),           normalvec };
			vertices[2] = { glm::vec3(start.x, start.y + v1, start.z),         color,      glm::vec2(0.0f, 0.0f),           normalvec };
			vertices[3] = { glm::vec3(start.x, start.y + v1, start.z + v2),    color,      glm::vec2(0.0f, 0.0f),           normalvec };
			break;
		case Quad::Y:
			                  //Position                                     //Color       //Texcoords                      //Normals
			vertices[0] = { glm::vec3(start.x + v1, start.y, start.z),         color,      glm::vec2(0.0f, 1.0f),           normalvec };
			vertices[1] = { glm::vec3(start.x, start.y, start.z),              color,      glm::vec2(0.0f, 0.0f),           normalvec };
			vertices[2] = { glm::vec3(start.x, start.y, start.z + v2),         color,      glm::vec2(0.0f, 0.0f),           normalvec };
			vertices[3] = { glm::vec3(start.x + v1, start.y, start.z + v2),    color,      glm::vec2(0.0f, 0.0f),           normalvec };
			break;
		case Quad::Z:
			                  //Position                                     //Color       //Texcoords                      //Normals
			vertices[0] = { glm::vec3(start.x, start.y + v2, start.z),         color,      glm::vec2(0.0f, 1.0f),           normalvec };
			vertices[1] = { glm::vec3(start.x, start.y, start.z),              color,      glm::vec2(0.0f, 0.0f),           normalvec };
			vertices[2] = { glm::vec3(start.x + v1, start.y, start.z),         color,      glm::vec2(0.0f, 0.0f),           normalvec };
			vertices[3] = { glm::vec3(start.x + v1, start.y + v2, start.z),    color,      glm::vec2(0.0f, 0.0f),           normalvec };
			break;

		}
		unsigned nrOfVertices = sizeof(vertices) / sizeof(Vertex);
		//std::cout << "PRIMITIVE::QUAD(....)::nrOfVertices = " << nrOfVertices<<"\n";

		GLuint indices[6];
		switch (normal.normalDirection) {
		case Quad::PLUS:
			indices[0] = 0; indices[1] = 1; indices[2] = 2;
			indices[3] = 0; indices[4] = 2; indices[5] = 3;
			break;
		case Quad::MINUS:
			indices[0] = 0; indices[1] = 2; indices[2] = 1;
			indices[3] = 0; indices[4] = 3; indices[5] = 2;
			break;
		}


		unsigned nrOfIndices = sizeof(indices) / sizeof(GLuint);

		this->set(vertices, nrOfVertices, indices, nrOfIndices);
	}

};

class Qube : public Primitive {
public:
	Qube(glm::vec3 color = glm::vec3(0.8f)) : Primitive() {
		Vertex vertices[] = {
			//Position                      //Color           //Texcoords                      //Normals
			glm::vec3(-0.5f,  0.5f, -0.5f),   color,      glm::vec2(0.0f, 1.0f),           glm::vec3(0.0f, 0.0f, 1.0f),
			glm::vec3(-0.5f, -0.5f, -0.5f),   color,      glm::vec2(0.0f, 0.0f),		   glm::vec3(0.0f, 0.0f, 1.0f),
			glm::vec3(0.5f,  -0.5f, -0.5f),   color,      glm::vec2(1.0f, 0.0f),		   glm::vec3(0.0f, 0.0f, 1.0f),
			glm::vec3(0.5f,   0.5f, -0.5f),   color,      glm::vec2(1.0f, 1.0f),		   glm::vec3(0.0f, 0.0f, 1.0f),

			glm::vec3(-0.5f,  0.5f, 0.5f),    color,      glm::vec2(0.0f, 1.0f),           glm::vec3(0.0f, 0.0f, 1.0f),
			glm::vec3(-0.5f, -0.5f, 0.5f),    color,      glm::vec2(0.0f, 0.0f),		   glm::vec3(0.0f, 0.0f, 1.0f),
			glm::vec3(0.5f,  -0.5f, 0.5f),    color,      glm::vec2(1.0f, 0.0f),		   glm::vec3(0.0f, 0.0f, 1.0f),
			glm::vec3(0.5f,   0.5f, 0.5f),    color,      glm::vec2(1.0f, 1.0f),		   glm::vec3(0.0f, 0.0f, 1.0f)

		};
		unsigned nrOfVertices = sizeof(vertices) / sizeof(Vertex);

		GLuint indices[] =
		{
			0, 2, 1,
			0, 3, 2,

			4, 5, 6,
			4, 6, 7,

			0, 1, 5,
			0, 5, 4,

			0, 4, 7,
			0, 7, 3,

			1, 6, 5,
			1, 2, 6,

			3, 7, 6,
			3, 6, 2
		};
		unsigned nrOfIndices = sizeof(indices) / sizeof(GLuint);

		this->set(vertices, nrOfVertices, indices, nrOfIndices);
	}
};





class Triangle : public Primitive {
public:
	Triangle() :Primitive() {
		Vertex vertices[] =
		{
			//Position                         //Color                           //Texcoords                      //Normals
			glm::vec3(-0.5f,  0.5f, 0.0f),   glm::vec3(1.0f, 0.0f, 0.0f),      glm::vec2(0.0f, 1.0f),          glm::vec3(0.0f, 0.0f, 1.0f),
			glm::vec3(-0.5f, -0.5f, 0.0f),   glm::vec3(0.0f, 1.0f, 0.0f),      glm::vec2(0.0f, 0.0f),		   glm::vec3(0.0f, 0.0f, 1.0f),
			glm::vec3(0.5f,  -0.5f, 0.0f),   glm::vec3(0.0f, 0.0f, 1.0f),      glm::vec2(1.0f, 0.0f),		   glm::vec3(0.0f, 0.0f, 1.0f)
		};
		unsigned nrOfVertices = sizeof(vertices) / sizeof(Vertex);

		GLuint indices[] =
		{
			0, 1, 2,
		};
		unsigned nrOfIndices = sizeof(indices) / sizeof(GLuint);

		this->set(vertices, nrOfVertices, indices, nrOfIndices);
	}
	Triangle(glm::vec3 p1, glm::vec3 p2, glm::vec3 p3, glm::vec3 color = glm::vec3(0.0f,0.0f,0.0f)) :Primitive() {
		glm::vec3 normal;
		normal.x = (p2.y - p1.y)*(p3.z - p1.z) - (p2.z - p1.z)*(p3.y - p1.y);
		normal.y = (p2.z - p1.z)*(p3.x - p1.x) - (p2.x - p1.x)*(p3.z - p1.z);
		normal.z = (p2.x - p1.x)*(p3.y - p1.y) - (p2.y - p1.y)*(p3.x - p1.x);
		normal = glm::normalize(normal);
		Vertex vertices[] =
		{
			//Position  //Color      //Texcoords                      //Normals
			p1,           color,      glm::vec2(0.0f, 1.0f),          normal,
			p2,           color,      glm::vec2(0.0f, 0.0f),		   normal,
			p3,           color,      glm::vec2(1.0f, 0.0f),		   normal
		};
		unsigned nrOfVertices = sizeof(vertices) / sizeof(Vertex);

		GLuint indices[] =
		{
			0, 1, 2
		};
		unsigned nrOfIndices = sizeof(indices) / sizeof(GLuint);

		this->set(vertices, nrOfVertices, indices, nrOfIndices);
	}
};

class Pyramid : public Primitive {
public:
	Pyramid() : Primitive() {
		Vertex vertices[] =
		{
			//Position                         //Color                           //Texcoords                      //Normals
			//Triangle front
			glm::vec3(  0.f,  0.5f,  0.f),   glm::vec3(1.0f, 0.0f, 0.0f),      glm::vec2(0.5f, 1.0f),          glm::vec3(0.0f, 0.0f, 1.0f),
			glm::vec3(-0.5f, -0.5f, 0.5f),   glm::vec3(0.0f, 1.0f, 0.0f),      glm::vec2(0.0f, 0.0f),		   glm::vec3(0.0f, 0.0f, 1.0f),
			glm::vec3( 0.5f, -0.5f,  0.5f),   glm::vec3(0.0f, 0.0f, 1.0f),      glm::vec2(1.0f, 0.0f),		   glm::vec3(0.0f, 0.0f, 1.0f),
			//Triangle left
			glm::vec3(  0.f,  0.5f,   0.f),   glm::vec3(1.0f, 1.0f, 0.0f),      glm::vec2(0.5f, 1.0f),         glm::vec3(-1.0f, 0.0f, 0.f),
			glm::vec3(-0.5f, -0.5f, -0.5f),   glm::vec3(0.0f, 1.0f, 1.0f),      glm::vec2(0.0f, 0.0f),		   glm::vec3(-1.0f, 0.0f, 0.f),
			glm::vec3(-0.5f, -0.5f,  0.5f),   glm::vec3(0.0f, 0.0f, 1.0f),      glm::vec2(1.0f, 0.0f),		   glm::vec3(-1.0f, 0.0f, 0.f),
			//Triangle back
			glm::vec3(  0.f,  0.5f,  0.f),    glm::vec3(1.0f, 1.0f, 0.0f),      glm::vec2(0.5f, 1.0f),         glm::vec3(0.0f, 0.0f, -1.f),
			glm::vec3( 0.5f, -0.5f, -0.5f),   glm::vec3(0.0f, 1.0f, 1.0f),      glm::vec2(0.0f, 0.0f),		   glm::vec3(0.0f, 0.0f, -1.f),
			glm::vec3(-0.5f, -0.5f, -0.5f),   glm::vec3(0.0f, 0.0f, 1.0f),      glm::vec2(1.0f, 0.0f),		   glm::vec3(0.0f, 0.0f, -1.f),
			//Triangle left
			glm::vec3(  0.f,  0.5f,   0.f),   glm::vec3(1.0f, 1.0f, 0.0f),      glm::vec2(0.5f, 1.0f),         glm::vec3(1.0f, 0.0f, 0.0f),
			glm::vec3( 0.5f, -0.5f,  0.5f),   glm::vec3(0.0f, 1.0f, 1.0f),      glm::vec2(0.0f, 0.0f),		   glm::vec3(1.0f, 0.0f, 0.0f),
			glm::vec3( 0.5f, -0.5f, -0.5f),   glm::vec3(0.0f, 0.0f, 1.0f),      glm::vec2(1.0f, 0.0f),		   glm::vec3(1.0f, 0.0f, 0.0f)
		};
		unsigned nrOfVertices = sizeof(vertices) / sizeof(Vertex);

		this->set(vertices, nrOfVertices, nullptr, 0);
	}
};