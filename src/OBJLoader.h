#pragma once

#include<iostream>
#include<string>
#include<fstream>
#include<vector>
#include<sstream>
#include<algorithm>

#include <GL\glew.h>
#include <GLFW\glfw3.h>

#include<glm.hpp>
#include<vec3.hpp>
#include<vec4.hpp>
#include<mat4x4.hpp>
#include<gtc/matrix_transform.hpp>
#include<gtc/type_ptr.hpp>

#include"Vertex.h"


static std::vector<Vertex> loadOBJ(const std::string& filePath) {
	//Vertex portion
	std::vector<glm::fvec3> vertex_positions;
	std::vector<glm::fvec2> vertex_texcoords;
	std::vector<glm::fvec3> vertex_normals;

	//Face vectors
	std::vector<GLint> vertex_position_indices;
	std::vector<GLint> vertex_texcoord_indices;
	std::vector<GLint> vertex_normal_indices;

	//Vertex array
	std::vector<Vertex> vertices;

	std::stringstream ss;
	std::ifstream in_file(filePath);
	std::string line = "";
	std::string prefix = "";
	glm::vec3 temp_vec3;
	glm::vec2 temp_vec2;
	GLint temp_glint = 0;

	//File open error check
	if (!in_file.is_open()) {
		throw "ERROR::OBJLOADER::Could not open file";
	}
	//Read one line at a time
	while (std::getline(in_file, line)) {
		//Get the prefix of the line
		ss.clear();
		ss.str(line);
		ss >> prefix;
		if (prefix == "#") {}
		else if (prefix == "o") {}
		else if (prefix == "s") {}
		else if (prefix == "use_mtl") {

		}
		else if (prefix == "v") { //Vertex position
			ss >> temp_vec3.x >> temp_vec3.y >> temp_vec3.z;
			vertex_positions.push_back(temp_vec3);
		}
		else if (prefix == "vt") { //Vertex texcoord
			ss >> temp_vec2.x >> temp_vec2.y;
			vertex_texcoords.push_back(temp_vec2);
		}
		else if (prefix == "vn") { //Vertex normals
			ss >> temp_vec3.x >> temp_vec3.y >> temp_vec3.z;
			vertex_normals.push_back(temp_vec3);
		}
		else if (prefix == "f") { //faces
			int counter = 0;
			while (ss >> temp_glint) {
				//Pushing indices into correct arrays
				if (counter == 0)
					vertex_position_indices.push_back(temp_glint);
				else if (counter == 1)
					vertex_texcoord_indices.push_back(temp_glint);
				else if (counter == 2)
					vertex_normal_indices.push_back(temp_glint);
				//Handling characters
				if (ss.peek() == '/') {
					++counter;
					ss.ignore(1, '/');
				}
				else if (ss.peek() == ' ') {
					++counter;
					ss.ignore(1, ' ');
				}

				if (counter > 2)
					counter = 0;
	
			}
		}
		else { }

		//Build final vertex array (mesh)
		vertices.resize(vertex_position_indices.size(), Vertex());
		//Load in all vertices
		for (size_t i = 0; i < vertices.size(); i++) {
			//std::cout << vertex_positions.size()<< vertex_position_indices[i]-1;
			vertices[i].position = vertex_positions[vertex_position_indices[i]-1];
			vertices[i].texcoord = vertex_texcoords[vertex_texcoord_indices[i]-1];
			vertices[i].normal = vertex_normals[vertex_normal_indices[i]-1];
			vertices[i].color = glm::vec3(1.f, 1.f, 1.f);
		}

		//DEBUG
		//std::cout << "Nr of vertices: " << vertices.size() << '\n';
	}

	//Loaded success
	std::cout << "OBJ file loaded. Nr of vertices: "<< vertices.size() << '\n';
	return vertices;
}