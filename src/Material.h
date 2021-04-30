#pragma once


#include <GL\glew.h>
#include <GLFW\glfw3.h>

#include<glm.hpp>
#include<vec2.hpp>
#include<vec3.hpp>
#include<vec4.hpp>
#include<mat4x4.hpp>
#include<gtc/matrix_transform.hpp>
#include<gtc/type_ptr.hpp>

#include "Shader.h"

class Material {
private:
	glm::vec3 ambient;
	glm::vec3 diffuse;
	glm::vec3 specular;
	GLint diffuseTex;
	GLint specularTex;
public:
	Material() {}
	Material(glm::vec3 p_ambient, glm::vec3 p_diffuse, glm::vec3 p_specular, GLint p_diffuseTex, GLint p_specularTex)
		:ambient(p_ambient), diffuse(p_diffuse), specular(p_specular), diffuseTex(p_diffuseTex), specularTex(p_specularTex) { }
	Material(Material* material){
		ambient = material->ambient;
		diffuse = material->diffuse;
		specular = material->specular;
		diffuseTex = material->diffuseTex;
		specularTex = material->specularTex;
	}
	~Material(){}

	void ChangeLighting(glm::vec3 p_ambient, glm::vec3 p_diffuse, glm::vec3 p_specular) {
		ambient = p_ambient;
		diffuse = p_diffuse;
		specular = p_specular;
	}
	//Functions
	void sendToShader(Shader& program) {
		program.setVec3f(this->ambient, "material.ambient");
		program.setVec3f(this->diffuse, "material.diffuse");
		program.setVec3f(this->specular, "material.specular");
		program.set1i(this->diffuseTex, "material.diffuseTex");
		program.set1i(this->specularTex, "material.specularTex");
	}
};
