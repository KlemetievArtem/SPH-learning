#pragma once

#include"Mesh.h"
#include"Texture.h"
#include"Shader.h"
#include"Material.h"
#include"OBJLoader.h"

class Model {
private:
	Material material;
	Texture* overrideTextureDiffuse;
	Texture* overrideTextureSpecular;
public:
	std::vector<Mesh*> meshes;
private:
	glm::vec3 position;
	void updateUniforms() {

	}
 
public:
	Model(glm::vec3 position, Material material, Texture* orTexDif,
		Texture* orTexSpec, std::vector<Mesh*> meshes) {
		this->position = position;
		this->material = material;
		this->overrideTextureDiffuse = orTexDif;
		this->overrideTextureSpecular = orTexSpec;

		for (auto*i : meshes) {
			this->meshes.push_back(new Mesh(*i));
		}

		for (auto&i : this->meshes) {
			i->move(this->position);
			i->setOrigin(this->position);
		}
	}
	Model(glm::vec3 position, Material material, Texture* orTexDif, Texture* orTexSpec, Mesh* mesh) {
		this->position = position;
		this->material = material;
		this->overrideTextureDiffuse = orTexDif;
		this->overrideTextureSpecular = orTexSpec;

		this->meshes.push_back(new Mesh(mesh));
	

		for (auto&i : this->meshes) {
			i->move(this->position);
			i->setOrigin(this->position);
		}
	}
	//OBJ file loaded model
	Model(glm::vec3 position, Material material, Texture* orTexDif, Texture* orTexSpec, const std::string& filePath) {
		this->position = position;
		this->material = material;
		this->overrideTextureDiffuse = orTexDif;
		this->overrideTextureSpecular = orTexSpec;

		std::vector<Vertex> mesh = loadOBJ(filePath);
		meshes.push_back(new Mesh(mesh.data(), mesh.size(), NULL, 0, glm::vec3(1.f, 0.f, 0.f), glm::vec3(0.f), glm::vec3(0.f), glm::vec3(1.f)));

		for (auto&i : this->meshes) {
			i->move(this->position);
			i->setOrigin(this->position);
		}
	}

	~Model() {
		for (auto*&i : this->meshes)
			delete i;
	}

	//Functions
	void rotate(glm::vec3 rotation) {
		for (auto& i : this->meshes)
			i->rotate(rotation);
	}

	void move(glm::vec3 vec) {
		for (auto& i : this->meshes)
			i->move(vec);
	}
	void scaleUp(glm::vec3 vec) {
		for (auto& i : this->meshes)
			i->scaleUp(vec);
	}

	void moveTo(glm::vec3 vec) {
		for (auto& i : this->meshes)
			i->moveTo(vec);
	}
	void scaleUpTo(glm::vec3 vec) {
		for (auto& i : this->meshes)
			i->changeScaleTo(vec);
	}
	

	void update() {

	}

	void render(Shader* shader) {
		this->updateUniforms();

		//Use a program
		shader->use();
		//Update the uniforms
		this->updateUniforms();
		//Update uniforms											
		this->material.sendToShader(*shader);
		//Activate texture
		//Draw
		for (auto& i : this->meshes) {
			this->overrideTextureDiffuse->bind(0);
			this->overrideTextureSpecular->bind(1);
			i->render(shader);
		}
	}

	glm::vec3 getPosition() { return this->position; }
	Material* getMaterial() { return &material; }

};
