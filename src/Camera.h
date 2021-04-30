#pragma once

#include <iostream>


#include <GL\glew.h>
#include <GLFW\glfw3.h>

#include<glm.hpp>
#include<vec3.hpp>
#include<mat4x4.hpp>
#include<gtc/matrix_transform.hpp>

enum DIRECTION { FORWARD = 0, BACKWARD, LEFT, RIGHT, UP, DOWN };
enum PROJECTION_TYPE { ORTOGRAPHIC, PERSPECTIVE };
enum STATIC_MODE { MOVEFREEDOM, XYFREEDOM,ON };

class Camera {
private:
	glm::mat4 ViewMatrix;
	glm::mat4 ProjectionMatrix;

	GLfloat movementspeed;
	GLfloat sensetivity;
	GLfloat wheel_sensetivity;

	glm::vec3 worldUp;

	glm::vec3 position;
	glm::vec3 front;
	glm::vec3 right;
	glm::vec3 up;



	GLfloat pitch;
	GLfloat yaw;
	GLfloat roll;

	float fov;
	float nearPlane;
	float farPlane;
	int width;
	int height;

	glm::vec2 orthoStart;
	glm::vec2 orthoSize;

	PROJECTION_TYPE ProjectionType = PERSPECTIVE;
	STATIC_MODE staticMode = MOVEFREEDOM;

	void updateCameraVectors() {
			if (staticMode == XYFREEDOM) {
				this->front.x = cos(glm::radians(this->yaw))*cos(glm::radians(this->pitch));
				this->front.y = sin(glm::radians(this->pitch));
				this->front.z = sin(glm::radians(this->yaw))*cos(glm::radians(this->pitch));
				this->front = glm::normalize(this->front);
				this->right = glm::normalize(glm::cross(this->front, this->worldUp));
				this->up = glm::normalize(glm::cross(this->right, this->front));
			}
			else if(staticMode == MOVEFREEDOM) {
				this->front.x = cos(glm::radians(this->yaw))*cos(glm::radians(this->pitch));
				this->front.y = sin(glm::radians(this->pitch));
				this->front.z = sin(glm::radians(this->yaw))*cos(glm::radians(this->pitch));
				this->front = glm::normalize(this->front);


				this->up.x = sin(glm::radians(this->roll));
				this->up.y = cos(glm::radians(this->roll))*cos(glm::radians(this->pitch));
				this->up.z = cos(glm::radians(this->roll))*sin(glm::radians(this->pitch));
				this->up = glm::normalize(this->up);

				this->right.x = cos(glm::radians(this->yaw))*cos(glm::radians(this->roll));
				this->right.y = cos(glm::radians(this->yaw))*sin(glm::radians(this->roll));
				this->right.z = sin(glm::radians(this->yaw));
				this->right = glm::normalize(this->right);
				this->right = glm::normalize(glm::cross(this->front, this->up));
				//Возможно не надо, а если надо, то придется подумать
				//this->worldUp = this->up; 
			}
			else {}
		//std::cout << "f:" << front.x << " " << front.y << " " << front.z << "\n";
		//std::cout << "r:" << right.x << " " << right.y << " " << right.z << "\n";
		//std::cout << "u:" << up.x << " " << up.y << " " << up.z << "\n";
		//std::cout << "\n";
	}

	void refreshProjections(PROJECTION_TYPE type = PERSPECTIVE) {
		this->ProjectionType = type;
		if (this->ProjectionType == PERSPECTIVE)
			setupPerspectiveProjection(this->fov, static_cast<float>(this->width) / this->height, this->nearPlane, this->farPlane);
		else
			setupOrthgraphicProjection(this->orthoStart, this->orthoSize, this->position.z - this->farPlane, this->position.z + this->farPlane);
		//std::cout << "Start: " << this->orthoStart.x << " " << this->orthoStart.x << "\n";
		//std::cout << "Size: " << this->orthoSize.x << " " << this->orthoSize.x << "\n";
	}

	inline void setupPerspectiveProjection(float fov, float aspectRatio, float zNear, float zFar) {
		ProjectionMatrix = glm::perspective(fov, aspectRatio, zNear, zFar);
	}
	inline void setupOrthgraphicProjection(glm::vec2 projectionStart, glm::vec2 projectionSize, float zNear, float zFar) {
		ProjectionMatrix = glm::ortho(projectionStart.x, projectionStart.x + projectionSize.x, projectionStart.y, projectionStart.y + projectionSize.y,zNear, zFar);
	}


public:
	Camera(glm::vec3 position, glm::vec3 direction, glm::vec3 worldUp) {
		this->ViewMatrix = glm::mat4(1.f);
		this->ProjectionMatrix = glm::mat4(1.f);
		this->movementspeed = 1.f;
		this->sensetivity = 5.f;
		this->wheel_sensetivity = 150.f;

		this->worldUp = worldUp;

		this->position = position;
		this->right = glm::vec3(0.f);
		this->up = worldUp;

		this->pitch = 0.f;
		this->yaw = -90.f;
		this->roll = 0.f;

		this->fov = 90.0;
		this->nearPlane = 0.01f;
		this->farPlane = 1000.f;

		this->width = 1920;
		this->height = 1080;

		this->orthoStart = glm::vec2(-this->width/2,-this->height/2);
		this->orthoSize = glm::vec2(this->width, this->height);

		this->updateCameraVectors();
	}
	~Camera() {

	}

	//Accessors
	const glm::mat4 getViewMatrix() {
		this->updateCameraVectors();

		this->ViewMatrix = glm::lookAt(this->position, this->position + this->front, this->up);

		return this->ViewMatrix;
	}

	const glm::mat4 getProjectionMatrix() {
		this->refreshProjections(this->ProjectionType);
		return this->ProjectionMatrix;
	}



	const glm::vec3 getPosition() const { return this->position; }
	//Functioms

	void move_p(const float &dt, const int direction) {

		if (staticMode == MOVEFREEDOM || staticMode == XYFREEDOM) {
			//Update position vector
			switch (direction)
			{
			case FORWARD:
				this->position += this->front * this->movementspeed * dt;
				break;
			case BACKWARD:
				this->position -= this->front * this->movementspeed * dt;
				break;
			case LEFT:
				this->position -= this->right * this->movementspeed * dt;
				break;
			case RIGHT:
				this->position += this->right * this->movementspeed * dt;
				break;
			case UP:
				this->position -= this->up * this->movementspeed * dt;
				break;
			case DOWN:
				this->position += this->up * this->movementspeed * dt;
				break;
			default:
				break;
			}
		}
		else {}
	}
	void move_o(const float &dt, const int direction) {
		if (staticMode == MOVEFREEDOM || staticMode == XYFREEDOM) {
			switch (direction)
			{
			case LEFT:
				this->orthoStart.x -= this->movementspeed * 1 * dt;
				break;
			case RIGHT:
				this->orthoStart.x += this->movementspeed * 1 * dt;
				break;
			case UP:
				this->orthoStart.y += this->movementspeed * 1 * dt;
				break;
			case DOWN:
				this->orthoStart.y -= this->movementspeed * 1 * dt;
				break;
			default:
				break;
			}
		}
		else {}
	}



	void updateMouseInput(const float &dt, const double offsetX, const double offsetY, const double wheelOffset) {
			//Update pitch, yam and roll 
			if (staticMode == XYFREEDOM) {
				this->pitch += static_cast<GLfloat>(offsetY) * this->sensetivity *dt;
				this->yaw += static_cast<GLfloat>(offsetX) * this->sensetivity *dt;
				if (this->ProjectionType == PERSPECTIVE)
					this->roll += static_cast<GLfloat>(wheelOffset) * this->wheel_sensetivity *dt;
				else {
					float ratio = static_cast<float>(this->width) / this->height;
					this->orthoSize.x += static_cast<int>(this->wheel_sensetivity*wheelOffset *ratio * dt);
					this->orthoSize.y += static_cast<int>(this->wheel_sensetivity*wheelOffset * dt);
					this->orthoStart.x -= static_cast<int>(this->wheel_sensetivity*wheelOffset *ratio / 2 * dt);
					this->orthoStart.y -= static_cast<int>(this->wheel_sensetivity*wheelOffset / 2 * dt);
				}
			}
			else if(staticMode == MOVEFREEDOM){

				this->pitch += (static_cast<GLfloat>(offsetY)*cos(glm::radians(this->roll)) - static_cast<GLfloat>(offsetX)*sin(glm::radians(this->roll))) * this->sensetivity *dt;
				this->yaw += (static_cast<GLfloat>(offsetX)*cos(glm::radians(this->roll)) + static_cast<GLfloat>(offsetY)*sin(glm::radians(this->roll))) * this->sensetivity *dt;
				if (this->ProjectionType == PERSPECTIVE)
					this->roll += static_cast<GLfloat>(wheelOffset) * this->wheel_sensetivity *dt;
				else {
					float ratio = static_cast<float>(this->width) / this->height;
					this->orthoSize.x += static_cast<int>(this->wheel_sensetivity*wheelOffset *ratio  * dt);
					this->orthoSize.y += static_cast<int>(this->wheel_sensetivity*wheelOffset * dt);
					this->orthoStart.x -= static_cast<int>(this->wheel_sensetivity*wheelOffset *ratio / 2 * dt);
					this->orthoStart.y -= static_cast<int>(this->wheel_sensetivity*wheelOffset / 2 * dt);


				}
			}
			else {}

			if (this->pitch > 90.f)
				this->pitch = 90.f;
			if (this->pitch < -90.f)
				this->pitch = -90.f;

			if (this->roll > 360.f || this->roll < -360.f)
				this->roll = 0.f;
			if (this->pitch > 360.f || this->pitch < -360.f)
				this->pitch = 0.f;
		
	}

	void updateInput(const float &dt, const int direction, const double& offsetX, const double offsetY, const double wheelOffset) {
		this->updateMouseInput(dt, offsetX, offsetY, wheelOffset);
	}
	PROJECTION_TYPE getProjectionType() {
		return this->ProjectionType;
	}
	STATIC_MODE getStaticMode() {
		return this->staticMode;
	}

	void changeProjType() {
		switch (ProjectionType) {
		case PERSPECTIVE:
			ProjectionType = ORTOGRAPHIC;
			break;
		case ORTOGRAPHIC:
			ProjectionType = PERSPECTIVE;
			break;
		default:
			break;
		}
	}
	void changeStaticMode() {
		switch (staticMode) {
		case MOVEFREEDOM:
			staticMode = XYFREEDOM;
			break;
		case XYFREEDOM:
			staticMode = ON;
			break;
		case ON:
			staticMode = MOVEFREEDOM;
			break;
		default:
			break;
		}
	}

	void setToDefault() {
		this->pitch = 0.f;
		this->yaw = -90.f;
		this->roll = 0.f;
		this->right = glm::vec3(0.f);
		this->up = this->worldUp;
		this->updateCameraVectors();
	}

};



