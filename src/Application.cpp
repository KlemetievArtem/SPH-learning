#include "Application.h"


#define BIND_EVENT_FN(x) std::bind(&Application::x, this, std::placeholders::_1)




//General functions

double Application::wheeloffset = 0;


//Private functions

void Application::initGLFW()
{
	//INIT GLFW
	if (glfwInit() == GLFW_FALSE) {
		std::cout << "ERROR::GLFW_INIT_FAILED\n";
		glfwTerminate();
	}
}

void Application::initWindow(const char*title,bool resizeable){
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, this->GL_VERSION_MAJOR);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, this->GL_VERSION_MINOR);
	glfwWindowHint(GLFW_RESIZABLE, resizeable);

	this->window = glfwCreateWindow(this->WINDOW_WIDTH, this->WINDOW_HEIGHT, title, NULL, NULL);
	glfwSetFramebufferSizeCallback(this->window, Application::framebuffer_resize_callback);

	if (this->window == nullptr) {
		std::cout << "ERROR::GLFW_WINDOW_INIT_FAILED\n";
		glfwTerminate();
	}

	glfwSetWindowUserPointer(window, &m_Data);
	this->SetEventCallback(BIND_EVENT_FN(EmptyEventFunction));


	glfwGetFramebufferSize(this->window, &this->framebufferWidth, &this->framebufferHeight);
	//glfwSetScrollCallback(this->window, Application::Mouse_scrool_callback);

	glfwSetScrollCallback(window, [](GLFWwindow* window, double xOffset, double yOffset) {
		WindowData& data = *(WindowData*)glfwGetWindowUserPointer(window);

		MouseScrolledEvent event((float)xOffset, (float)yOffset);
		data.EventCallback(event);
	});

	glfwSetCursorPosCallback(window, [](GLFWwindow* window, double xPos, double yPos)
	{
		WindowData& data = *(WindowData*)glfwGetWindowUserPointer(window);

		MouseMovedEvent event((float)xPos, (float)yPos);
		data.EventCallback(event);
	});

	glfwSetKeyCallback(window, [](GLFWwindow* window, int key, int scancode, int action, int mods)
	{
		WindowData& data = *(WindowData*)glfwGetWindowUserPointer(window);

		switch (action)
		{
		case GLFW_PRESS:
		{
			KeyPressedEvent event(key, 0);
			data.EventCallback(event);
			break;
		}
		case GLFW_RELEASE:
		{
			KeyReleasedEvent event(key);
			data.EventCallback(event);
			break;
		}
		case GLFW_REPEAT:
		{
			KeyPressedEvent event(key, 1);
			data.EventCallback(event);
			break;
		}
		}
	});

	glfwSetMouseButtonCallback(window, [](GLFWwindow* window, int button, int action, int mods)
	{
		WindowData& data = *(WindowData*)glfwGetWindowUserPointer(window);

		switch (action)
		{
		case GLFW_PRESS:
		{
			MouseButtonPressedEvent event(button);
			data.EventCallback(event);
			break;
		}
		case GLFW_RELEASE:
		{
			MouseButtonReleasedEvent event(button);
			data.EventCallback(event);
			break;
		}
		}
	});

/*
	glfwSetScrollCallback(this->window, [](GLFWwindow* window, double xOffset, double yOffset) {
		WindowData& data = *(WindowData*)glfwGetWindowUserPointer(window);
		
		MouseScrolledEvent event((float)xOffset, (float)yOffset);
		std::cout << event.ToString();
		std::cout << "\n&data" << &data << "\n&(data.EventCallback)" << &(data.EventCallback) << "\n&event"<< &event << "\n&window" << &window << "\nglfwGetWindowUserPointer(window)" << glfwGetWindowUserPointer(window);
		//data.EventCallback(event);
		//data.EventCallback(event);
	});
*/

	glViewport(0, 0, framebufferWidth, framebufferHeight);

	glfwMakeContextCurrent(this->window);
}

void Application::initGLEW() {
	glewExperimental = GL_TRUE;

	if (glewInit() != GLEW_OK) {
		std::cout << "ERROR:MAIN.CPP::GLEW_INIT_FAILED\n";
		glfwTerminate();
	}

}

void Application::initOpenGLOptions() {//OPENGL OPTIONS
	glEnable(GL_DEPTH_TEST);

	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
	glFrontFace(GL_CCW);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	//Input
	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
}

void Application::initMatrices() {

	this->ViewMatrix = glm::mat4(1.0f);
	this->ViewMatrix = glm::lookAt(this->camPosition, this->camPosition + this->camFront, this->worldUp);

	this->ProjectionMatrix = glm::mat4(1.f);
	this->ProjectionMatrix = glm::perspective(glm::radians(this->fov), static_cast<float>(this->framebufferWidth) / this->framebufferHeight, this->nearPlane, this->farPlane);
}

void Application::initShaders() {
	this->shaders.push_back(new Shader(this->GL_VERSION_MAJOR, this->GL_VERSION_MINOR, "res/vertex_core.glsl", "res/fragment_core.glsl"));
	//this->shaders.push_back(new Shader(this->GL_VERSION_MAJOR, this->GL_VERSION_MINOR, "res/vertex_core.glsl", "res/fragment_core.glsl"));
}

void Application::initTextures() {
	//TEXTURE 0
	this->textures.push_back(new Texture("res/textures/GarfildThiccSexy.png", GL_TEXTURE_2D));
	this->textures.push_back(new Texture("res/textures/GarfildThiccSexy_specular.png", GL_TEXTURE_2D));
	//TEXTURE 1
	this->textures.push_back(new Texture("res/textures/GarfildThiccSexyUpsideDown.PNG", GL_TEXTURE_2D));
	this->textures.push_back(new Texture("res/textures/GarfildThiccSexyUpsideDown_specular.PNG", GL_TEXTURE_2D));
}

void Application::initMaterials() {
	//MATERIAL0
	this->materials.push_back(new Material(glm::vec3(0.3f), glm::vec3(0.9f), glm::vec3(1.f), 0, 1));

}

void Application::initOBJModels() {
	//std::vector<Vertex> temp;
	//temp = loadOBJ("res/OBJFiles/figny.obj");

}

void Application::initCursorMesh() {
	//this->EssentialMeshes.push_back(new Quad())
	this->EssentialMeshes.push_back(new Mesh(&Triangle(glm::vec3(0.f), glm::vec3(0.1f, 0.0f, 0.0f), glm::vec3(0.0f, 0.1f, 0.0f), glm::vec3(1.0f, 0.0f, 0.0f)), glm::vec3(0.f, 0.f, 0.f), glm::vec3(0.f), glm::vec3(0.f), glm::vec3(0.1f)));
	this->EssentialMeshes.push_back(new Mesh(&Triangle(glm::vec3(0.f), glm::vec3(0.0f, 0.1f, 0.0f), glm::vec3(0.0f, 0.0f, 0.1f), glm::vec3(0.0f, 1.0f, 0.0f)), glm::vec3(0.f, 0.f, 0.f), glm::vec3(0.f), glm::vec3(0.f), glm::vec3(0.1f)));
	this->EssentialMeshes.push_back(new Mesh(&Triangle(glm::vec3(0.f), glm::vec3(0.0f, 0.0f, 0.1f), glm::vec3(0.1f, 0.0f, 0.0f), glm::vec3(0.0f, 0.0f, 1.0f)), glm::vec3(0.f, 0.f, 0.f), glm::vec3(0.f), glm::vec3(0.f), glm::vec3(0.1f)));
}

void Application::initModels() {


	this->models.push_back(new Model(glm::vec3(0.f), this->materials[0], this->textures[TEX_GTSud], this->textures[TEX_GTSud_SPECULAR], EssentialMeshes));
	std::vector<Mesh*> meshes;
	//meshes.push_back(new Mesh(&Pyramid(), glm::vec3(1.f, 0.f, 0.f), glm::vec3(0.f), glm::vec3(0.f), glm::vec3(1.f)));
	//meshes.push_back(new Mesh(&Quad(), glm::vec3(0.f, 0.f, 0.f), glm::vec3(0.f), glm::vec3(0.f), glm::vec3(1.f)));
	//meshes.push_back(new Mesh(&Something(), glm::vec3(0.f, 0.f, 0.f), glm::vec3(0.f), glm::vec3(0.f), glm::vec3(1.f)));

	//this->meshes.push_back(new Mesh(&Quad(), glm::vec3(0.f), glm::vec3(0.f), glm::vec3(1.f)));

	//this->models.push_back(new Model(glm::vec3(0.f), this->materials[0], this->textures[TEX_GTSud], this->textures[TEX_GTSud_SPECULAR], meshes));
	//this->models.push_back(new Model(glm::vec3(0.f,1.f,1.f), this->materials[0], this->textures[TEX_GTS], this->textures[TEX_GTS_SPECULAR], meshes));
	//this->models.push_back(new Model(glm::vec3(2.f, 0.f, 2.f), this->materials[0], this->textures[TEX_GTSud], this->textures[TEX_GTSud_SPECULAR], meshes));
	//this->models.push_back(new Model(glm::vec3(4.f, 0.f, 4.f), this->materials[0], this->textures[TEX_GTSud], this->textures[TEX_GTSud_SPECULAR], "res/OBJFiles/10053_Walrus_v1_L3.obj"));

	this->models.push_back(new Model(glm::vec3(0.f), this->materials[0], this->textures[TEX_GTSud], this->textures[TEX_GTSud_SPECULAR], meshes));


	for (auto*& i : meshes)
		delete i;
}

void Application::initLights() {
	//LIGHTS
	//this->ligts.push_back(new glm::vec3(sin(glm::radians(light_fi))*light_radius, light_radius, cos(glm::radians(light_fi))*light_radius));
	this->ligts.push_back(new glm::vec3(2.f,2.f,2.f));
}

void Application::initUniforms()
{
	
	this->shaders[SHADER_CORE_PROGRAM]->use();         
	this->shaders[SHADER_CORE_PROGRAM]->setMat4fv(ViewMatrix, "ViewMatrix");                 
	this->shaders[SHADER_CORE_PROGRAM]->setMat4fv(ProjectionMatrix, "ProjectionMatrix");   
	this->shaders[SHADER_CORE_PROGRAM]->setVec3f(*this->ligts[0], "lightPos0");      


	this->shaders[SHADER_CORE_PROGRAM]->unuse();
}

void Application::updateUniforms() {
    //Update ViewMatrix (camera)      
	this->ViewMatrix = this->camera.getViewMatrix();
	this->shaders[SHADER_CORE_PROGRAM]->setMat4fv(this->ViewMatrix, "ViewMatrix");
	this->shaders[SHADER_CORE_PROGRAM]->setVec3f(this->camera.getPosition(), "cameraPos");

	//Update framebuffer size and projection matrix
	glfwGetFramebufferSize(this->window, &this->framebufferWidth, &this->framebufferHeight);
	
	this->ProjectionMatrix = this->camera.getProjectionMatrix();

	//this->ProjectionMatrix = glm::perspective(glm::radians(this->fov), static_cast<float>(this->framebufferWidth) / this->framebufferHeight, this->nearPlane, this->farPlane);
	this->shaders[SHADER_CORE_PROGRAM]->setMat4fv(this->ProjectionMatrix, "ProjectionMatrix");

	//light_fi += 1.f;
	//*this->ligts[0] = glm::vec3(sin(glm::radians(light_fi))*light_radius, light_radius, cos(glm::radians(light_fi))*light_radius);
	//this->shaders[SHADER_CORE_PROGRAM]->setVec3f(glm::vec3(*this->ligts[0]), "lightPos0");



}

//Constructord and distructors
Application::Application(const char* title,
	const int WINDOW_WIDTH, const int WINDOW_HEIGHT,
	const int GL_VERSION_MAJOR, const int GL_VERSION_MINOR, bool resizable)
	: WINDOW_WIDTH(WINDOW_WIDTH),WINDOW_HEIGHT(WINDOW_HEIGHT), GL_VERSION_MAJOR(GL_VERSION_MAJOR), GL_VERSION_MINOR(GL_VERSION_MINOR),
	camera(glm::vec3(0.5f, 0.5f, 0.5f), glm::vec3(0.f, 0.f, -1.f), glm::vec3(0.f, 1.f, 0.f))
{
	//Init variables
	this->window = nullptr;
	this->framebufferWidth = this->WINDOW_WIDTH;
	this->framebufferHeight = this->WINDOW_HEIGHT;


	this->camPosition = glm::vec3(0.f, 0.f, 2.f);
	this->worldUp = glm::vec3(0.f, 1.f, 0.f);
	this->camFront = glm::vec3(0.f, 0.f, -1.f);

	this->fov = 90.0f;
	this->nearPlane = 0.1f;
	this->farPlane = 1000.f;

	this->dt = 0.f;
	this->curTime = 0.f;
	this->lastTime = 0.f;

	this->lastMouseX = 0.f;
	this->lastMouseY = 0.f;
	this->mouseX = 0.f;
	this->mouseY = 0.f;
	this->mouseOffsetX = 0.f;
	this->mouseOffsetY = 0.f;
	this->firstMouse = true;



	this->cursorPosition = glm::vec3(0.f);





	this->initGLFW();
	this->initWindow(title, resizable);
	this->initGLEW();
	this->initOpenGLOptions();
	this->initMatrices();
	this->initShaders();
	this->initTextures();
	this->initMaterials();
	this->initOBJModels();
	this->initCursorMesh();
	this->initModels();
	this->initLights();
	this->initUniforms();



	this->initIMGUI();
	this->optionsRegistration();


	this->CompDomainInit();



	if(CalculateParallelToRender){
		this->initThread();
	}
}
void Application::initThread() {
	//mainCalculationThread = std::make_unique<std::thread>([]()->void {
	//	std::cout << "Entering thread: " << std::this_thread::get_id() << "\n";
	//});
	mainCalculationThread = std::make_unique<std::thread>([this]() -> void {
		while (this->MCThreadLoopCondition) {
			for (auto*&i : this->ComputationalDomains) {
				//std::cout << "first\n";
				i->timeStep_thread(i->getDeltaTime(),std::ref(this->dataReadyForRender), std::ref(this->dataIsRendering));
			}
		}
	});	
}
void Application::calculationLoop() {
	std::cout << "Entering thread: " << std::this_thread::get_id() << "\n";
}

void Application::joinThread() {
	MCThreadLoopCondition = false;
	mainCalculationThread->join();
}


Application::~Application()
{

	if (CalculateParallelToRender) {
		this->joinThread();
	}




	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();

	glfwDestroyWindow(this->window);
	glfwTerminate();

	for (size_t i = 0;i < this->shaders.size();i++)
		safeDelete(&this->shaders[i]);
		//delete this->shaders[i];
	for (size_t i = 0;i < this->textures.size();i++)
		safeDelete(&this->textures[i]);
		//delete this->textures[i];
	for (size_t i = 0;i < this->materials.size();i++)
		safeDelete(&this->materials[i]);
		//delete this->materials[i];
	//for (size_t i = 0;i < this->meshes.size();i++)
	//	delete this->meshes[i];
	for (auto*&i : this->models)
		safeDelete(&i);
		//delete i;
	for (auto*&i : this->RenewableModels)
		safeDelete(&i);
		//delete i;
	for (size_t i = 0;i < this->ligts.size();i++)
		safeDelete(&this->ligts[i]);
		//delete this->ligts[i];
	//for (size_t i = 0;i < this->eventQueue.size();i++)
	//	delete this->eventQueue[i];


}
//Accessor
int Application::getWindowShouldClouse()
{
	return glfwWindowShouldClose(this->window);
}
//Modifires
void Application::setWindowShouldClose()
{
	glfwSetWindowShouldClose(this->window, GLFW_TRUE);
}

//Functions

void Application::updateDt() {
	this->curTime = static_cast<float>(glfwGetTime());
	this->dt = this->curTime - this->lastTime;
	this->lastTime = this->curTime;

}

void Application::updateMouseInput() {
	glfwGetCursorPos(this->window, &this->mouseX, &this->mouseY);

	if (this->firstMouse) {
		this->lastMouseX = this->mouseX;
		this->lastMouseY = this->mouseY;
		this->lastWheel = 0;
		this->firstMouse = false;
	}
	//Calc offset
	this->mouseOffsetX = this->mouseX - this->lastMouseX;
	this->mouseOffsetY = this->lastMouseY - this->mouseY;

	this->lastWheel = this->lastWheel - this->wheeloffset;
	//Set last X ad Y
	this->lastMouseX = this->mouseX;
	this->lastMouseY = this->mouseY;

	this->wheeloffset = this->lastWheel;
}

//MACROS
/*
#define RESIZE_IF_NEEDED(array,key) if (array.size()<key) array.resize(key+1);

#define KeyInputPress(key,body) if ((glfwGetKey(this->window, key) == GLFW_PRESS)) {\
									RESIZE_IF_NEEDED(s_Presscount, key)\
									body\
									++s_Presscount[key];\
								}
#define KeyInputRelease(key,body) if ((glfwGetKey(this->window, key) == GLFW_RELEASE)) {\
									RESIZE_IF_NEEDED(s_Presscount, key)\
									body\
									s_Presscount[key] = 0;\
								}
*/
void Application::updateKeyboardInput() {

	/*
	//Program
	KeyInputPress(GLFW_KEY_ESCAPE, glfwSetWindowShouldClose(this->window, GLFW_TRUE);)
	KeyInputRelease(GLFW_KEY_ESCAPE, {})
	//Save point to constructor
	KeyInputPress(GLFW_KEY_ENTER, if (s_Presscount[GLFW_KEY_ENTER] > 20 || s_Presscount[GLFW_KEY_ENTER] == 0) {
		this->toolConstructor.remember(cursorPosition.x, cursorPosition.y, cursorPosition.z);
		if (this->toolConstructor.checkSize() % 3 == 0) {

			std::vector<Mesh*> meshes;

			this->toolConstructor.addTriangleTo(&meshes, &EssentialMeshes);

			this->models.push_back(new Model(glm::vec3(0.f), this->materials[0], this->textures[TEX_GTSud], this->textures[TEX_GTSud_SPECULAR], meshes));

			for (auto*& i : meshes)
				delete i;
		}
	})
	KeyInputRelease(GLFW_KEY_ENTER, {})
	//Save primitive to a file
	KeyInputPress(GLFW_KEY_P, if (s_Presscount[GLFW_KEY_P] > 4 || s_Presscount[GLFW_KEY_P] == 0) {
			this->toolConstructor.savePrimitive("src/Something.h", "Something");
		})
	KeyInputRelease(GLFW_KEY_P, {})
	KeyInputPress(GLFW_KEY_C, this->camera.changeProjType();)
	KeyInputRelease(GLFW_KEY_C, {})
	KeyInputPress(GLFW_KEY_SPACE, if (s_Presscount[GLFW_KEY_SPACE] > 4 || s_Presscount[GLFW_KEY_SPACE] == 0) {
			this->camera.changeStaticMode();
			switch (this->camera.getStaticMode()) {
			case ON:
				glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
				break;
			default:
				glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
				break;
			}
		})
	KeyInputRelease(GLFW_KEY_SPACE, {})
	//Camera and Constructor
	/*
	KeyInputPress(GLFW_KEY_W, if (s_Presscount[GLFW_KEY_W] > 8 || s_Presscount[GLFW_KEY_W] == 0) {
			if (this->camera.getProjectionType() == PERSPECTIVE)
				this->camera.move_p(this->dt, FORWARD);
			else
				this->camera.move_o(this->dt, UP);
			if (this->camera.getStaticMode() == ON) {
				//this->ModelMatrix = glm::translate(this->ModelMatrix, glm::vec3(0.f, 1.f, 0.f));
				for (auto &i : this->models)
					i->move(this->dt*glm::vec3(0.f, 1.f, 0.f));
				this->EssentialMeshes[0]->move(this->dt*glm::vec3(0.f, 1.f, 0.f));
			}
		})
	KeyInputRelease(GLFW_KEY_W, {})
	KeyInputPress(GLFW_KEY_S, if (s_Presscount[GLFW_KEY_S] > 8 || s_Presscount[GLFW_KEY_S] == 0) {
			if (this->camera.getProjectionType() == PERSPECTIVE)
				this->camera.move_p(this->dt, BACKWARD);
			else
				this->camera.move_o(this->dt, DOWN);
			if (this->camera.getStaticMode() == ON) {
				//this->ModelMatrix = glm::translate(this->ModelMatrix, glm::vec3(0.f, -1.f, 0.f));
				for (auto &i : this->models)
					i->move(this->dt*glm::vec3(0.f, -1.f, 0.f));
				this->EssentialMeshes[0]->move(this->dt*glm::vec3(0.f, -1.f, 0.f));
			}
		})
	KeyInputRelease(GLFW_KEY_S, {})
	KeyInputPress(GLFW_KEY_A, if (s_Presscount[GLFW_KEY_A] > 8 || s_Presscount[GLFW_KEY_A] == 0) {
			if (this->camera.getProjectionType() == PERSPECTIVE)
				this->camera.move_p(this->dt, LEFT);
			else
				this->camera.move_o(this->dt, LEFT);
			if (this->camera.getStaticMode() == ON) {
				//this->ModelMatrix = glm::translate(this->ModelMatrix, glm::vec3(-1.f, 0.f, 0.f));
				for (auto &i : this->models)
					i->move(this->dt*glm::vec3(-1.f, 0.f, 0.f));
				this->EssentialMeshes[0]->move(this->dt*glm::vec3(-1.f, 0.f, 0.f));
			}
		})
	KeyInputRelease(GLFW_KEY_A, {})
	KeyInputPress(GLFW_KEY_D, if (s_Presscount[GLFW_KEY_D] > 8 || s_Presscount[GLFW_KEY_D] == 0) {
			if (this->camera.getProjectionType() == PERSPECTIVE)
				this->camera.move_p(this->dt, RIGHT);
			else
				this->camera.move_o(this->dt, RIGHT);
			if (this->camera.getStaticMode() == ON) {
				//this->ModelMatrix = glm::translate(this->ModelMatrix, glm::vec3(1.f, 0.f, 0.f));
				for (auto &i : this->models)
					i->move(this->dt*glm::vec3(1.f, 0.f, 0.f));
				this->EssentialMeshes[0]->move(this->dt*glm::vec3(1.f, 0.f, 0.f));
			}
		})
	KeyInputRelease(GLFW_KEY_D, {})
	

	KeyInputPress(GLFW_KEY_R, this->camera.setToDefault();)
	KeyInputRelease(GLFW_KEY_R, {})
	KeyInputPress(GLFW_KEY_Q, this->camera.move_p(this->dt, DOWN);)
	KeyInputRelease(GLFW_KEY_Q, {})
	KeyInputPress(GLFW_KEY_E, this->camera.move_p(this->dt, UP);)
	KeyInputRelease(GLFW_KEY_E, {})

	*/
	//glfwSetScrollCallback(this->window, UserCallBackScroll);
}

void Application::updateInput() {
	glfwPollEvents();


	ImGui_ImplOpenGL3_NewFrame();
	ImGui_ImplGlfw_NewFrame();
	ImGui::NewFrame();

	this->updateKeyboardInput();
	this->updateMouseInput();
	this->camera.updateInput(dt, -1, this->mouseOffsetX, this->mouseOffsetY, this->wheeloffset);


	if (this->camera.getStaticMode() == ON) {
		//this->ModelMatrix = glm::translate(this->ModelMatrix, glm::vec3(1.f, 0.f, 0.f));
		for (auto &i : this->models)
			i->move(this->dt*glm::vec3(wheeloffset));
			this->EssentialMeshes[0]->move(this->dt*glm::vec3(wheeloffset));
	}




	//std::cout << this->wheeloffset << lastWheel<<'\n';
}

void Application::update() {
	//UPDATE COMPUTATIONAL DOMAIN

	//ПАРАЛЛЕЛИМ

	if (CalculateParallelToRender) { }
	else {
		for (auto*&i : this->ComputationalDomains) {
			//std::cout << "first\n";
			i->timeStep(i->getDeltaTime());
		}
	}



	//UPDATE INPUT
	
	//if (this->currentOption->getFlag()) {
	//	this->currentOption->setFlag(false);
	//	this->currentOption->outside();
	//
	//
	//	std::vector<Mesh*> meshes;
	//
	//	if (this->currentOption) {
	//		this->currentOption->InitModels(&meshes, &models);
	//	}
	//	//std::cout << meshes.size();
	//	this->models.push_back(new Model(glm::vec3(0.f), this->materials[0], this->textures[TEX_GTSud], this->textures[TEX_GTSud_SPECULAR], meshes));
	//	for (auto*& i : meshes)
	//		delete i;
	//}
	

	this->CursorUpdate();


	this->updateDt();
	this->updateInput();

	//this->models[0]->rotate(glm::vec3(0.f, 1.f, 0.f));
	//this->models[1]->rotate(glm::vec3(0.f, 1.f, 0.f));
	//this->models[3]->rotate(glm::vec3(0.f, 1.f, 0.f));
	//std::cout << this->dt<<"  " <<this->mouseOffsetX << " , " << this->mouseOffsetY << "\n";


}

void Application::render() {
	if ((application_mode == MODE::DEBUG_WITHOUT_RENDERING) or (application_mode == MODE::DEBUG_WITH_RENDERING)) {
		std::cout << "Application::render\n";
	}





	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	if (this->currentOption)
	{
		this->currentOption->OnUpdateBeg(cursorPosition.x, cursorPosition.y, cursorPosition.z);
	}

	//DRAW
	//Clear
	glClearColor(0.35f, 0.35f, 0.35f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
	//Use a program

	//ПАРАЛЛЕЛИТЬ


	if (CalculateParallelToRender) {
		if (this->dataReadyForRender) {
			this->dataIsRendering.store(true);
			if ((application_mode == MODE::RUNNING) or (application_mode == MODE::DEBUG_WITH_RENDERING)) {
				for (auto*&i : this->ComputationalDomains) {
					i->punctualColorChange(localnumber, glm::vec3(0.f));
					i->UpdateRendering(&models);
				}
			}
			this->dataIsRendering.store(false);
		}
	}
	else {
		if((application_mode== MODE::RUNNING) or (application_mode == MODE::DEBUG_WITH_RENDERING)){
			for (auto*&i : this->ComputationalDomains) {
				i->punctualColorChange(localnumber, glm::vec3(0.f));
				i->UpdateRendering(&models);
			}
		}
		std::cout << this->dataReadyForRender << "\n";
	}









	if ((application_mode == MODE::RUNNING) or (application_mode == MODE::DEBUG_WITH_RENDERING)) {
		this->shaders[SHADER_CORE_PROGRAM]->use();
		//this->shaders[RENEWABLE_SHADER]->use();
		//Update the uniforms
		this->updateUniforms();
		//Update uniforms											
		//this->materials[MAT_1]->sendToShader(*this->shaders[SHADER_CORE_PROGRAM]);
		//Activate texture
		//this->textures[TEX_GTSud]->bind(0);
		//this->textures[TEX_GTSud_SPECULAR]->bind(1);
		//Draw
		//this->meshes[MESH_QUAD]->render(this->shaders[SHADER_CORE_PROGRAM]);
		for (auto&i : this->models)
			i->render(this->shaders[SHADER_CORE_PROGRAM]);
		//for (auto&i : this->RenewableModels)
		//	i->render(this->shaders[RENEWABLE_SHADER]);
	}

	ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
	ImGui::Text("Press 'SPASE' to fix camera, press R to reset view");

	for (auto*&i : this->ComputationalDomains) {
		ImGui::Text("Particles in simulation: Real %.0f , Boundary %.0f , Virtual %.0f ", i->getGlobalStats(0), i->getGlobalStats(1), i->getGlobalStats(2));
	}
	for (auto*&i : this->ComputationalDomains) {
		ImGui::Text("Time %.5f ", i->getCurrentTime());
	}
	ImGui::InputText(" - Coloring parameter", app_coloringParam, 10);

	for (auto*&i : this->ComputationalDomains) {
		i->changeColorParamTo(app_coloringParam);

		if(i->ColoringUnknow == true) ImGui::Text("Coloring parameter: id");
		else ImGui::Text("Coloring parameter: %s", (i->ColoringParam).c_str());
		ImGui::Text("MinColor Value: %.6f, MaxColor Value: %.6f", i->m_minVal, i->m_maxVal);
	}
	for (auto*&i : this->ComputationalDomains) {
		ImGui::Text("Particle number");
		ImGui::InputInt("", &localnumber);
		ImGui::Text(i->getLocalStats(localnumber).c_str());
	}


	
	ImGui::SliderFloat3("Mouse ", &cursorPosition.x, 0.0f, 1.0f,"%.3f",0.5f);

	if (this->currentOption)
	{
		this->currentOption->OnUpdateEnd();
		ImGui::Begin("Menu");
		if (this->currentOption != this->optionMenu && ImGui::Button("<-")) {
			delete this->currentOption;
			this->currentOption = this->optionMenu;
		}
		this->currentOption->OnImGuiRender();
		ImGui::End();
	}
	





	ImGui::Render();
	ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

	//End Draw
	glfwSwapBuffers(window);
	glFlush();
	//Reset
	glBindVertexArray(0);
	glUseProgram(0);
	glActiveTexture(0);
	glBindTexture(GL_TEXTURE_2D, 0);



	//for (int i = 1; i < RenewableModels.size(); i++) {
	//	delete RenewableModels[i];
	//}
	//delete shaders[RENEWABLE_SHADER];
	//for (auto*&i : this->RenewableModels)
	//	delete i;

	//ПАРАЛЛЕЛИТЬ
	if (CalculateParallelToRender) {
		if (this->dataReadyForRender) {
			this->dataIsRendering.store(true);
			if ((application_mode == MODE::RUNNING) or (application_mode == MODE::DEBUG_WITH_RENDERING)) {
				for (auto*&i : this->ComputationalDomains) {
					i->AfterRendering(&models);
				}
			}
			this->dataIsRendering.store(false);
		}
	}
	else {
		if ((application_mode == MODE::RUNNING) or (application_mode == MODE::DEBUG_WITH_RENDERING)) {
			for (auto*&i : this->ComputationalDomains) {
				i->AfterRendering(&models);
			}
		}
	}



	if ((application_mode == MODE::DEBUG_WITHOUT_RENDERING) or (application_mode == MODE::DEBUG_WITH_RENDERING)) {
		std::cout << "\n";
	}


}

//Static functions
void Application::framebuffer_resize_callback(GLFWwindow* window, int fbW, int fbH) {
	glViewport(0, 0, fbW, fbH);
}

/*
void Application::Mouse_scrool_callback(GLFWwindow* window, double xOffset, double yOffset) {
	wheeloffset = yOffset;
}
*/

bool Application::framebuffer_resize_event(MouseScrolledEvent& e) {
	return false;
}

bool Application::Mouse_scrool_event(MouseScrolledEvent& e) {
	wheeloffset = e.GetYOffset();
	return false;
}

bool Application::Keyboard_event(KeyPressedEvent& e) {


	//Program
	if (e.GetKey() == GLFW_KEY_ESCAPE) { glfwSetWindowShouldClose(this->window, GLFW_TRUE); }
	//Save point to constructor
	if (e.GetKey() == GLFW_KEY_ENTER){
		this->toolConstructor.remember(cursorPosition.x, cursorPosition.y, cursorPosition.z);
		if (this->toolConstructor.checkSize() % 3 == 0) {

			std::vector<Mesh*> meshes;

			this->toolConstructor.addTriangleTo(&meshes, &EssentialMeshes);

			this->models.push_back(new Model(glm::vec3(0.f), this->materials[0], this->textures[TEX_GTSud], this->textures[TEX_GTSud_SPECULAR], meshes));

			for (auto*& i : meshes)
				delete i;
		}
	}
	//Save primitive to a file
	if (e.GetKey() == GLFW_KEY_P) { this->toolConstructor.savePrimitive("src/Something.h", "Something"); }
	if (e.GetKey() == GLFW_KEY_C) { this->camera.changeProjType(); }
	if (e.GetKey() == GLFW_KEY_SPACE) {
		this->camera.changeStaticMode();
		switch (this->camera.getStaticMode()) {
		case ON:
			glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
			break;
		default:
			glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
			break;
		}
	}

	//Camera and Constructor
	if (e.GetKey() == GLFW_KEY_W) {
		if (this->camera.getProjectionType() == PERSPECTIVE)
			this->camera.move_p(this->dt, FORWARD);
		else
			this->camera.move_o(this->dt, UP);
		if (this->camera.getStaticMode() == ON) {
			//this->ModelMatrix = glm::translate(this->ModelMatrix, glm::vec3(0.f, 1.f, 0.f));
			for (auto &i : this->models)
				i->move(this->dt*glm::vec3(0.f, 1.f, 0.f));
			this->EssentialMeshes[0]->move(this->dt*glm::vec3(0.f, 1.f, 0.f));
		}
	}
	if (e.GetKey() == GLFW_KEY_S) {
		if (this->camera.getProjectionType() == PERSPECTIVE)
			this->camera.move_p(this->dt, BACKWARD);
		else
			this->camera.move_o(this->dt, DOWN);
		if (this->camera.getStaticMode() == ON) {
			//this->ModelMatrix = glm::translate(this->ModelMatrix, glm::vec3(0.f, -1.f, 0.f));
			for (auto &i : this->models)
				i->move(this->dt*glm::vec3(0.f, -1.f, 0.f));
			this->EssentialMeshes[0]->move(this->dt*glm::vec3(0.f, -1.f, 0.f));
		}
	}
	if (e.GetKey() == GLFW_KEY_A){
		if (this->camera.getProjectionType() == PERSPECTIVE)
			this->camera.move_p(this->dt, LEFT);
		else
			this->camera.move_o(this->dt, LEFT);
		if (this->camera.getStaticMode() == ON) {
			//this->ModelMatrix = glm::translate(this->ModelMatrix, glm::vec3(-1.f, 0.f, 0.f));
			for (auto &i : this->models)
				i->move(this->dt*glm::vec3(-1.f, 0.f, 0.f));
			this->EssentialMeshes[0]->move(this->dt*glm::vec3(-1.f, 0.f, 0.f));
		}
	}
	if (e.GetKey() == GLFW_KEY_D){
		if (this->camera.getProjectionType() == PERSPECTIVE)
			this->camera.move_p(this->dt, RIGHT);
		else
			this->camera.move_o(this->dt, RIGHT);
		if (this->camera.getStaticMode() == ON) {
			//this->ModelMatrix = glm::translate(this->ModelMatrix, glm::vec3(1.f, 0.f, 0.f));
			for (auto &i : this->models)
				i->move(this->dt*glm::vec3(1.f, 0.f, 0.f));
			this->EssentialMeshes[0]->move(this->dt*glm::vec3(1.f, 0.f, 0.f));
		}
	}

	if (e.GetKey() == GLFW_KEY_R) { this->camera.setToDefault(); }
	if (e.GetKey() == GLFW_KEY_Q) { this->camera.move_p(this->dt, DOWN); }
	if (e.GetKey() == GLFW_KEY_E) { this->camera.move_p(this->dt, UP); }



	return false;
}


void Application::Keyboard_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {

	WindowData& data = *(WindowData*)glfwGetWindowUserPointer(window);
	//EventCallbackFn& EventCallback = *(EventCallbackFn*)glfwGetWindowUserPointer(window);
	switch (action)
	{
	case GLFW_PRESS: {
		KeyPressedEvent event = KeyPressedEvent(key, 0);
		data.EventCallback(event);
		break;
	}
	case GLFW_RELEASE: {
		KeyReleasedEvent event(key);
		//data.EventCallback(event);
		break;
	}
	case GLFW_REPEAT: {
		break;
	}
	}
}




//IMGUI Functions

void Application::initIMGUI() {

	//Setup IMGUI
	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImGuiIO& io = ImGui::GetIO(); (void)io;
	ImGui::StyleColorsDark();
	ImGui_ImplGlfw_InitForOpenGL(this->window, true);
	ImGui_ImplOpenGL3_Init((char *)glGetString(GL_NUM_SHADING_LANGUAGE_VERSIONS));
	}

void Application::optionsRegistration() {

	this->currentOption = this->optionMenu;
	this->optionMenu->RegisterTest<options::OptionColor>("Option Color");

}

void Application::CursorUpdate() {
	this->cursorPosition = this->EssentialMeshes[0]->getPosition();
	//std::cout << cursorPosition.x << " " << cursorPosition.y << " " << cursorPosition.z << "\n";
}




//COMPUTATIONAL DOMAIN
void Application::CompDomainInit() {
	SPH_OPTIONS options;
	options.nrOfParticles[PARTICLETYPE::REAL] = 5000;
	options.nrOfParticles[PARTICLETYPE::BOUNDARY] = 280;
	options.nrOfParticles[PARTICLETYPE::VIRTUAL] = 0;
	this->ComputationalDomains.push_back(new SPH_CD(options, D2));

	for (auto*&i : this->ComputationalDomains) {
		if (application_mode == MODE::RUNNING) { i->setModeToRUNNING();	}
		if ((application_mode == MODE::DEBUG_WITHOUT_RENDERING) or (application_mode == MODE::DEBUG_WITH_RENDERING)) { i->setModeToDEBUG(); }
	}

	for (auto*&i : this->ComputationalDomains) {

		i->setInitialDeltaTime(0.000005);

		std::vector<СD_Boundary*> Boundaries;
		std::vector<Mesh*> meshes;
		//meshes.push_back(new Mesh(&Quad(glm::vec3(i->getXmin(), i->getYmin(), i->getZmin()), Quad::QuadNormal(Quad::X, Quad::PLUS), (i->getYmax() - i->getYmin()), (i->getZmax() - i->getZmin())), glm::vec3(0.f, 0.f, 0.f), glm::vec3(0.f), glm::vec3(0.f), glm::vec3(1.f)));
		//meshes.push_back(new Mesh(&Quad(glm::vec3(i->getXmin(), i->getYmin(), i->getZmin()), Quad::QuadNormal(Quad::Y, Quad::PLUS), (i->getXmax() - i->getXmin()), (i->getZmax() - i->getZmin())), glm::vec3(0.f, 0.f, 0.f), glm::vec3(0.f), glm::vec3(0.f), glm::vec3(1.f)));

		meshes.push_back(new Mesh(&Quad(glm::vec3(i->getXmin(), i->getYmin(), -0.005f), Quad::QuadNormal(Quad::X, Quad::PLUS), (i->getYmax() - i->getYmin()), 0.01f), glm::vec3(0.f, 0.f, 0.f), glm::vec3(0.f), glm::vec3(0.f), glm::vec3(1.f)));
		meshes.push_back(new Mesh(&Quad(glm::vec3(i->getXmin(), i->getYmin(), -0.005f), Quad::QuadNormal(Quad::Y, Quad::PLUS), (i->getXmax() - i->getXmin()), 0.01f), glm::vec3(0.f, 0.f, 0.f), glm::vec3(0.f), glm::vec3(0.f), glm::vec3(1.f)));

		//meshes.push_back(new Mesh(&Quad(glm::vec3(i->getXmin(), i->getYmin(), i->getZmin()), Quad::QuadNormal(Quad::Z, Quad::PLUS), (i->getXmax() - i->getXmin()), (i->getYmax() - i->getYmin())), glm::vec3(0.f, 0.f, 0.f), glm::vec3(0.f), glm::vec3(0.f), glm::vec3(1.f)));

		//meshes.push_back(new Mesh(&Quad(glm::vec3(i->getXmax(), i->getYmax(), i->getZmax()), Quad::QuadNormal(Quad::X, Quad::MINUS), (i->getYmin() - i->getYmax()), (i->getZmin() - i->getZmax())), glm::vec3(0.f, 0.f, 0.f), glm::vec3(0.f), glm::vec3(0.f), glm::vec3(1.f)));
		//meshes.push_back(new Mesh(&Quad(glm::vec3(i->getXmax(), i->getYmax(), i->getZmax()), Quad::QuadNormal(Quad::Y, Quad::MINUS), (i->getXmin() - i->getXmax()), (i->getZmin() - i->getZmax())), glm::vec3(0.f, 0.f, 0.f), glm::vec3(0.f), glm::vec3(0.f), glm::vec3(1.f)));

		meshes.push_back(new Mesh(&Quad(glm::vec3(i->getXmax(), i->getYmax(), 0.005f), Quad::QuadNormal(Quad::X, Quad::MINUS), (i->getYmin() - i->getYmax()), -0.01f), glm::vec3(0.f, 0.f, 0.f), glm::vec3(0.f), glm::vec3(0.f), glm::vec3(1.f)));
		meshes.push_back(new Mesh(&Quad(glm::vec3(i->getXmax(), i->getYmax(), 0.005f), Quad::QuadNormal(Quad::Y, Quad::MINUS), (i->getXmin() - i->getXmax()), -0.01f), glm::vec3(0.f, 0.f, 0.f), glm::vec3(0.f), glm::vec3(0.f), glm::vec3(1.f)));

		//meshes.push_back(new Mesh(&Quad(glm::vec3(i->getXmax(), i->getYmax(), i->getZmax()), Quad::QuadNormal(Quad::Z, Quad::MINUS), (i->getXmin() - i->getXmax()), (i->getYmin() - i->getYmax())), glm::vec3(0.f, 0.f, 0.f), glm::vec3(0.f), glm::vec3(0.f), glm::vec3(1.f)));
		

		BOUNDARYCONDITION BCtype_List[] = { BOUNDARYCONDITION::FIRSTORDER_BC ,BOUNDARYCONDITION::FIRSTORDER_BC ,BOUNDARYCONDITION::FIRSTORDER_BC ,BOUNDARYCONDITION::FIRSTORDER_BC ,BOUNDARYCONDITION::FIRSTORDER_BC ,BOUNDARYCONDITION::FIRSTORDER_BC };
		BCPARAMETER BCparam_List[] = { BCPARAMETER::VELOCITY_X,BCPARAMETER::VELOCITY_X,BCPARAMETER::VELOCITY_X,BCPARAMETER::VELOCITY_X,BCPARAMETER::VELOCITY_X,BCPARAMETER::VELOCITY_X };
		float BCval_List[] = { 0.f, 0.0f, 0.f, 0.1f, 0.0f, 0.0f,};
		int count = 0;

		//Boundaries.push_back(new СD_Boundary(meshes[0], meshes[2]));
		Boundaries.push_back(new СD_Boundary(meshes[0], BOUNDARYCONDITION::FIRSTORDER_BC, BCPARAMETER::VELOCITY_X, 0.0f));   //ЛЕВО

		Boundaries.push_back(new СD_Boundary(meshes[1], BOUNDARYCONDITION::FIRSTORDER_BC, BCPARAMETER::VELOCITY_X, 0.0));   //НИЗ

		//Boundaries.push_back(new СD_Boundary(meshes[2], meshes[0]));
		Boundaries.push_back(new СD_Boundary(meshes[2], BOUNDARYCONDITION::FIRSTORDER_BC, BCPARAMETER::VELOCITY_X, 0.0f));   //ПРАВО
		
		Boundaries.push_back(new СD_Boundary(meshes[3], BOUNDARYCONDITION::FIRSTORDER_BC, BCPARAMETER::VELOCITY_X, 1.0f));   //ВЕРХ


		for (auto* j : meshes) {
			//Boundaries.push_back(new СD_Boundary(j, BCtype_List[count], BCparam_List[count], BCval_List[count]));
			//count++;
			//j->printAll();
			delete j;
		}


		for (int j = 0;j < 4;j++) {
			i->addBoundary(Boundaries[j]);
		}

		if ((application_mode == MODE::RUNNING) or (application_mode == MODE::DEBUG_WITH_RENDERING)) {
			BoundaryMeshing(i);
		}
		i->Initilization();

		//if ((application_mode == MODE::RUNNING) or (application_mode == MODE::DEBUG_WITH_RENDERING)) {
		//	CDMeshing(i);
		//}
	}


}

void Application::BoundaryMeshing(CompDomain* CDptr) {

	//BOUNDARYMESHING
	std::vector<Mesh*> B_meshes;
	CDptr->InitialBoundaryRendering(&B_meshes);
	this->models.push_back(new Model(glm::vec3(0.f), this->materials[0], this->textures[TEX_GTSud], this->textures[TEX_GTSud_SPECULAR], B_meshes));
	CDptr->assignBoundaryModel(models.size() - 1);
	for (auto*& i : B_meshes)
		delete i;


}

/*
void Application::CDMeshing(CompDomain* CDptr) {
	//models.resize(models.size() - 2);
	//models.erase(models.begin()+1,models.begin()+ models.size());
	std::vector<Mesh*> meshes;
	CDptr->InitialRendering(&meshes);
	this->models.push_back(new Model(glm::vec3(0.f), this->materials[0], this->textures[TEX_GTSud], this->textures[TEX_GTSud_SPECULAR], meshes));
	CDptr->assignModel(models.size() - 1);
	for (auto*& i : meshes)
		delete i;
}
*/

