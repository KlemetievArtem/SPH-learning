#pragma once

#include "Camera.h"
#include "Events/MouseEvent.h"
#include "Events/KeyEvent.h"

#include "options/OptionColor.h"

#include "Something.h"


#include "SPH/ComputationalDomain.h"



#define GLCORE_BIND_EVENT_FN(fn) std::bind(&fn, this, std::placeholders::_1)



//ENUMERATIONS
enum shader_enums{SHADER_CORE_PROGRAM=0, RENEWABLE_SHADER};
enum texture_enum{TEX_GTS=0, TEX_GTS_SPECULAR, TEX_GTSud, TEX_GTSud_SPECULAR};
enum material_enum{MAT_1=0};
enum MESH_ENUM { MESH_QUAD = 0 };

enum MODE {RUNNING, DEBUG_WITH_RENDERING, DEBUG_WITHOUT_RENDERING};

class Application
{
private:
//
	MODE application_mode = DEBUG_WITH_RENDERING;



//Variables
	//Window
	GLFWwindow* window;
	const int WINDOW_WIDTH;
	const int WINDOW_HEIGHT;
	int framebufferWidth;
	int framebufferHeight;
	//OpenGL Context
	const int GL_VERSION_MAJOR;
	const int GL_VERSION_MINOR;

	//Delta time
	float dt;
	float curTime;
	float lastTime;
	//Mouse input
	double lastMouseX;
	double lastMouseY;
	double mouseX;
	double mouseY;
	double mouseOffsetX;
	double mouseOffsetY;
	bool firstMouse;
	double lastWheel;

	//Camera
	Camera camera;

	//Matrices
	glm::mat4 ViewMatrix;
	glm::vec3 camPosition;
	glm::vec3 worldUp;
	glm::vec3 camFront;
	glm::mat4 ProjectionMatrix;
	float fov;
	float nearPlane;
	float farPlane;

	glm::mat4 ModelMatrix;

	//Shaders
	std::vector<Shader*> shaders;

	//Textures
	std::vector<Texture*> textures;

	//Materials
	std::vector<Material*> materials;

	//Meshes
	std::vector<Mesh*> EssentialMeshes;

	//Models
	std::vector<Model*> models;
	std::vector<Model*> RenewableModels;

	//Lights
	std::vector<glm::vec3*> ligts;
	//LIGHT MOVEMENT
	float light_fi=0.f;
	float light_radius=10.f;





	std::vector<int> s_Presscount;




//Private functions
	void initGLFW();
	void initWindow(const char*title, bool resizeable);
	void initGLEW(); //AFTER CONTEXT CREATION
	void initOpenGLOptions();
	void initMatrices();

	void initShaders();
	void initTextures();
	void initMaterials();
	void initOBJModels();
	void initCursorMesh();
	void initModels();
	void initLights();

	void initUniforms();
	void updateUniforms();
//Static variables
	static double wheeloffset;



//Cursor
	glm::vec3 cursorPosition;


public:
//Constructord and distructors
	Application(const char*title,
		const int WINDOW_WIDTH, const int WINDOW_HEIGHT,
		const int GL_VERSION_MAJOR, const int GL_VERSION_MINOR, bool resizable);
	virtual ~Application();

//Accessors
	int getWindowShouldClouse();
	GLFWwindow* getWindow() { return window; };
	glm::vec3 getAllMouseOffSet() { return glm::vec3(mouseOffsetX, mouseOffsetY, wheeloffset); }
//Modifiers
	void setWindowShouldClose();
//Functions
	void updateDt();
	void updateMouseInput();
	void updateKeyboardInput();
	void updateInput();
	void update();
	void render();

//Static functions
	static void framebuffer_resize_callback(GLFWwindow* window, int fbW, int fbH);
	//static void Mouse_scrool_callback(GLFWwindow* window, double xOffset, double yOffset);
	static void Keyboard_callback(GLFWwindow* window, int key, int scancode, int action, int mods);
	/*
	static void updateInput(GLFWwindow * window);
	static void updateInput(GLFWwindow * window, Mesh & mesh);
	*/

	bool framebuffer_resize_event(MouseScrolledEvent& e);
	bool Mouse_scrool_event(MouseScrolledEvent& e);

	bool Keyboard_event(KeyPressedEvent & e);

	

	using EventCallbackFn = std::function<void(Event&)>;

	struct WindowData
	{
		EventCallbackFn EventCallback;
	};
	WindowData m_Data;


	void SetEventCallback(const EventCallbackFn& callback)  { m_Data.EventCallback = callback; }
	
	
	//IMGUI PARAMS
	int localnumber = -1;




	//IMGUI Functions


	void initIMGUI();

	//Options
	options::Option* currentOption = nullptr;
	options::OptionMenu* optionMenu = new options::OptionMenu(currentOption);


	//std::unique_ptr<options::Option> currentOption;
	//std::unique_ptr<options::OptionMenu> optionMenu;

	void optionsRegistration();



	void CursorUpdate();


	ToolConstructor toolConstructor;


	void EmptyEventFunction(Event& e) { 
		std::cout << e.ToString(); 
		EventDispatcher dispatcher(e);
		dispatcher.Dispatch<MouseScrolledEvent>(GLCORE_BIND_EVENT_FN(Application::Mouse_scrool_event));
		dispatcher.Dispatch<KeyPressedEvent>(GLCORE_BIND_EVENT_FN(Application::Keyboard_event));

	}




//COMPUTATIONAL DOMAIN
	std::vector<CompDomain*> ComputationalDomains;

	void CompDomainInit();


	void BoundaryMeshing(CompDomain* CDptr);

	void CDMeshing(CompDomain* CDptr);


















};


