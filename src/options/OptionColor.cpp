#include "OptionColor.h"

void options::OptionColor::outside() {

	//mouse = std::make_unique<Mouse>();
}

void options::OptionColor::InitModels(std::vector<Mesh*>* meshes, std::vector<Model*>* models) {

}

void options::OptionColor::OnUpdateBeg(double offsetX,double offsetY,double wheelOffset) {

}

void options::OptionColor::OnUpdateEnd() {
	//glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	//glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

}

void options::OptionColor::OnImGuiRender() {
	//ImGui::Text("Mouse.x %d", static_cast<const int> (mouse->position.x));
}
