#pragma once

#include "libs.h"
#include <functional>
#include "imgui/imgui.h"

namespace options {
	class Option {
	public:
		Option() {}
		virtual ~Option() {}

		virtual void outside() {}
		virtual void InitModels(std::vector<Mesh*>* meshes, std::vector<Model*>* models) {}
		virtual void OnUpdateBeg(double offsetX, double offsetY, double wheelOffset) {}
		virtual void OnUpdateEnd() {}
		virtual void OnImGuiRender() {}
		virtual bool getFlag() { return true; }
		virtual void setFlag(bool flag) {}
	};

	class OptionMenu : public Option {
	public:
		OptionMenu(Option*& currentOptionPointer) 
			:m_CurrentOption(currentOptionPointer) {}

		void OnImGuiRender() override {
			for (auto& option : m_Options) {
				if (ImGui::Button(option.first.c_str()))
					m_CurrentOption = option.second();
			}
		}
		template<class T>
		void RegisterTest(const std::string& name) {
			std::cout << "Registering option " << name << '\n';
			m_Options.push_back(std::make_pair(name, []() { return new T(); }));
		}
	private:
		Option*& m_CurrentOption;
		std::vector<std::pair<std::string, std::function<Option*()>>> m_Options;
	};

}