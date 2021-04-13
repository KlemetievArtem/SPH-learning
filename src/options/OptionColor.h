#pragma once
#include "option.h"

namespace options {
	class OptionColor : public Option {
	public:
		OptionColor() {}
		~OptionColor() {}

		void outside() override;
		void InitModels(std::vector<Mesh*>* meshes, std::vector<Model*>* models) override;
		void OnUpdateBeg(double offsetX, double offsetY, double wheelOffset) override;
		void OnUpdateEnd() override;
		virtual void OnImGuiRender() override;

		bool getFlag() override { return firstcycle; }
		void setFlag(bool flag) override { firstcycle = flag; }

		bool firstcycle = true;
		std::unique_ptr<Mouse> mouse;
	};
}
