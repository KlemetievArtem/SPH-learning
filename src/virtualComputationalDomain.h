#pragma once

#include "Boundary.h"

enum DIMENSIONS {
	D0,
	D1 = 1,
	D2,
	D3
};


/*
class CD_boundary {
private:
	std::vector<Primitive*> PrimitivesList;
public:
	CD_boundary() {}
	void addBoundary(std::vector<Primitive*> Plist) {
		for (auto*i : newPlist) {
			PrimitivesList.push_back(new Primitive(*i));
		}
	}
	~CD_boundary() {
		for (auto*&i : this->PrimitivesList)
			delete i;

	}
};
*/


enum MICROCOSME
{
	MC_SPH
};

struct CompDomainStatystic {
	float maxVelocity;
	float minScale;

	float aveScale;
};

enum MODE_CD{CD_RUNNING,CD_DEBUG};

class CompDomain {
private:
	

	float xmin = 0.f;
	float xmax = 1.f;
	float ymin = 0.f;
	float ymax = 1.f;
	float zmin = 0.f;
	float zmax = 1.f;


	cd_prec initTime = 0.00;
	cd_prec deltaTime;
	cd_prec initialdeltaTime;
	cd_prec finshTime;
	cd_prec currentTime = 0.00;

	MICROCOSME MCType;

	CompDomainStatystic m_statystic;


protected:
	std::vector<ÑD_Boundary*> Boundaries;
	int ModelId;
	int CD_ModelId;
	int MaterialId;
	int BoundaryModelId;
	MODE_CD computationalDomain_mode;
	DIMENSIONS nrOfDim;

public:



	//COLORING
	float m_maxVal;
	float m_minVal;
	std::string ColoringParam;
	bool ColoringUnknow = true;
	void changeColorParamTo(std::string name) {
		ColoringParam = name;
	}
	std::string getColorParam() {
		return ColoringParam;
	}


	void setXminTo(float const val) { xmin = val; }
	void setXmaxTo(float const val) { xmax = val; }
	void setYminTo(float const val) { ymin = val; }
	void setYmaxTo(float const val) { ymax = val; }
	void setZminTo(float const val) { zmin = val; }
	void setZmaxTo(float const val) { zmax = val; }


	void setTypeTo(MICROCOSME type) { MCType = type; }
	void setModeToDEBUG() { computationalDomain_mode = MODE_CD::CD_DEBUG; }
	void setModeToRUNNING() { computationalDomain_mode = MODE_CD::CD_RUNNING; }

	float getXmin() { return xmin; }
	float getXmax() { return xmax; }
	float getYmin() { return ymin; }
	float getYmax() { return ymax; }
	float getZmin() { return zmin; }
	float getZmax() { return zmax; }


	cd_prec getDeltaTime() { return deltaTime; }
	void setDeltaTime(cd_prec dt) { deltaTime = dt; }
	cd_prec getInitialDeltaTime() { return initialdeltaTime; }
	void setInitialDeltaTime(cd_prec dt) { initialdeltaTime = dt; deltaTime = dt; }
	cd_prec getCurrentTime() { return currentTime; }
	void timeStepEnd() { currentTime += deltaTime; }


	void addBoundary(ÑD_Boundary* boundary) {
		Boundaries.push_back(boundary);
	}
	void InitialBoundaryRendering(std::vector<Mesh*>* meshes) {

		//std::cout << "\n" << Boundaries.size() << "\n";
		for (auto i : Boundaries) {
			Mesh* testMesh = i->ReturnMesh();
			//testMesh->printAll();
			meshes->push_back(new Mesh(i->ReturnMesh()));
		}
	}
	void assignBoundaryModel(int id) { BoundaryModelId = id; }

	virtual void Initilization() = 0;
	//virtual void InitialRendering(std::vector<Mesh*>* meshes) = 0;

	virtual void assignPresetModels(int id) { ModelId = id; }
	virtual void assignNewModels(int id) { CD_ModelId = id; }
	virtual void assignPresetMaterials(int id) { MaterialId = id; }
	//virtual void UpdateRendering(std::vector<Model*>* models) = 0;
	virtual void UpdateRendering(std::vector<Model*>* models, Texture* tex, Texture* tex_specualar, std::vector<Material*>* materials) = 0;
	virtual void AfterRendering(std::vector<Model*>* models) = 0;


	virtual void timeStep(cd_prec dt) = 0;
	virtual void timeStep_thread(cd_prec dt, std::atomic<bool>& dataReadyForRender, std::atomic<bool>& dataIsRendering) = 0;





	void setStatMaxVel(float val) { m_statystic.maxVelocity = val; }
	void setStatMinScale(float val) { m_statystic.minScale = val; }
	cd_prec getStatMinTime() { return m_statystic.minScale / m_statystic.maxVelocity; }

	void setStatAveScale(float val) { m_statystic.aveScale = val; }

	virtual float getGlobalStats(int number) = 0;
	virtual std::string getLocalStats(int number) = 0;

	virtual ~CompDomain() {

	}


	virtual void punctualColorChange(int number, glm::vec3 color){}

};
