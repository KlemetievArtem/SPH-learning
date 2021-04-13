#pragma once

/*
class KernelPart {
private:
	std::map<int, double> PolinominalPart;
	std::map<int, double> ExponentialPart;
	double minRange, maxRange;
public:
	KernelPart(double min = 0.0, double max =1.0): minRange(min), maxRange(max) {
		PolinominalPart.insert({ 0,0 });
		ExponentialPart.insert({ 0,0 });
	}
	void addPolinominalPart(int degree, double factor) {
		PolinominalPart.insert({ degree,factor });
	}
	void addExponentialPart(int degree, double factor) {
		ExponentialPart.insert({ degree,factor });
	}

	void clearPolinominalPart() { PolinominalPart.clear(); }
	void clearExponentialPart() { ExponentialPart.clear(); }
	void changeRanges(double min, double max) { minRange = min; maxRange = max; }

	std::vector<int> returnPolinominalDegrees() {
		std::vector<int> retVal;
		for (auto i : PolinominalPart) {
			retVal.push_back(i.first);
		}
		return retVal;
	}
	std::vector<int> returnExponentialDegrees() {
		std::vector<int> retVal;
		for (auto i : ExponentialPart) {
			retVal.push_back(i.first);
		}
		return retVal;
	}

	~KernelPart() {
		clearPolinominalPart();
		clearExponentialPart();
	}

};


class Kernel {
private:
	std::vector<KernelPart> kernelParts;
	double maxRange;
public:
	Kernel(int r): maxRange(r){}
	void addKernelPart(double minr, double maxr) {
		kernelParts.push_back(KernelPart(minr, maxr));
	}
	~Kernel() { kernelParts.clear(); }

};
*/

enum DiffAxis {
	X,
	Y,
	Z,
	R
};

