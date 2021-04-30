#pragma once

template<class T>
void inline safeDelete(T** pointer)
{
	if (*pointer != nullptr)
	{
		delete *pointer;
		*pointer = nullptr;
	}
}

struct LinearFunc {
	float k,b;
};

LinearFunc inline LinearFuncCoefficients(glm::vec2 p1, glm::vec2 p2) {
	LinearFunc temp;
	if ((p2.x - p1.x) ==0)
		temp.k = (p2.y - p1.y) / (p2.x - p1.x*1.000001f);
	else
		temp.k = (p2.y - p1.y) / (p2.x - p1.x);
	temp.b = p1.y- temp.k*p1.x;
	return temp;
}

