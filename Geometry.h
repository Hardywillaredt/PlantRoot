#pragma once

#include <iostream>
#include "boost/python.hpp"


struct Point3d
{
public:
	
	float p[6];

	int id;

	float lambdaOdd = 0.5;
	float lambdaEven = 1 / (0.1 - 1 / 0.5);

	float curvature=-1;

	Point3d(float ax, float ay, float az, float thickness, float width, int index = 0);

	Point3d();

	Point3d(const Point3d& other);

	friend Point3d operator-(Point3d &first, Point3d &second);

	friend Point3d operator+(Point3d &first, Point3d &second);

	friend Point3d operator+=(Point3d &first, Point3d &second);

	friend std::ostream& operator<<(std::ostream& out, const Point3d& point);

	friend std::istream& operator>>(std::istream& in, Point3d& point);

	friend Point3d operator/(Point3d &p, float div);

	friend Point3d operator*(Point3d &p, float mul);

	bool operator==(Point3d &other)
	{
		for (int i = 0; i < 3; ++i)
		{
			if (other.p[i] != p[i])
			{
				return false;
			}
		}
		return true;
	}

	bool operator!=(Point3d &other)
	{
		return !operator==(other);
	}

	bool operator<(const Point3d &other) const
	{
		return id < other.id;
	}
	float& operator[](size_t i)
	{
		return p[i];
	}
	float x()
	{
		return p[0];
	}
	//void fairX(float x1, float x2, float x, int i) {
	//	if (i % 2 == 0) {
	//		p[0] = (1-lambdaEven)*x+lambdaEven*(x1 + x2) / 2;
	//	}
	//	else {
	//		p[0]= (1 - lambdaOdd)*x + lambdaOdd * (x1 + x2) / 2;
	//	}
	//	
	//}
	//void fairY(float y1, float y2,float y, int i) {
	//	if (i % 2 == 0) {
	//		p[1] = (1 - lambdaEven)*y + lambdaEven * (y1 +y2) / 2;
	//	}
	//	else {
	//		p[1] = (1 - lambdaOdd)*y + lambdaOdd * (y1 + y2) / 2;
	//	}
	//}
	//void fairZ(float z1, float z2, float z, int i) {
	//	if (i % 2 == 0) {
	//		p[2] = (1 - lambdaEven)*z + lambdaEven * (z1 + z2) / 2;
	//	}
	//	else {
	//		p[2] = (1 - lambdaOdd)*z + lambdaOdd * (z1 +z2) / 2;
	//	}
	//}
	float y()
	{
		return p[1];
	}
	float z()
	{
		return p[2];
	}
	float thickness()
	{
		return p[3];
	}
	float width()
	{
		return p[4];
	}
	float ratio()
	{
		return p[5];
	}
	float getCurvature() {
		return curvature;
	}
	void setCurvature(float cur) {
		curvature = cur;
	}
	float mag();
private:

};

struct Log
{
	static std::ofstream out;
	static bool isLogOpen;
	Log();
	static void WriteLine(std::string line);
	static void CloseLog();
	
	~Log();
};
