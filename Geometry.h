#pragma once
#include <iostream>
#include "boost/python.hpp"
//#include "json\json.h"
//#include "Serializable.h"


struct Point3d
{
public:
	//serializable implementation >> //
	//Point3d(Json::Value json);
	//Json::Value ToJson();
	// << serializable implementation//
	float p[3];
	int id;

	Point3d(float *aData, int index = 0);

	Point3d(float ax, float ay, float az, int index = 0);

	Point3d();

	Point3d(const Point3d& other);

	friend Point3d operator-(Point3d &first, Point3d &second);

	friend Point3d operator+(Point3d &first, Point3d &second);

	friend Point3d operator+=(Point3d &first, Point3d &second);

	friend std::ostream& operator<<(std::ostream& out, const Point3d& point);

	friend std::istream& operator>>(std::istream& in, Point3d& point);

	friend Point3d operator/(Point3d &p, float div);

	friend Point3d operator*(Point3d &p, float mul);

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
	float y()
	{
		return p[1];
	}
	float z()
	{
		return p[2];
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
