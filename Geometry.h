#pragma once

#include <iostream>
#include "boost/python.hpp"


struct Point3d
{
public:
	
	float p[6];

	int id;

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
