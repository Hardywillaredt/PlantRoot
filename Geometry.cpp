#include "Geometry.h"
#include <fstream>
#include <iostream>
#include <string>
#include <time.h>

std::string logFile = "D:/dev/school/rootsProject/python/log.txt";
std::ofstream Log::out = std::ofstream();
bool Log::isLogOpen = false;

Point3d::Point3d(float *aData, int index)
{
	x = aData[0];
	y = aData[1];
	z = aData[2];
	id = index;
}


Point3d::Point3d(float ax, float ay, float az, int index)
{
	x = ax;
	y = ay;
	z = az;
	id = index;
}

Point3d::Point3d()
{
	x = 0;
	y = 0;
	z = 0;
	id = 0;
}

Point3d::Point3d(const Point3d &other)
{

	x = other.x;
	y = other.y;
	z = other.z;
	id = other.id;
}

//Json::Value Point3d::ToJson()
//{
//	Json::Value edge = Json::Value(Json::ValueType::arrayValue);
//	edge.resize(3);
//
//	edge[0] = data[0];
//	edge[1] = data[1];
//	edge[2] = data[2];
//
//	return edge;
//}



Point3d operator-(Point3d &first, Point3d &second)
{
	return Point3d(first.x - second.x, first.y - second.y, first.z - second.z);
}

Point3d operator+(Point3d &first, Point3d &second)
{

	return Point3d(first.x + second.x, first.y + second.y, first.z + second.z);
}

Point3d operator+=(Point3d &first, Point3d &second)
{
	first = first + second;
	return first;
}

std::ostream& operator<<(std::ostream& out, const Point3d& point)
{
	out << point.x << " " << point.y << " " << point.z;
	return out;
}

std::istream& operator>>(std::istream& in, Point3d& point)
{
	in >> point.x;
	in >> point.y;
	in >> point.z;
	return in;
}

Point3d operator/(Point3d & p, float div)
{
	return Point3d(p.x / div, p.y / div, p.z / div);
}

Point3d operator*(Point3d & p, float mul)
{
	return Point3d(p.x*mul, p.y*mul, p.z*mul);
}

float Point3d::mag()
{
	float result = 0.0;
	result += x*x + y*y + z*z;
	return std::sqrt(result);
}

Log::Log()
{
}

Log::~Log()
{
	out.close();
}

void Log::WriteLine(std::string line)
{
	if (!isLogOpen)
	{
		out.open(logFile, std::ios_base::app);
		isLogOpen = true;

		time_t theTime = time(NULL);
		struct tm *aTime = localtime(&theTime);


		out << aTime->tm_hour << ":" << aTime->tm_min << std::endl;
	}
	out << line << std::endl;
}

void Log::CloseLog()
{
	isLogOpen = false;
	out.close();
}
