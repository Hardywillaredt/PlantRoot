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
	p[0] = aData[0];
	p[1] = aData[1];
	p[2] = aData[2];
	id = index;
}


Point3d::Point3d(float ax, float ay, float az, int index)
{
	p[0] = ax;
	p[1] = ay;
	p[2] = az;
	id = index;
}

Point3d::Point3d()
{
	p[0] = 0;
	p[1] = 0;
	p[2] = 0;
	id = 0;
}

Point3d::Point3d(const Point3d &other)
{
	for (int i = 0; i < 3; ++i)
	{
		p[i] = other.p[i];
	}
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
	return Point3d(first[0] - second[0], first[1] - second[1], first[2] - second[2]);
}

Point3d operator+(Point3d &first, Point3d &second)
{

	return Point3d(first[0] + second[0], first[1] + second[1], first[2] + second[2]);
}

Point3d operator+=(Point3d &first, Point3d &second)
{
	first = first + second;
	return first;
}

std::ostream& operator<<(std::ostream& out, const Point3d& point)
{
	out << point.p[0] << " " << point.p[1] << " " << point.p[2];
	return out;
}

std::istream& operator>>(std::istream& in, Point3d& point)
{
	in >> point[0];
	in >> point[1];
	in >> point[2];
	return in;
}

Point3d operator/(Point3d & p, float div)
{
	return Point3d(p[0] / div, p[1] / div, p[2] / div);
}

Point3d operator*(Point3d & p, float mul)
{
	return Point3d(p[0]*mul, p[1]*mul, p[2]*mul);
}



float Point3d::mag()
{
	float result = 0.0;
	result += p[0]*p[0] + p[1]*p[1] + p[2]*p[2];
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
