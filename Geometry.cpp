#include "Geometry.h"



Point3d::Point3d(double *aData)
{
	data[0] = aData[0];
	data[1] = aData[1];
	data[2] = aData[2];
}

Point3d::Point3d(double x, double y, double z)
{
	data[0] = x;
	data[1] = y;
	data[2] = z;
}

Point3d::Point3d()
{
	data[0] = 0;
	data[1] = 0;
	data[2] = 0;
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

double * Point3d::getData()
{
	return data;
}

double Point3d::x() 
{ 
	return data[0]; 
}
double Point3d::y()
{
	return data[1]; 
}
double Point3d::z() 
{
	return data[2]; 
}

Point3d operator-(Point3d &first, Point3d &second)
{
	double resultData[3];
	for (int i = 0; i < 2; ++i)
	{
		resultData[i] = first.data[i] - second.data[i];
	}
	return Point3d(resultData);
}

Point3d operator+(Point3d &first, Point3d &second)
{
	double resultData[3];
	for (int i = 0; i < 2; ++i)
	{
		resultData[i] = first.data[i] + second.data[i];
	}
	return Point3d(resultData);
}

std::ostream& operator<<(std::ostream& out, const Point3d& point)
{
	out << point.data[0] << " " << point.data[1] << " " << point.data[2] << std::endl;
	return out;
}

std::istream& operator>>(std::istream& in, Point3d& point)
{
	in >> point.data[0];
	in >> point.data[1];
	in >> point.data[2];
	return in;
}

