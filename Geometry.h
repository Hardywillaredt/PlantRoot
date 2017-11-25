#pragma once
#include <iostream>
//#include "json\json.h"
//#include "Serializable.h"


struct Point3d
{
public:
	//serializable implementation >> //
	//Point3d(Json::Value json);
	//Json::Value ToJson();
	// << serializable implementation//


	Point3d(double *aData);

	Point3d(double x, double y, double z);

	Point3d();

	double *getData();

	double x(); 
	double y();
	double z();

	friend Point3d operator-(Point3d &first, Point3d &second);

	friend Point3d operator+(Point3d &first, Point3d &second);

	friend std::ostream& operator<<(std::ostream& out, const Point3d& point);

	friend std::istream& operator>>(std::istream& in, Point3d& point);

private:
	double data[3];
};


