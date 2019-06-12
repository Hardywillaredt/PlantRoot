#pragma once

#include <iostream>
#include <vector>
#include <sstream>
#include "Geometry.h"


namespace Roots
{

	enum Attributes
	{
		Thickness = 0,
		Width=1,
		Length=2,
		NumAttributes=3
	};

	struct RootAttributes
	{

	public:
		int v0id, v1id;
		float euclidLength;

		friend std::ostream& operator<<(std::ostream &out, const RootAttributes &attribs);
		friend std::istream& operator>>(std::istream &in, RootAttributes &attribs);

		RootAttributes();
		RootAttributes(float *attributeData, Point3d v0 = Point3d(), Point3d v1 = Point3d());
		RootAttributes(std::vector<float> attributeData, Point3d v0 = Point3d(), Point3d v1 = Point3d());
		RootAttributes(float aThickness, float aWidth, Point3d v0 = Point3d(), Point3d v1 = Point3d());
		RootAttributes(int id0, int id1, Point3d &v0, Point3d &v1);

		bool operator==(RootAttributes& second);
	};
}
