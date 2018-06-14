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
		float euclidLength;
		int v0id, v1id;
		float thickness, width, length;
		int glV0, glV1;
		//boost::python::list data;



		friend std::ostream& operator<<(std::ostream &out, const RootAttributes &attribs);
		friend std::istream& operator>>(std::istream &in, RootAttributes &attribs);

		RootAttributes();
		RootAttributes(float *attributeData, Point3d v0 = Point3d(), Point3d v1 = Point3d());
		RootAttributes(std::vector<float> attributeData, Point3d v0 = Point3d(), Point3d v1 = Point3d());
		RootAttributes(float aThickness, float aWidth, Point3d v0 = Point3d(), Point3d v1 = Point3d());

		bool operator==(RootAttributes& second);

		float operator[](const int index);


	
	};
}
