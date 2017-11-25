#pragma once
#include "Geometry.h"
////////////////////////////NO LONGER A THING///////////////////////////
namespace Roots
{
	class Node
	{
		Node(Point3d aNodeCenter = Point3d(), int degree = 0);

		void IncrementDegree();

	private:
		Point3d mNodeCenter;
		int mDegree;
	};
}