#include "Skeleton.h"

namespace Roots
{
	namespace
	{
		struct WeightedEdge
		{
			SkeletonEdge edge;
			double AStarWeight;
			
			WeightedEdge(SkeletonEdge edge, Point3d from, Point3d to, double AStarprevious)
			{

			}
		};
	}
}