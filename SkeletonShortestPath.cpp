#include "Skeleton.h"

namespace Roots
{
	namespace
	{
		struct WeightedEdge
		{
			SkeletonEdge edge;
			float AStarWeight;
			
			WeightedEdge(SkeletonEdge edge, Point3d from, Point3d to, float AStarprevious)
			{

			}
		};
	}
}