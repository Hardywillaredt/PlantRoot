#pragma once
#include "RootAttributes.h"
#include "Listener.h"
#include <set>
#include <string>

namespace Roots
{
	struct SkeletonEdge
	{
	public:
		
		int v0, v1;
		RootAttributes attributes;
		

		SkeletonEdge();

		SkeletonEdge(int firstVert, int secondVert, RootAttributes aAttributes = RootAttributes());

		SkeletonEdge(int firstVert, int secondVert, double * attributeData);

		SkeletonEdge(const SkeletonEdge& other);

		//Json::Value ToJson()
		//{
		//	Json::Value edge = Json::Value(Json::ValueType::arrayValue);
		//	edge.resize(2 + RootAttributes::NumAttributes);
		//	edge[0] = v0;
		//	edge[1] = v1;
		//	for (int att = 0; att < RootAttributes::NumAttributes; ++att)
		//	{
		//		edge[att + 2] = attributes[att];
		//	}
		//	return edge;
		//}
		
		friend std::ostream& operator<<(std::ostream& out, const SkeletonEdge& edge);
		friend std::istream& operator>>(std::istream& in, SkeletonEdge& edge);

	/*	std::ostream& operator<<(std::ostream& out)
		{
			out << v0 << " " << v1 << " ";
			
			out << attributes;
		}

		std::istream& operator>>(std::istream& in)
		{
			in >> v0;
			in >> v1;
			for (int att = 0; att < RootAttributes::NumAttributes; ++att)
			{
				in >> attributes[att];
			}
		}*/

		/*SkeletonEdge(Json::Value jsonEdge)
		{
			v0 = jsonEdge[0].asInt();
			v1 = jsonEdge[1].asInt();

			int tempV;
			if (v0 > v1)
			{
				tempV = v1;
				v1 = v0;
				v0 = tempV;
			}

			double *attributeData = new double[RootAttributes::NumAttributes];

			for (int i = 0; i < RootAttributes::NumAttributes; ++i)
			{
				attributeData[i] = jsonEdge[i + 2].asDouble();
			}
			attributes = RootAttributes(attributeData);

			delete[] attributeData;

		}*/

		~SkeletonEdge();

		static std::string ClassName();

		void RegisterListener(Listener *listener);

		void RemoveListener(Listener *listener);

		void update(RootAttributes updatedAttributes);

		bool operator<(SkeletonEdge &other);

		friend bool edgePtrLessThan(SkeletonEdge *first, SkeletonEdge *second);

		//this operator ignores the attribution information for convenience in the edge finding and sorting process
		bool operator==(SkeletonEdge &other);

		bool operator!=(SkeletonEdge &other);


	private:
		std::set<Listener*> mListeners;
	};

}
