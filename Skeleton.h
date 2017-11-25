#pragma once

#include "Geometry.h"
#include "RootAttributes.h"
#include "SkeletonEdge.h"


#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <sstream>
#include <map>
#include <vector>

#include "boost/algorithm/string.hpp"
#include "boost/lexical_cast.hpp"


namespace Roots
{

	typedef std::vector<SkeletonEdge> edgeList;
	typedef edgeList::iterator edgeIter;

	typedef std::vector<SkeletonEdge*> edgePtrList;
	typedef edgePtrList::iterator edgePtrIter;

	typedef std::vector<Point3d> vertList;
	typedef vertList::iterator vertIter;



	class Skeleton
	{
	public:
		//Serializable Implementation >> //
		//Skeleton(Json::Value json);
		//Json::Value ToJson();

		static std::string beginSkeletonString;
		static std::string endSkeletonString;
		static std::string vertexString;
		static std::string edgeString;
		static std::string endHeaderString ;

		friend std::ostream& operator<<(std::ostream& out, const Skeleton &skel);
		friend std::istream& operator>>(std::istream& in, Skeleton &skel);

		// << Serializable Implementation//

		Skeleton();
		Skeleton(vertList aVertices, std::vector<edgePtrList> aEdges);

		void LoadFromTextFile(char *filename);
		void LoadWenzhenStyleFile(std::istream& in);
		//void LoadFromJson(Json::Value skelJson);
		//void LoadFromJsonFile(char *filename);
		

		vertList getVertices();
		std::vector<edgePtrList> getEdges();


		void AddEdge(SkeletonEdge toAdd);

		void RemoveEdges(edgeList toRemove, bool careAboutAttributes = false);
		void RemoveEdge(SkeletonEdge toRemove, bool careAboutAttributes = false);
		void RemoveEdge(int v0, int v1);

		SkeletonEdge* GetEdge(int v0, int v1);

		
		//static Skeleton FromJson();


	private:

		//vertex list, potentially sorted in some way, but for now, nada
		vertList mVerts;

		//list of edges sorted by the id of the first vertex, then the second
		std::vector<edgePtrList> mEdges;
		
	public:
		std::vector<std::vector<int>> mNeighbors;
		int mNumVertices;
		int mNumEdges;
	};

	/*
	This function will take a single root (edgePtrList) and order that list of edgePtrs from beginning to end.  If no edge with a degree 1 vertex is provided as
	the forward end, then the returned edge-list will be ordered by the lower degree 1 vertex id.  Will return the 
	*/
	edgePtrList tidyUpRoot(edgePtrList toTidy, bool &valid, int &end1, int &end2, SkeletonEdge *beginning = nullptr);

}