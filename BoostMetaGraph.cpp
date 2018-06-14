#include "BoostMetaGraph.h"
#include "boost/graph/connected_components.hpp"
#include "boost/graph/kruskal_min_spanning_tree.hpp"
#include <iostream>
#include <fstream>
namespace
{
	GLfloat defaultColor[4] = { 1.0, 0.0, 0.0, 1.0 };
}

namespace Roots
{
	drawing::VBOSphere BMetaGraph::drawSphere = drawing::VBOSphere();
	drawing::VBOCube BMetaGraph::drawCube = drawing::VBOCube();
	float MinMaxStruct::minThickness = 10000;
	float MinMaxStruct::maxThickness = -1;

	float MinMaxStruct::minWidth = 10000;
	float MinMaxStruct::maxWidth = -1;
	
	float MinMaxStruct::minRatio = 10000;
	float MinMaxStruct::maxRatio = -1;
	
	float MinMaxStruct::minComponent = 10000;
	float MinMaxStruct::maxComponent = -1;

	float MinMaxStruct::minDegree = 0;
	float MinMaxStruct::maxDegree = -1;

	std::vector<std::vector<float>> ColorTable::colors = {
		{0.902f, 0.098f, 0.294f, 1.0f},
			{0.235f, 0.706f, 0.294f, 1.0f},
			{1.0f, 0.882f, 0.098f, 1.0f},
			{0.0f, 0.51f, 0.784f, 1.0f},
			{0.961f, 0.51f, 0.188f, 1.0f},
			{0.569f, 0.118f, 0.706f, 1.0f},
			{0.275f, 0.941f, 0.941f, 1.0f},
			{0.941f, 0.196f, 0.902f, 1.0f},
			{0.824f, 0.961f, 0.235f, 1.0f},
			{0.98f, 0.745f, 0.745f, 1.0f},
			{0.0f, 0.502f, 0.502f, 1.0f},
			{0.902f, 0.745f, 1.0f, 1.0f},
			{0.667f, 0.431f, 0.157f, 1.0f},
			{1.0f, 0.98f, 0.784f, 1.0f},
			{0.502f, 0.0f, 0.0f, 1.0f},
			{0.667f, 1.0f, 0.765f, 1.0f},
			{0.502f, 0.502f, 0.0f, 1.0f},
			{1.0f, 0.843f, 0.706f, 1.0f},
			{0.0f, 0.0f, 0.502f, 1.0f},
			{0.502f, 0.502f, 0.502f, 1.0f},
			{0.902f, 0.098f, 0.294f, 1.0f},
			{0.235f, 0.706f, 0.294f, 1.0f},
			{1.0f, 0.882f, 0.098f, 1.0f},
			{0.0f, 0.51f, 0.784f, 1.0f},
			{0.961f, 0.51f, 0.188f, 1.0f},
			{0.569f, 0.118f, 0.706f, 1.0f},
			{0.275f, 0.941f, 0.941f, 1.0f},
			{0.941f, 0.196f, 0.902f, 1.0f},
			{0.824f, 0.961f, 0.235f, 1.0f},
			{0.98f, 0.745f, 0.745f, 1.0f},
			{0.0f, 0.502f, 0.502f, 1.0f},
			{0.902f, 0.745f, 1.0f, 1.0f},
			{0.667f, 0.431f, 0.157f, 1.0f},
			{1.0f, 0.98f, 0.784f, 1.0f},
			{0.502f, 0.0f, 0.0f, 1.0f},
			{0.667f, 1.0f, 0.765f, 1.0f},
			{0.502f, 0.502f, 0.0f, 1.0f},
			{1.0f, 0.843f, 0.706f, 1.0f},
			{0.0f, 0.0f, 0.502f, 1.0f},
			{0.502f, 0.502f, 0.502f, 1.0f}
	};

	/////////////////////////////////////////////// BMetaNode /////////////////////////////////////
	BMetaNode::BMetaNode()
	{
		mSrcVert = 0;
		mSrcSkeleton = nullptr;
		connectedComponent = -1;
		nodeThickness = 0.0;
		nodeWidth = 0.0;
		hasGeom = false;
		p[0] = 0;
		p[1] = 0;
		p[2] = 0;
		degree = 0;
		for (int i = 0; i < 4; ++i)
		{
			glThicknessColor[i] = defaultColor[i];
			glWidthColor[i] = defaultColor[i];
			glDegreeColor[i] = defaultColor[i];
			glComponentColor[i] = defaultColor[i];
			currentColor = glThicknessColor;
		}
	}
	BMetaNode::BMetaNode(SkelVert srcId, BSkeleton *skel)
	{
		mSrcVert = srcId;
		mSrcSkeleton = skel;
		//std::cout << "number of vertices in pointed skeleton " << boost::num_vertices(*skel) << std::endl;
		connectedComponent = -1;
		nodeThickness = 0.0;
		nodeWidth = 0.0;
		hasGeom = true;
		p[0] = skel->operator[](srcId)[0];
		p[1] = skel->operator[](srcId)[1];
		p[2] = skel->operator[](srcId)[2];
		float degree = boost::degree(mSrcVert, *skel);
		MinMaxStruct::minDegree = std::min(degree, MinMaxStruct::minDegree);
		MinMaxStruct::maxDegree = std::max(degree, MinMaxStruct::maxDegree);
		for (int i = 0; i < 4; ++i)
		{
			glThicknessColor[i] = defaultColor[i];
			glWidthColor[i] = defaultColor[i];
			glDegreeColor[i] = defaultColor[i];
			glComponentColor[i] = defaultColor[i];
			currentColor = glThicknessColor;
		}
	}

	void BMetaNode::updateColors(NodeVisualizationOptions options)
	{
		int length = options.heatmap.size() - 1;

		//std::cout << "Update Colors node thickness " << nodeThickness << " node width " << nodeWidth << std::endl;

		float thickRatio = ((nodeThickness - MinMaxStruct::minThickness) / (MinMaxStruct::maxThickness - MinMaxStruct::minThickness));
		float widthRatio = ((nodeWidth - MinMaxStruct::minWidth) / (MinMaxStruct::maxWidth - MinMaxStruct::minWidth));
		float degreeRatio = ((1.0*degree - MinMaxStruct::minDegree) / (MinMaxStruct::maxDegree - MinMaxStruct::minDegree));

		thickRatio = (thickRatio - options.minColorCutoff) / (options.maxColorCutoff - options.minColorCutoff);
		widthRatio = (widthRatio - options.minColorCutoff) / (options.maxColorCutoff - options.minColorCutoff);
		degreeRatio = (degreeRatio - options.minColorCutoff) / (options.maxColorCutoff - options.minColorCutoff);

		thickRatio = std::min(thickRatio, 1.0f);
		thickRatio = std::max(thickRatio, 0.0f);
		widthRatio = std::min(widthRatio, 1.0f);
		widthRatio = std::max(widthRatio, 0.0f);
		degreeRatio = std::min(degreeRatio, 1.0f);
		degreeRatio = std::max(degreeRatio, 0.0f);

		int thickpos = thickRatio * length;
		int widthpos = widthRatio * length;
		int degreepos = degreeRatio * length;

		thickpos = std::min(thickpos, length);
		widthpos = std::min(widthpos, length);
		degreepos = std::min(degreepos, length);

		thickpos = std::max(thickpos, 0);
		widthpos = std::max(widthpos, 0);
		degreepos = std::max(degreepos, 0);



		

		if (thickpos > length || widthpos > length || degreepos > length)
		{
			std::cout << "Value out of range : thickpos " << thickpos << " width pos : " << widthpos << " ratiopos : " << degreepos << std::endl;
		}
		int p = 3;

		

		std::vector<GLfloat> thicknessColor = options.heatmap[thickpos];
		std::vector<GLfloat> widthColor = options.heatmap[widthpos];
		std::vector<GLfloat> degreeColor = options.heatmap[degreepos];
		for (int c = 0; c < 3; ++c)
		{
			glThicknessColor[c] = thicknessColor[c];

			glWidthColor[c] = widthColor[c];
			
			glDegreeColor[c] = degreeColor[c];
		}
		currentColor = glThicknessColor;
		updateComponentColor();
	}

	void BMetaNode::updateComponentColor()
	{

		std::vector<float> componentColor = ColorTable::getComponentColor(connectedComponent);

		for (int c = 0; c < 3; ++c)
		{
			glComponentColor[c] = componentColor[c];
		}
	}

	////////////////////////////////////////////// BMetaEdge ///////////////////////////////////////

	BMetaEdge::BMetaEdge()
	{
		mVertices = {};
		mEdges = {};
		indicesList = {};
		localThicknessColors = {};
		localWidthColors = {};
		localRatioColors = {};
		localComponentColors = {};
		averageThickness = 0.0;
		averageWidth = 0.0;
		mLength = 0.0;
		isBridge = false;
		isSelected = false;
		connectedComponent = -1;
		//mSrcSkeleton = nullptr;

	}

	BMetaEdge::BMetaEdge(std::vector<SkelVert> vertices, BSkeleton* srcSkeleton, std::vector<GLfloat> &glVertices)
	{
		mVertices = vertices;
		mEdges = {};
		indicesList = {};
		localThicknessColors = {};
		localWidthColors = {};
		localRatioColors = {};
		localComponentColors = {};
		float weightedAvgThickness = 0.0;
		float weightedAvgWidth = 0.0;
		float skelEdgeLength = 0.0;
		averageThickness = 0.0;
		averageWidth = 0.0;
		mLength = 0.0;
		isBridge = false;
		isSelected = false;
		connectedComponent = -1;
		for (int i = 0; i < mVertices.size() - 1; ++i)
		{
			SkelVert v0 = mVertices[i];
			SkelVert v1 = mVertices[i + 1];
			SkelEdge e;
			bool exists;
			boost::tie(e, exists) = boost::edge(v0, v1, *srcSkeleton);
			if (exists)
			{
				RootAttributes skelRootAttributes = srcSkeleton->operator[](e);

				float thickness = skelRootAttributes.thickness;
				float width = skelRootAttributes.width;
				float ratio = thickness / width;

				MinMaxStruct::minThickness = std::min(thickness, MinMaxStruct::minThickness);
				MinMaxStruct::maxThickness = std::max(thickness, MinMaxStruct::maxThickness);

				MinMaxStruct::minWidth = std::min(width, MinMaxStruct::minWidth);
				MinMaxStruct::maxWidth = std::max(width, MinMaxStruct::maxWidth);

				MinMaxStruct::minRatio = std::min(ratio, MinMaxStruct::minRatio);
				MinMaxStruct::maxRatio = std::max(ratio, MinMaxStruct::maxRatio);

				mEdges.push_back(skelRootAttributes);
				float length = skelRootAttributes.euclidLength;
				weightedAvgThickness += thickness * length;
				weightedAvgWidth += width * length;
				skelEdgeLength += length;
				indicesList.push_back(skelRootAttributes.glV0);
				indicesList.push_back(skelRootAttributes.glV1);
			}
			else
			{
				std::cout << "Attempting to create a metaedge over skeleton vertices that are not already an edge" << std::endl;
				exit(1);
				//this is okay - this just means that this edge does not exist at the skeleton
				//level, eg it is an added edge at the metagraph level
				//do nothing here
				/*float cumThickness = 0.0;
				float cumWidth = 0.0;
				int numEdges = 0;
				std::cout << "Not seeing a skeleton edge " << std::endl;

				BSkeleton::adjacency_iterator adjIt, adjEnd;
				boost::tie(adjIt, adjEnd) = boost::adjacent_vertices(v0, *srcSkeleton);
				for (; adjIt != adjEnd; ++adjIt)
				{
					SkelVert leadVert = *adjIt;
					boost::tie(e, exists) = boost::edge(v0, leadVert, *srcSkeleton);
					RootAttributes skelRootAttributes = srcSkeleton->operator[](e);
					cumThickness += skelRootAttributes[Thickness];
					cumWidth += skelRootAttributes[Width];
					++numEdges;
				}
				boost::tie(adjIt, adjEnd) = boost::adjacent_vertices(v1, *srcSkeleton);
				for (; adjIt != adjEnd; ++adjIt)
				{
					SkelVert leadVert = *adjIt;
					boost::tie(e, exists) = boost::edge(v1, leadVert, *srcSkeleton);
					RootAttributes skelRootAttributes = srcSkeleton->operator[](e);
					cumThickness += skelRootAttributes[Thickness];
					cumWidth += skelRootAttributes[Width];
					++numEdges;
				}
				if (numEdges == 0)
				{
					numEdges = 1;
				}
				weightedAvgThickness += cumThickness / numEdges;
				weightedAvgThickness += cumWidth / numEdges;
				Point3d dif = srcSkeleton->operator[](v0) - srcSkeleton->operator[](v1);
				skelEdgeLength = dif.mag();*/
			}
		}
		if (skelEdgeLength == 0)
		{
			skelEdgeLength = 1;
		}
		averageThickness = weightedAvgThickness / skelEdgeLength;
		
		if (isnan(averageThickness))
		{
			averageThickness = 1.0;
			std::cout << "Note : improper MetaEdge thickness, skeleton data may be missing thickness data on one or more roots" << std::endl;
		}

		averageWidth = weightedAvgWidth / skelEdgeLength;

		if (isnan(averageWidth))
		{
			averageWidth = 1.0;
			std::cout << "Note : improper MetaEdge width, skeleton data may be missing width data on one or more roots" << std::endl;
		}

		mLength = skelEdgeLength;

		for (int i = 0; i < mVertices.size() - 1; ++i)
		{
			for (int goTwice = 0; goTwice < 2; ++goTwice)
			{
				for (int c = 0; c < 3; ++c)
				{
					localThicknessColors.push_back(defaultColor[c]);
					localWidthColors.push_back(defaultColor[c]);
					localRatioColors.push_back(defaultColor[c]);
					localComponentColors.push_back(defaultColor[c]);
				}
			}
		}

		//if (addToVertices)
		//{
		//	indicesList = {};
		//	int vertexIndex = edgeVertices.size() / 3;
		//	/*glVertices = edgeVertices.end()._Ptr;
		//	glThicknessColors = edgeColors[0].data() + edgeColors[0].size();
		//	glWidthColors = edgeColors[1].data() + edgeColors[1].size();
		//	glRatioColors = edgeColors[2].data() + edgeColors[2].size();
		//	glComponentColors = edgeColors[3].data() + edgeColors[3].size();;*/
		//	for (int i = 0; i < mVertices.size() - 1; ++i)
		//	{
		//		Point3d p0 = srcSkeleton->operator[](mVertices[i]);
		//		edgeVertices.push_back(p0.x());
		//		edgeVertices.push_back(p0.y());
		//		edgeVertices.push_back(p0.z());
		//		for (int coloringType = 0; coloringType < 4; ++coloringType)
		//		{
		//			for (int c = 0; c < 3; ++c)
		//			{
		//				edgeColors[coloringType].push_back(defaultColor[c]);
		//			}
		//		}
		//		
		//		indicesList.push_back(vertexIndex);
		//		++vertexIndex;

		//		Point3d p1 = srcSkeleton->operator[](mVertices[i + 1]);
		//		edgeVertices.push_back(p1.x());
		//		edgeVertices.push_back(p1.y());
		//		edgeVertices.push_back(p1.z());
		//		for (int coloringType = 0; coloringType < 4; ++coloringType)
		//		{
		//			for (int c = 0; c < 3; ++c)
		//			{
		//				edgeColors[coloringType].push_back(defaultColor[c]);
		//			}
		//		}
		//		indicesList.push_back(vertexIndex);
		//		++vertexIndex;
		//	}
		//}

		//if (addToIndices)
		//{
		//	addToIndicesList(edgeIndices);
		//}

		//mSrcSkeleton = srcSkeleton;
	}

	void BMetaEdge::addToIndicesList(std::vector<GLuint> &edgeIndices)
	{
		edgeIndices.insert(edgeIndices.end(), indicesList.begin(), indicesList.end());
	}

	BMetaEdge BMetaEdge::join(BMetaEdge &other, BSkeleton *srcSkeleton)
	{
		int myMaxId = mVertices.size() - 1;
		int otherMaxId = other.mVertices.size() - 1;
		std::vector<SkelVert> joinedVertices = {};
		if (mVertices[myMaxId] != other.mVertices[0] &&
			mVertices[myMaxId] != other.mVertices[otherMaxId] &&
			mVertices[0] != other.mVertices[0] &&
			mVertices[0] != other.mVertices[otherMaxId])
		{
			//std::cout << "Meta edges are not successfully joined.  No endpoints match between either edge." << std::endl;
			return BMetaEdge();
		}
		//in the case that the last vertex in this edge is the starting vertex in the 
		//other edge, simply append all but the first vertex of that edge to the vertex list
		//of this
		if (mVertices[myMaxId] == other.mVertices[0])
		{
			for (int i = 0; i < mVertices.size(); ++i)
			{
				joinedVertices.push_back(mVertices[i]);
			}
			for (int i = 1; i < other.mVertices.size(); ++i)
			{
				joinedVertices.push_back(other.mVertices[i]);
			}
		}
		else if (mVertices[myMaxId] == other.mVertices[otherMaxId])
		{
			for (int i = 0; i < mVertices.size(); ++i)
			{
				joinedVertices.push_back(mVertices[i]);
			}
			for (int i = otherMaxId - 1; i >= 0; --i)
			{
				joinedVertices.push_back(other.mVertices[i]);
			}
		}
		else if (mVertices[0] == other.mVertices[0])
		{
			for (int i = myMaxId; i >= 0; --i)
			{
				joinedVertices.push_back(mVertices[i]);
			}
			for (int i = 1; i <= otherMaxId; ++i)
			{
				joinedVertices.push_back(other.mVertices[i]);
			}
		}
		else
		{
			//to avoid repetitive code, simply flip the of the ordering for the same result
			//if matches are not made on this end.
			return other.join(*this, srcSkeleton);
		}

		std::vector<GLfloat> fakeVerts = {};
		BMetaEdge result = BMetaEdge(joinedVertices, srcSkeleton, fakeVerts);
		result.connectedComponent = other.connectedComponent;
		return result;

	}

	std::pair<BMetaEdge, BMetaEdge> BMetaEdge::split(SkelVert toSplitOn, BSkeleton *srcSkeleton)
	{
		std::vector<SkelVert> leftEdge = {}, rightEdge = {};
		std::vector<GLuint> leftEdgeIndices = {}, rightEdgeIndices = {};
		bool onLeft = true;
		for (int i = 0; i < mVertices.size(); ++i)
		{
			
			if (mVertices[i] == toSplitOn)
			{
				onLeft = false;
				rightEdge.push_back(mVertices[i]);
				leftEdge.push_back(mVertices[i]);
				continue;
			}
			if (onLeft)
			{
				leftEdge.push_back(mVertices[i]);
				//leftEdgeIndices.push_back(indicesList[i * 2]);
				//leftEdgeIndices.push_back(indicesList[i * 2 + 1]);
			}
			if (!onLeft)
			{
				rightEdge.push_back(mVertices[i]);
				//rightEdgeIndices.push_back(indicesList[i * 2 - 2]);
				//rightEdgeIndices.push_back(indicesList[i * 2 - 1]);
			}
		}
		std::cout << std::endl << "Left edge verts";
		for (SkelVert v : leftEdge)
		{
			std::cout << " " <<v;
		}
		std::cout << std::endl << "left edge indices";
		for (GLuint idx : leftEdgeIndices)
		{
			std::cout << " " <<idx;
		}
		std::cout << std::endl << "right edge verts";
		for (SkelVert v : rightEdge)
		{
			std::cout << " " << v;
		}
		std::cout << std::endl << "right edge indices";
		for (GLuint idx : rightEdgeIndices)
		{
			std::cout << " " << idx;
		}
		

		std::vector<GLfloat> fakeverts = {};
		BMetaEdge left = BMetaEdge(leftEdge, srcSkeleton, fakeverts);
		BMetaEdge right = BMetaEdge(rightEdge, srcSkeleton, fakeverts);

		return std::make_pair(left, right);
	}

	float BMetaEdge::getAvgThickness()
	{
		return averageThickness;
	}

	SkelVert BMetaEdge::start()
	{
		return mVertices[0];
	}
	SkelVert BMetaEdge::end()
	{
		return mVertices[mVertices.size() - 1];
	}

	void BMetaEdge::updateColors(EdgeVisualizationOptions options, std::vector<std::vector<GLfloat>> &edgeColors)
	{
		int length = options.heatmap.size() - 1;

		localThicknessColors.resize(mEdges.size() * 2 * 3);
		localWidthColors.resize(mEdges.size() * 2 * 3);
		localRatioColors.resize(mEdges.size() * 2 * 3);
		localComponentColors.resize(mEdges.size() * 2 * 3);

		for (int i = 0; i < mEdges.size(); ++i)
		{
			RootAttributes att = mEdges[i];

			float thickRatio = ((att.thickness - MinMaxStruct::minThickness) / (MinMaxStruct::maxThickness - MinMaxStruct::minThickness));
			float widthRatio = ((att.width - MinMaxStruct::minWidth) / (MinMaxStruct::maxWidth - MinMaxStruct::minWidth));
			float ratio = att.thickness / att.width;
			float ratioRatio = ((ratio - MinMaxStruct::minRatio) / (MinMaxStruct::maxRatio - MinMaxStruct::minRatio));


			thickRatio = (thickRatio - options.minColorCutoff) / (options.maxColorCutoff - options.minColorCutoff);
			widthRatio = (widthRatio - options.minColorCutoff) / (options.maxColorCutoff - options.minColorCutoff);
			ratioRatio = 1.0 - (ratioRatio - options.minColorCutoff) / (options.maxColorCutoff - options.minColorCutoff);

			thickRatio = std::min(1.0f, thickRatio);
			thickRatio = std::max(0.0f, thickRatio);

			widthRatio = std::min(1.0f, widthRatio);
			widthRatio = std::max(0.0f, widthRatio);

			ratioRatio = std::min(1.0f, ratioRatio);
			ratioRatio = std::max(0.0f, ratioRatio);

			
			int thickpos = thickRatio * length;
			int widthpos = widthRatio * length;
			int ratiopos = ratioRatio * length;

			if (thickpos > length || widthpos > length || ratiopos > length)
			{
				std::cout << "Value out of range : thickpos " << thickpos << " width pos : " << widthpos << " ratiopos : " << ratiopos << std::endl;
			}

			thickpos = std::min(thickpos, length);
			widthpos = std::min(widthpos, length);
			ratiopos = std::min(ratiopos, length);
			thickpos = std::max(thickpos, 0);
			widthpos = std::max(widthpos, 0);
			ratiopos = std::max(ratiopos, 0);

			std::vector<GLfloat> thicknessColor = options.heatmap[thickpos];
			std::vector<GLfloat> widthColor = options.heatmap[widthpos];
			std::vector<GLfloat> ratioColor = options.heatmap[ratiopos];
			
			for (int c = 0; c < 3; ++c)
			{
				/*glThicknessColors[i * 6 + c] = boost::python::extract<float>(thicknessColor[c]);
				glThicknessColors[i * 6 + c + 3] = boost::python::extract<float>(thicknessColor[c]);*/

				localThicknessColors[i * 6 + c] = thicknessColor[c];
				localThicknessColors[i * 6 + c + 3] = thicknessColor[c];

				/*glWidthColors[i * 6 + c] = boost::python::extract<float>(widthColor[c]);
				glWidthColors[i * 6 + c + 3] = boost::python::extract<float>(widthColor[c]);*/

				localWidthColors[i * 6 + c] = widthColor[c];
				localWidthColors[i * 6 + c + 3] = widthColor[c];

				/*glRatioColors[i * 6 + c] = boost::python::extract<float>(ratioColor[c]);
				glRatioColors[i * 6 + c + 3] = boost::python::extract<float>(ratioColor[c]);*/

				localRatioColors[i * 6 + c] = ratioColor[c];
				localRatioColors[i * 6 + c + 3] = ratioColor[c];

			}

		}

		updateComponentColor(edgeColors);

		updateGraphColors(edgeColors);
	}

	void BMetaEdge::updateComponentColor(std::vector<std::vector<GLfloat>> &edgeColors)
	{

		std::vector<float> componentColor = ColorTable::getComponentColor(connectedComponent);

		for (int i = 0; i < mEdges.size(); ++i)
		{
			for (int c = 0; c < 3; ++c)
			{

				localComponentColors[i * 6 + c] = componentColor[c];
				localComponentColors[i * 6 + c + 3] = componentColor[c];
			}
		}
		updateGraphColors(edgeColors);
	}

	void BMetaEdge::select(GLfloat *selectionColor, std::vector<std::vector<GLfloat>> &edgeColors)
	{
		isSelected = true;
		std::cout << "Calling select" << std::endl;
		std::cout << "edge count " << mEdges.size() << std::endl;
		for (int i = 0; i < indicesList.size(); ++i)
		{
			for (int coloringType = 0; coloringType < 4; ++coloringType)
			{
				for (int c = 0; c < 3; ++c)
				{
					edgeColors[coloringType][indicesList[i] * 3 + c] = selectionColor[c];
				}
			}
		}
		std::cout << "ending select" << std::endl;
	}

	void BMetaEdge::unselect(std::vector<std::vector<GLfloat>> &edgeColors)
	{
		isSelected = false;
		updateGraphColors(edgeColors);
	}

	void BMetaEdge::updateGraphColors(std::vector<std::vector<GLfloat>> &edgeColors)
	{
		int maxIndex = 0;
		for (int i = 0; i < indicesList.size(); ++i)
		{
			maxIndex = std::max(maxIndex, (int)indicesList[i]);
		}
		if (maxIndex * 3 + 3 > edgeColors[0].size())
		{
			for (int i = 0; i < 4; ++i)
			{
				edgeColors[i].resize(maxIndex * 3 + 3);
			}
		}
		for (int i = 0; i < indicesList.size(); ++i)
		{
			for (int c = 0; c < 3; ++c)
			{
				edgeColors[0][indicesList[i] * 3 + c] = localThicknessColors[i * 3 + c];
				edgeColors[1][indicesList[i] * 3 + c] = localWidthColors[i * 3 + c];
				edgeColors[2][indicesList[i] * 3 + c] = localRatioColors[i * 3 + c];
				edgeColors[3][indicesList[i] * 3 + c] = localComponentColors[i * 3 + c];
			}
		}
	}


	BoundingBox::BoundingBox()
	{
		minx = 0;
		maxx = 0;
		miny = 0;
		maxy = 0;
		minz = 0;
		maxz = 0;

		hasPoints = false;


		for (int i = 0; i < 8; ++i)
		{
			corners.push_back({ 0, 0, 0 });
		}
	}

	void BoundingBox::addPoint(float *p)
	{
		if (!hasPoints)
		{
			hasPoints = true;
			minx = p[0];
			maxx = p[0];
			miny = p[1];
			maxy = p[1];
			minz = p[2];
			maxz = p[2];
		}
		else
		{
			minx = std::min(minx, p[0]);
			maxx = std::max(maxx, p[0]);
			miny = std::min(miny, p[1]);
			maxy = std::max(maxy, p[1]);
			minz = std::min(minz, p[2]);
			maxz = std::max(maxz, p[2]);
		}

		corners[0][0] = minx;
		corners[1][0] = minx;
		corners[2][0] = minx;
		corners[3][0] = minx;
		corners[4][0] = maxx;
		corners[5][0] = maxx;
		corners[6][0] = maxx;
		corners[7][0] = maxx;

		corners[0][1] = miny;
		corners[1][1] = miny;
		corners[2][1] = maxy;
		corners[3][1] = maxy;
		corners[4][1] = miny;
		corners[5][1] = miny;
		corners[6][1] = maxy;
		corners[7][1] = maxy;

		corners[0][2] = minz;
		corners[1][2] = maxz;
		corners[2][2] = maxz;
		corners[3][2] = minz;
		corners[4][2] = minz;
		corners[5][2] = maxz;
		corners[6][2] = maxz;
		corners[7][2] = minz;
	}


	void BoundingBox::draw(std::vector<GLfloat> componentColor, float lineWidth)
	{
		glColor3f(componentColor[0], componentColor[1], componentColor[2]);
		glLineWidth(lineWidth);
		glBegin(GL_LINES);

		for (int i = 0; i < 3; ++i)
		{
			glVertex3f(corners[i][0], corners[i][1], corners[i][2]);
			glVertex3f(corners[i + 1][0], corners[i + 1][1], corners[i + 1][2]);

			glVertex3f(corners[i+4][0], corners[i+4][1], corners[i+4][2]);
			glVertex3f(corners[i + 5][0], corners[i + 5][1], corners[i + 5][2]);
		}

		int i = 3;

		glVertex3f(corners[i][0], corners[i][1], corners[i][2]);
		glVertex3f(corners[0][0], corners[0][1], corners[0][2]);

		glVertex3f(corners[i + 4][0], corners[i + 4][1], corners[i + 4][2]);
		glVertex3f(corners[i+1][0], corners[i+1][1], corners[i+1][2]);

		for (int i = 0; i < 4; ++i)
		{
			glVertex3f(corners[i][0], corners[i][1], corners[i][2]);
			glVertex3f(corners[i + 4][0], corners[i + 4][1], corners[i + 4][2]);
		}


		glEnd();
	}
	///////////////////////////////////////// MetaNode3d //////////////////////////////////////////

	//MetaNode3d::MetaNode3d()
	//{
	//	order = -1;
	//	connectedComponent = -1;
	//	nodeThickness = -1;
	//	nodeWidth = -1;
	//	degree = 0;
	//}

	//MetaNode3d::MetaNode3d(BMetaNode src, Roots::BSkeleton *srcSkeleton, int aOrder, int aDegree)
	//	:Point3d(srcSkeleton->operator[](src.mSrcVert))
	//{
	//	if (src.hasGeom)
	//	{
	//		p[0] = src[0];
	//		p[1] = src[1];
	//		p[2] = src[2];
	//	}
	//	order = aOrder;
	//	connectedComponent = src.connectedComponent;
	//	nodeThickness = src.nodeThickness;
	//	nodeWidth = src.nodeWidth;
	//	degree = aDegree;
	//}



	/////////////////////////////////////////// MetaEdge3d //////////////////////////////////////////

	//MetaEdge3d::MetaEdge3d()
	//{
	//	node0 = 0;
	//	node1 = 0;
	//	avgThickness = 0;
	//	avgWidth = 0;
	//	connectedComponent = -1;
	//	isBridge = false;
	//	//edges = boost::python::list();

	//}

	//MetaEdge3d::MetaEdge3d(int n0, int n1, float thickness, float width, int component, int aOrder, std::vector<RootAttributes> aEdges, bool aIsBridge)
	//{
	//	node0 = n0;
	//	node1 = n1;
	//	avgThickness = thickness;
	//	avgWidth = width;
	//	connectedComponent = component;
	//	order = aOrder;
	//	edges = toPythonList<RootAttributes>(aEdges);
	//	isBridge = aIsBridge;
	//}


	///////////////////////////////////////// BMetaGraph //////////////////////////////////////////

	BMetaGraph::BMetaGraph()
	{
		Initialize();
	}

	void BMetaGraph::Initialize()
	{
		mSkeleton = BSkeleton();
		vertNodeMap = std::map<SkelVert, MetaV>();
		mComponentMap = {};
		mComponentSizeMap = {};
		mMinimumSpanningTree = {};
		numComponents = -1;
		mDuplicateEdgeMap = std::map<MetaE, MetaE>();
		mDuplicateNodeMap = std::map<MetaV, MetaV>();
		componentBounds = {};
		isLoaded = false;

		edgeOptions = EdgeVisualizationOptions();
		nodeOptions = NodeVisualizationOptions();
		mode = OperationMode::None;

		
		edgeVertices = {};
		edgeColors = { {},{},{},{} };
		//nonBridgeEdgeIndices = {};
		//bridgeEdgeIndices = {};
		//baseEdgeIndices = {};
		//componentIndices = { {} };
		for (int i = 0; i < 3; ++i)
		{
			selectionColor[i] = 0.0;
		}
		selectionColor[3] = 1.0;
		selectionColor[2] = 1.0;

		onlyDisplaySelectedComponents = false;
		selectedComponent1 = 0, selectedComponent2 = 0;
		showBoundingBoxes = false;

		selectNode1Valid = false;
		selectNode2Valid = false;
		selectNode1 = 99999999;
		selectNode2 = 99999999;

		breakEdgeValid = false;

		splitEdgeValid = false;
		splitNeighbors = {};

		pythonUpToDate = false;
		clear();
		arcball_reset();
	}

	//BMetaGraph::BMetaGraph(std::string filename)
	//{
	//	std::cout << "Creating new metagraph" << std::endl;;
	//	loadFromFile(filename);
	//	//std::cout << "Finished loading file" << std::endl;
	//}


	void BMetaGraph::loadFromFile(std::string filename)
	{
		Initialize();

		if (filename.length() == 0)
		{
			std::cout << "Provided filename is empty" << std::endl;
			return;
		}

		

		std::ifstream filestream;
		filestream.open(filename);
		std::vector<std::string> lines = {};
		std::string line;
		while (!filestream.eof())
		{
			std::getline(filestream, line);
			lines.push_back(line);
		}
		filestream.close();
		//for each(std::string line in lines)
		//{
		//	std::cout << line << std::endl;
		//}
		loadFromLines(lines, 0);
		findAndLabelConnectedComponents();
		isLoaded = true;
	}

	//void BMetaGraph::loadFromFile(char *filename)
	//{
	//	std::string wrappedFile = std::string(filename);
	//	loadFromFileString(wrappedFile);
	//}

	int BMetaGraph::loadFromLines(std::vector<std::string> lines, int startingLine)
	{
		int result = loadSkeletonFromLines(lines, startingLine);
		//std::cout << "sucessfully loaded skeleton from file" << std::endl;
		pythonUpToDate = false;
		//std::cout << "successfully updated python " << std::endl;
		mSkeleton.findBoundingSphere();
		if (useArcball)
		{
			mSkeleton = mSkeleton.recenterSkeleton(Point3d(), edgeVertices);
		}
		initializeFromSkeleton();
		findBridges();
		return result;
	}

	void BMetaGraph::writeToStream(std::ostream & out)
	{
		mSkeleton.writeToStream(out);
	}

	void BMetaGraph::writeToFile(std::string filename)
	{
		std::ofstream filestream;
		filestream.open(filename);
		writeToStream(filestream);
	}

	int BMetaGraph::loadSkeletonFromLines(std::vector<std::string> lines, int & startingLine)
	{
		return mSkeleton.loadFromLines(lines, startingLine, edgeVertices);
	}

	

	void BMetaGraph::initializeFromSkeleton()
	{
		isLoaded = true;
		std::vector<SkelVert> nodesOfInterest = mSkeleton.GetNotableVertices();
		std::vector<bool> visitedNodes = std::vector<bool>(mSkeleton.m_vertices.size(), false);

		for each(SkelVert node in nodesOfInterest)
		{
			MetaV startV = addNode(node, &mSkeleton);
			BSkeleton::adjacency_iterator adjIt, adjEnd;
			boost::tie(adjIt, adjEnd) = boost::adjacent_vertices(node, mSkeleton);

			for (; adjIt != adjEnd; ++adjIt)
			{
				SkelVert leadVert = *adjIt;
				buildMetaEdge(startV, leadVert, visitedNodes, true);
			}
		}
		
		findAndLabelConnectedComponents();


		pythonUpToDate = false;
	}

	MetaE BMetaGraph::buildMetaEdge(MetaV startV, SkelVert &lead, std::vector<bool> &visitedSkelVerts, bool isLoading)
	{
		//if the metaedge leading out of this node towards the lead has already been visited, then
		//we do not need to create a new edge
		if (visitedSkelVerts[lead])
		{
			return MetaE();
		}

		std::vector<SkelVert> edgeVerts = {};
		SkelVert startSkelVert = operator[](startV).mSrcVert;
		edgeVerts.push_back(startSkelVert);
		visitedSkelVerts[startSkelVert] = true;
		SkelVert last = startSkelVert;
		while (boost::degree(lead, mSkeleton) == 2 && lead != startSkelVert)
		{
			edgeVerts.push_back(lead);
			visitedSkelVerts[lead] = true;

			BSkeleton::adjacency_iterator adjIt, adjEnd;
			boost::tie(adjIt, adjEnd) = boost::adjacent_vertices(lead, mSkeleton);
			if (*adjIt != last)
			{
				last = lead;
				lead = *adjIt;
				continue;
			}
			++adjIt;
			if (*adjIt != last)
			{
				last = lead;
				lead = *adjIt;
				continue;
			}

			break;
		}

		edgeVerts.push_back(lead);
		MetaV endV = addNode(lead, &mSkeleton);
		return addEdge(startV, endV, edgeVerts, &mSkeleton, isLoading);
	}

	MetaV BMetaGraph::addNode(SkelVert srcVert, BSkeleton *srcSkel, bool inheretEdgeParams)
	{
		if (vertNodeMap.count(srcVert) == 0)
		{
			BMetaNode node = BMetaNode(srcVert, srcSkel);
			MetaV v = boost::add_vertex(*this);


			if (inheretEdgeParams)
			{
				BSkeleton::adjacency_iterator adjIt, adjEnd;
				boost::tie(adjIt, adjEnd) = boost::adjacent_vertices(srcVert, mSkeleton);

				double thickness = 0;
				double width = 0;
				int count = 0;

				for (; adjIt != adjEnd; ++adjIt)
				{
					SkelEdge e;
					bool exists = false;
					boost::tie(e, exists) = boost::edge(srcVert, *adjIt, mSkeleton);
					thickness += mSkeleton[e].thickness;
					width += mSkeleton[e].width;
					++count;
				}
				count = std::max(0, count);
				node.nodeThickness = thickness / count;
				node.nodeWidth = width / count;
				std::cout << "node thickenss " << node.nodeThickness << " node width " << node.nodeWidth << std::endl;

			}

			operator[](v) = node;
			vertNodeMap[srcVert] = v;
		}

		pythonUpToDate = false;
		updateNodeDegrees();
		


		return vertNodeMap[srcVert];

		
	}

	MetaE BMetaGraph::addEdge(MetaV v0, MetaV v1, std::vector<SkelVert> skelVerts, BSkeleton *srcSkel, bool isLoading, bool isJoining)
	{
		BMetaEdge edge;
		if (!isJoining)
		{
			edge = BMetaEdge(skelVerts, srcSkel, edgeVertices);
		}
		else
		{
			edge = BMetaEdge(skelVerts, srcSkel, edgeVertices);
		}
		MetaE e;
		bool success;
		boost::tie(e, success) = boost::add_edge(v0, v1, *this);


		operator[](v0).nodeThickness = std::max(operator[](v0).nodeThickness, edge.averageThickness);
		operator[](v1).nodeThickness = std::max(operator[](v1).nodeThickness, edge.averageThickness);

		operator[](v0).nodeWidth = std::max(operator[](v0).nodeWidth, edge.averageWidth);
		operator[](v1).nodeWidth = std::max(operator[](v1).nodeWidth, edge.averageWidth);
		
		if (success)
		{
			operator[](e) = edge;
		}

		if (!isLoading)
		{
			findBridges();
		}

		pythonUpToDate = false;
		updateNodeDegrees();
		return e;
	}

	int BMetaGraph::getNumGLEdges()
	{
		return 100;
	}
	int BMetaGraph::getNumGLVertices()
	{
		return 100;
	}



	void BMetaGraph::duplicateEdge(MetaE e0, BSkeleton *srcSkel, MetaV &dupV0, MetaV &dupV1, MetaE &dupE)
	{
		BMetaEdge edge = operator[](e0);

		std::vector<SkelVert> addedVerts = {};
		
		//iterate over source vertices, and add copies of their information to the skeleton
		if (edge.mVertices.size() > 0)
		{
			for each(SkelVert v in edge.mVertices)
			{
				addedVerts.push_back(srcSkel->addVertex(srcSkel->getVertData(v)));
			}

			bool exists = false;
			std::vector<SkelEdge> srcEdges = {};
			//find all the edge data associated with the original set of vertices
			for (int i = 0; i < edge.mVertices.size() - 1; ++i)
			{
				SkelEdge e = srcSkel->getEdge(edge.mVertices[i], edge.mVertices[i + 1], exists);
				if (exists)
				{
					srcEdges.push_back(e);
				}
			}

			//add edges between all of the new added vertices containing the info from the original edge data
			for (int i = 0; i < edge.mVertices.size() - 1; ++i)
			{
				srcSkel->addEdge(addedVerts[i], addedVerts[i + 1], srcSkel->getEdgeData(srcEdges[i]), edgeVertices);
			}
		}


		//add metanodes and edges for the augmented elements in the skeleton
		dupV0 = addNode(addedVerts[0], srcSkel);
		dupV1 = addNode(addedVerts[addedVerts.size() - 1], srcSkel);
		MetaV srcNode0 = vertNodeMap[edge.mVertices[0]];
		MetaV srcNode1 = vertNodeMap[edge.mVertices[edge.mVertices.size() - 1]];
		dupE = addEdge(dupV0, dupV1, addedVerts, srcSkel, false);

		
		mDuplicateEdgeMap[dupE] = e0;
		mDuplicateNodeMap[dupV0] = srcNode0;
		mDuplicateNodeMap[dupV1] = srcNode1;
		return;
	}

	void BMetaGraph::removeNode(MetaV nodeToRemove)
	{
		BMetaNode toRemove = operator[](nodeToRemove);
		vertNodeMap.erase(toRemove.mSrcVert);
		boost::remove_vertex(nodeToRemove, *this);
		updateVertNodeMap();
		updateNodeDegrees();
		pythonUpToDate = false;
	}

	void BMetaGraph::removeEdge(MetaE edgeToRemove)
	{
		BMetaEdge descriptor = operator[](edgeToRemove);

		/*MetaV v0 = boost::source(edgeToRemove, *this);
		MetaV v1 = boost::target(edgeToRemove, *this)*/;

		boost::remove_edge(edgeToRemove, *this);

		MetaV v0 = vertNodeMap[descriptor.end()];
		MetaV v1 = vertNodeMap[descriptor.start()];

		BMetaGraph::adjacency_iterator adjIt, adjEnd;

		BMetaGraph::out_edge_iterator edgeIt, edgeEnd;

		boost::tie(adjIt, adjEnd) = boost::adjacent_vertices(v0, *this);

		float maxThickness = -1;
		float maxWidth = -1;
		while (adjIt != adjEnd)
		{
			bool exists;
			MetaE adjEdge;
			boost::tie(adjEdge, exists) = boost::edge(v0, *adjIt, *this);

			maxThickness = std::max(maxThickness, operator[](adjEdge).averageThickness);
			maxWidth = std::max(maxWidth, operator[](adjEdge).averageWidth);
			++adjIt;
		}
		if (maxThickness <= 0)
		{
			maxThickness = 1;
		}
		if (maxWidth <= 0)
		{
			maxWidth = 1;
		}
		operator[](v0).nodeThickness = maxThickness;
		operator[](v0).nodeWidth = maxWidth;


		boost::tie(adjIt, adjEnd) = boost::adjacent_vertices(v1, *this);

		maxThickness = -1;
		maxWidth = -1;
		while (adjIt != adjEnd)
		{
			bool exists;
			MetaE adjEdge;
			boost::tie(adjEdge, exists) = boost::edge(v1, *adjIt, *this);

			maxThickness = std::max(maxThickness, operator[](adjEdge).averageThickness);
			maxWidth = std::max(maxWidth, operator[](adjEdge).averageWidth);
			++adjIt;
		}
		if (maxThickness <= 0)
		{
			maxThickness = 1;
		}
		if (maxWidth <= 0)
		{
			maxWidth = 1;
		}
		operator[](v1).nodeThickness = maxThickness / 2.0;

		
		//baseEdgeIndices = {};

		//metaEdgeIter mei = boost::edges(*this);

		//for (; mei.first != mei.second; ++mei)
		//{
		//	operator[](*mei.first).addToIndicesList(baseEdgeIndices);
		//}

		findBridges();

		//updateNodeDegrees();
		std::cout << "Edge removal successfull" << std::endl;
		pythonUpToDate = false;
	}

	void BMetaGraph::removeEdgeNoBridging(MetaE edgeToRemove)
	{
		MetaV v0 = boost::source(edgeToRemove, *this);
		MetaV v1 = boost::target(edgeToRemove, *this);

		boost::remove_edge(edgeToRemove, *this);


		BMetaGraph::adjacency_iterator adjIt, adjEnd;

		BMetaGraph::out_edge_iterator edgeIt, edgeEnd;

		boost::tie(adjIt, adjEnd) = boost::adjacent_vertices(v0, *this);

		float maxThickness = -1;
		float maxWidth = -1;
		while (adjIt != adjEnd)
		{
			bool exists;
			MetaE adjEdge;
			boost::tie(adjEdge, exists) = boost::edge(v0, *adjIt, *this);

			maxThickness = std::max(maxThickness, operator[](adjEdge).averageThickness);
			maxWidth = std::max(maxWidth, operator[](adjEdge).averageWidth);
			++adjIt;
		}
		if (maxThickness <= 0)
		{
			maxThickness = 1;
		}
		if (maxWidth <= 0)
		{
			maxWidth = 1;
		}
		operator[](v0).nodeThickness = maxThickness;
		operator[](v0).nodeWidth = maxWidth;


		boost::tie(adjIt, adjEnd) = boost::adjacent_vertices(v1, *this);

		maxThickness = -1;
		maxWidth = -1;
		while (adjIt != adjEnd)
		{
			bool exists;
			MetaE adjEdge;
			boost::tie(adjEdge, exists) = boost::edge(v1, *adjIt, *this);

			maxThickness = std::max(maxThickness, operator[](adjEdge).averageThickness);
			maxWidth = std::max(maxWidth, operator[](adjEdge).averageWidth);
			++adjIt;
		}
		if (maxThickness <= 0)
		{
			maxThickness = 1;
		}
		if (maxWidth <= 0)
		{
			maxWidth = 1;
		}
		operator[](v1).nodeThickness = maxThickness / 2.0;


		//baseEdgeIndices = {};

		//metaEdgeIter mei = boost::edges(*this);

		//for (; mei.first != mei.second; ++mei)
		//{
		//	operator[](*mei.first).addToIndicesList(baseEdgeIndices);
		//}

		updateNodeDegrees();
		pythonUpToDate = false;
	}

	void BMetaGraph::updateVertNodeMap()
	{
		vertNodeMap.clear();
		
		metaVertIter mvi = boost::vertices(*this);
		for (; mvi.first != mvi.second; ++mvi)
		{
			MetaV  node = *mvi.first;
			SkelVert vert = operator[](node).mSrcVert;
			vertNodeMap[vert] = node;
		}
	}

	void BMetaGraph::updateNodeDegrees()
	{
		metaVertIter mvi = boost::vertices(*this);
		for (; mvi.first != mvi.second; ++mvi)
		{
			MetaV  node = *mvi.first;
			int degree = boost::degree(node, *this);
			//operator[](node).degree = degree;

			MinMaxStruct::minDegree = std::min((float)degree, MinMaxStruct::minDegree);
			MinMaxStruct::maxDegree = std::max((float)degree, MinMaxStruct::maxDegree);

		}
		
	}


	void BMetaGraph::findAndLabelConnectedComponents()
	{
		std::cout << "Computing connected components" << std::endl;

		//build the boost component map vector<int> mapping vertices to components
		mComponentMap.resize(vertNodeMap.size(), 0);
		size_t numVertices = boost::num_vertices(*this);

		numComponents = boost::connected_components(*this, &mComponentMap[0]);
		mComponentSizeMap = {};
		mComponentSizeMap.resize(numComponents, 0);

		//build the component size map by determining for each edge the assigned component 
		//of its nodes, and adding its length to the cumulative size for each component
		metaEdgeIter mei = boost::edges(*this);
		for (; mei.first != mei.second; ++mei)
		{
			SkelVert v0 = operator[](*mei.first).start();
			MetaV node0 = vertNodeMap[v0];
			int component = mComponentMap[node0];
			mComponentSizeMap[component] += operator[](*mei.first).mLength;
		}

		//sort the components by their size (cumulative edge length)
		std::vector<float> allSizes = {};
		allSizes.insert(allSizes.begin(), mComponentSizeMap.begin(), mComponentSizeMap.end());
		

		std::sort(allSizes.begin(), allSizes.end());

		//map from current component Ids to prioritized component Ids - the highest length components
		//will have the highest priorty (iterate backwards on allsizes)
		std::map<int, int> componentPriorityMap = std::map<int, int>();
		int priority = 0;
		for (int i = allSizes.size() -1 ; i >= 0; --i)
		{
			for (int component = 0; component < mComponentSizeMap.size(); ++component)
			{
				if (mComponentSizeMap[component] == allSizes[i])
				{
					componentPriorityMap[component] = priority;
					++priority;
				}
			}
		}

		//remap the existing component map to the prioritized component map
		for (int i = 0; i < m_vertices.size(); ++i)
		{
			mComponentMap[i] = componentPriorityMap[mComponentMap[i]];
		}
		//update the component size map to match the sorted sizes of the priority map
		for (int i = 0, j = allSizes.size() - 1; i < allSizes.size(), j >= 0; ++i, --j)
		{
			mComponentSizeMap[i] = allSizes[j];
		}


		//remap the components of the nodes and edges locally
		for (MetaV node = 0; node < m_vertices.size(); ++node)
		{
			operator[](node).connectedComponent = mComponentMap[node];
			operator[](node).updateComponentColor();
		}

		//componentIndices = {};
		//componentIndices.resize(numComponents, {});
		mei = boost::edges(*this);
		for (; mei.first != mei.second; ++mei)
		{
			SkelVert v0 = operator[](*mei.first).start();
			MetaV node0 = vertNodeMap[v0];
			int component = mComponentMap[node0];
			operator[](*mei.first).connectedComponent = component;
			//operator[](*mei.first).addToIndicesList(componentIndices[component]);
			operator[](*mei.first).updateComponentColor(edgeColors);
		}
		MinMaxStruct::minComponent = 0;
		MinMaxStruct::maxComponent = numComponents - 1;

		componentBounds = std::vector<BoundingBox>(numComponents, BoundingBox());

		metaVertIter mvi = boost::vertices(*this);
		for (; mvi.first != mvi.second; ++mvi)
		{
			BMetaNode *node = &operator[](*mvi.first);
			componentBounds[node->connectedComponent].addPoint(node->p);
		}
	}

	boost::python::list BMetaGraph::getComponentSizes()
	{
		return toPythonList<float>(mComponentSizeMap);
	}

	void BMetaGraph::ConnectComponents(int compOne, int compTwo)
	{
		int lowerComponent = std::min(compOne, compTwo);
		int higherComponent = std::max(compOne, compTwo);
		metaVertIter mvi = boost::vertices(*this);

		for (; mvi.first != mvi.second; ++mvi)
		{
			int nodeComponent = operator[](*mvi.first).connectedComponent;
			if (nodeComponent == compOne || nodeComponent == compTwo)
			{
				operator[](*mvi.first).connectedComponent = lowerComponent;
			}
			if (nodeComponent > higherComponent)
			{
				operator[](*mvi.first).connectedComponent--;
			}
		}
	}

	void BMetaGraph::FixUpComponents()
	{
		metaVertIter mvi = boost::vertices(*this);

		for (; mvi.first != mvi.second; ++mvi)
		{
			int nodeComponent = operator[](*mvi.first).connectedComponent;
			BMetaGraph::adjacency_iterator adjIt, adjEnd;
			boost::tie(adjIt, adjEnd) = boost::adjacent_vertices(*mvi.first, *this);
			for (; adjIt != adjEnd; ++adjIt)
			{
				int neighborComponent = operator[](*adjIt).connectedComponent;
				if (neighborComponent != nodeComponent)
				{
					ConnectComponents(nodeComponent, neighborComponent);
				}
			}
		}
	}

	void BMetaGraph::FindMinimumSpanningTree()
	{
		mMinimumSpanningTree = {};

		kruskal_minimum_spanning_tree(*this, std::back_inserter(mMinimumSpanningTree),
			weight_map(get(&BMetaEdge::averageThickness, *this)));
	}

	int BMetaGraph::GetNumEdgesToBreak()
	{
		FindMinimumSpanningTree();
		return m_edges.size() - mMinimumSpanningTree.size();
	}

	//void BMetaGraph::JoinOperation(Point3d p1picked, Point3d p2picked)
	//{
	//	
	//	//check if both picked points are in the node map
	//	if (vertNodeMap.count(p1picked.id) == 0 ||
	//		vertNodeMap.count(p2picked.id) == 0)
	//	{
	//		//std::cout << "Neither point that has been picked is a MetaNode.  Fear ye who enter here" << std::endl;
	//		return;
	//	}
	//	
	//	BSkeleton::adjacency_iterator adjIt, adjEnd;
	//	boost::tie(adjIt, adjEnd) = boost::adjacent_vertices(p1picked.id, mSkeleton);
	//	float avgThick1 = 0.0, avgWidth1 = 0.0;
	//	int count1 = 0;
	//	for (; adjIt != adjEnd; ++adjIt)
	//	{
	//		SkelEdge e;
	//		bool exists;
	//		boost::tie(e, exists) = boost::edge(p1picked.id, *adjIt, mSkeleton);
	//		if (exists)
	//		{
	//			avgThick1 += mSkeleton[e].thickness;
	//			avgWidth1 += mSkeleton[e].width;
	//		}
	//		++count1;
	//	}
	//	avgThick1 /= count1;
	//	avgWidth1 /= count1;
	//	boost::tie(adjIt, adjEnd) = boost::adjacent_vertices(p2picked.id, mSkeleton);
	//	float avgThick2 = 0.0, avgWidth2 = 0.0;
	//	int count2 = 0;
	//	for (; adjIt != adjEnd; ++adjIt)
	//	{
	//		SkelEdge e;
	//		bool exists;
	//		boost::tie(e, exists) = boost::edge(p2picked.id, *adjIt, mSkeleton);
	//		if (exists)
	//		{
	//			avgThick2 += mSkeleton[e].thickness;
	//			avgWidth2 += mSkeleton[e].width;
	//		}
	//		++count2;
	//	}
	//	avgThick2 /= count2;
	//	avgWidth2 /= count2;

	//	RootAttributes newAttributes = RootAttributes((avgThick1 + avgThick2) / 2.0, (avgWidth1 + avgWidth2) / 2.0, p1picked, p2picked);



	//	SkelEdge newSkelEdge = mSkeleton.addEdge(p1picked.id, p1picked.id, newAttributes);

	//	//create an edge between these nodes and add it to the metagraph
	//	std::vector<SkelVert> bridgeVerts = {};
	//	bridgeVerts.push_back(p1picked.id);
	//	bridgeVerts.push_back(p2picked.id);
	//	MetaV node1 = vertNodeMap[p1picked.id];
	//	MetaV node2 = vertNodeMap[p2picked.id];
	//	MetaE newMetaEdge = addEdge(node1, node2, bridgeVerts, &mSkeleton, false);

	//	
	//	//deprecated
	//	//operator[](newMetaEdge).mEdges.push_back(newAttributes);
	//	//deprecated

	//	ConnectComponents(operator[](node1).connectedComponent, operator[](node2).connectedComponent);

	//	

	//	//join the components of the two nodes that have been picked
	//	//int componentOne = operator[](node1).connectedComponent;
	//	//int componentTwo = operator[](node2).connectedComponent;
	//	//if (componentOne != componentTwo)
	//	//{
	//	//	ConnectComponents(componentOne, componentTwo);
	//	//}


	//	//attempt to bridge both of the adjoining nodes (eg. if they are degree 2, then remove them
	//	//from the metagraph to be replaced by a single edge
	//	BridgeNode(node1);
	//	BridgeNode(node2);

	//	pythonUpToDate = false;
	//	return;
	//}


	void BMetaGraph::JoinOperation()
	{

		if (!selectNode1Valid || !selectNode2Valid)
		{
			return;
		}

		Point3d p1picked = mSkeleton[operator[](selectNode1).mSrcVert];
		Point3d p2picked = mSkeleton[operator[](selectNode2).mSrcVert];
			
		BSkeleton::adjacency_iterator adjIt, adjEnd;
		boost::tie(adjIt, adjEnd) = boost::adjacent_vertices(p1picked.id, mSkeleton);
		float avgThick1 = 0.0, avgWidth1 = 0.0;
		int count1 = 0;
		for (; adjIt != adjEnd; ++adjIt)
		{
			SkelEdge e;
			bool exists;
			boost::tie(e, exists) = boost::edge(p1picked.id, *adjIt, mSkeleton);
			if (exists)
			{
				avgThick1 += mSkeleton[e].thickness;
				avgWidth1 += mSkeleton[e].width;
			}
			++count1;
		}
		avgThick1 /= count1;
		avgWidth1 /= count1;
		boost::tie(adjIt, adjEnd) = boost::adjacent_vertices(p2picked.id, mSkeleton);
		float avgThick2 = 0.0, avgWidth2 = 0.0;
		int count2 = 0;
		for (; adjIt != adjEnd; ++adjIt)
		{
			SkelEdge e;
			bool exists;
			boost::tie(e, exists) = boost::edge(p2picked.id, *adjIt, mSkeleton);
			if (exists)
			{
				avgThick2 += mSkeleton[e].thickness;
				avgWidth2 += mSkeleton[e].width;
			}
			++count2;
		}
		avgThick2 /= count2;
		avgWidth2 /= count2;

		RootAttributes newAttributes = RootAttributes((avgThick1 + avgThick2) / 2.0, (avgWidth1 + avgWidth2) / 2.0, p1picked, p2picked);



		SkelEdge newSkelEdge = mSkeleton.addEdge(p1picked.id, p2picked.id, newAttributes, edgeVertices);

		//create an edge between these nodes and add it to the metagraph
		std::vector<SkelVert> bridgeVerts = {};
		bridgeVerts.push_back(operator[](selectNode1).mSrcVert);
		bridgeVerts.push_back(operator[](selectNode2).mSrcVert);
		MetaV node1 = vertNodeMap[p1picked.id];
		MetaV node2 = vertNodeMap[p2picked.id];

		if (boost::degree(node1, *this) == 0)
		{
			PromoteOperation(p1picked.id);
		}
		if (boost::degree(node2, *this) == 0)
		{
			PromoteOperation(p2picked.id);
		}


		MetaE newMetaEdge = addEdge(node1, node2, bridgeVerts, &mSkeleton, false, false);

		unselectAll();


		//deprecated
		//operator[](newMetaEdge).mEdges.push_back(newAttributes);
		//deprecated

		//ConnectComponents(operator[](node1).connectedComponent, operator[](node2).connectedComponent);

		findAndLabelConnectedComponents();
		operator[](newMetaEdge).updateColors(edgeOptions, edgeColors);
		
		if (node1 < node2)
		{
			BridgeNode(node2);
			BridgeNode(node1);
		}
		else
		{
			BridgeNode(node1);
			BridgeNode(node2);
		}

		//join the components of the two nodes that have been picked
		//int componentOne = operator[](node1).connectedComponent;
		//int componentTwo = operator[](node2).connectedComponent;
		//if (componentOne != componentTwo)
		//{
		//	ConnectComponents(componentOne, componentTwo);
		//}


		//attempt to bridge both of the adjoining nodes (eg. if they are degree 2, then remove them
		//from the metagraph to be replaced by a single edge
		/*BridgeNode(std::max(node1, node2));
		BridgeNode(std::max(node1, node2));*/
		//BridgeNodeSweeper();
		findAndLabelConnectedComponents();
		pythonUpToDate = false;
		std::cout << "Finished joining " << std::endl;
		return;
	}



	//void BMetaGraph::BreakOperation(MetaEdge3d toBreak) 
	//{
	//	MetaE e;
	//	bool exists;
	//	boost::tie(e, exists) = boost::edge(toBreak.node0, toBreak.node1, *this);
	//	
	//	if (exists)
	//	{
	//		std::cout << "Matching edge found to break" << std::endl;
	//		removeEdge(e);
	//	}
	//	else
	//	{
	//		std::cout << "No matching edge found to break" << std::endl;
	//	}
	//	BridgeNode(toBreak.node0);
	//	BridgeNode(toBreak.node1);

	//	pythonUpToDate = false;
	//	return;
	//}

	void BMetaGraph::BreakOperation()
	{
		std::cout << "C++ break operation " << std::endl;
		if (!breakEdgeValid)
		{
			return;
		}

		MetaV node0 = breakEdge.first;
		MetaV node1 = breakEdge.second;

		MetaE breakEdgeE;
		bool exists = false;

		boost::tie(breakEdgeE, exists) = boost::edge(node0, node1, *this);

		if (!exists)
		{
			std::cout << "Break edge doesnt exist " << std::endl;
			return;
		}

		removeEdge(breakEdgeE);
		breakEdgeValid = false;


		findAndLabelConnectedComponents();
		std::vector<MetaV> orderedNodes = { node0, node1 };
		std::sort(orderedNodes.begin(), orderedNodes.end());

		for (int i = orderedNodes.size() - 1; i >= 0; --i)
		{
			if (boost::degree(orderedNodes[i], *this) == 0)
			{
				removeNode(orderedNodes[i]);
			}
			else
			{
				BridgeNode(orderedNodes[i]);
			}
		}

		findAndLabelConnectedComponents();
		return;
	}

	//void BMetaGraph::SplitOperation(MetaEdge3d toSplit, std::vector<MetaEdge3d> connectedEdges)
	//{
	//	MetaE dupE;
	//	MetaV dupV0, dupV1;
	//	MetaE srcE;
	//	bool exists = false;
	//	boost::tie(srcE, exists) = boost::edge(toSplit.node0, toSplit.node1, *this);
	//	if (exists)
	//	{
	//		duplicateEdge(srcE, &mSkeleton, dupV0, dupV1, dupE);
	//		for (int i = 0; i < connectedEdges.size(); ++i)
	//		{
	//			MetaEdge3d srcEdge =connectedEdges[i];
	//			int srcNode, oldTarget, newTarget;
	//			//if the srcEdge node0 is equal to either of the nodes on the splitting edge, then 
	//			//we will break the connection between srcEdge node0 and node1, and create a new edge
	//			//between node0 and the new target node (determined by mapping of the duplicated nodes
	//			//to the old target node of the srcEdge
	//			if (srcEdge.node0 == toSplit.node0 || srcEdge.node0 == toSplit.node1)
	//			{
	//				srcNode = srcEdge.node1;
	//				oldTarget = srcEdge.node0;
	//				newTarget = mDuplicateNodeMap[dupV0] == oldTarget ? dupV0 : dupV1;
	//			}
	//			else if (srcEdge.node1 == toSplit.node0 || srcEdge.node1 == toSplit.node1)
	//			{
	//				srcNode = srcEdge.node0;
	//				oldTarget = srcEdge.node1;
	//				newTarget = mDuplicateNodeMap[dupV0] == oldTarget ? dupV0 : dupV1;
	//			}

	//			MetaE removeE;
	//			boost::tie(removeE, exists) = boost::edge(srcNode, oldTarget, *this);

	//			if (exists)
	//			{
	//				std::vector<SkelVert> verts = operator[](removeE).mVertices;
	//				removeEdge(removeE);
	//				addEdge(srcNode, newTarget, verts, &mSkeleton, false);
	//			}
	//			BMetaNode &newNode = operator[](newTarget);
	//			BMetaNode &goTowardsNode = operator[](srcNode);
	//			float deltaX = goTowardsNode.x() - newNode.x();
	//			float deltaY = goTowardsNode.y() - newNode.y();
	//			float deltaZ = goTowardsNode.z() - newNode.z();
	//			float deltaNorm = sqrt(deltaX*deltaX + deltaY * deltaY + deltaZ * deltaZ);
	//			float xProp = deltaX / deltaNorm;
	//			float yProp = deltaY / deltaNorm;
	//			float zProp = deltaZ / deltaNorm;

	//			float xDist = std::min(abs(xProp *toSplit.avgThickness * 2.0), abs(deltaX / 3.0));
	//			float yDist = std::min(abs(yProp *toSplit.avgThickness * 2.0), abs(deltaY / 3.0));
	//			float zDist = std::min(abs(zProp *toSplit.avgThickness * 2.0), abs(deltaZ / 3.0));

	//			newNode[0] += xProp < 0 ? -xDist : xDist;
	//			newNode[1] += yProp < 0 ? -yDist : yDist;
	//			newNode[2] += zProp < 0 ? -zDist : zDist;
	//			
	//		}

	//		//attempt to bridge the nodes at either end of both sides of the split
	//		//eg. if both the duplicated edge and the src edge will now have a single edge
	//		//connected at 1 end, then join those two edges together
	//		BridgeNode(toSplit.node0);
	//		BridgeNode(toSplit.node1);

	//		BridgeNode(dupV0);
	//		BridgeNode(dupV1);

	//		pythonUpToDate = false;
	//	}
	//	return;
	//}

	void BMetaGraph::SplitOperation()
	{
		if (!splitEdgeValid)
		{
			return;
		}
		bool sourceCovered = false, targetCovered = false;
		for (std::pair<MetaV, MetaV> edge : splitNeighbors)
		{
			if (edge.first == splitEdge.second || edge.second == splitEdge.second)
			{
				targetCovered = true;
			}
			if (edge.first == splitEdge.first || edge.second == splitEdge.first)
			{
				sourceCovered = true;
			}
		}
		if (!sourceCovered || !targetCovered)
		{
			return;
		}
		unselectEdge(splitEdge);
		for (std::pair<MetaV, MetaV> edge : splitNeighbors)
		{
			unselectEdge(edge);
		}


		MetaV node0 = splitEdge.first;
		MetaV node1 = splitEdge.second;
		//MetaE dupE;
		//MetaV dupV0, dupV1;
		MetaE srcE;

		bool exists = false;

		boost::tie(srcE, exists) = boost::edge(splitEdge.first, splitEdge.second, *this);

		if (!exists)
		{
			return;
		}

		BMetaEdge toSplit = operator[](srcE);

		float shiftX = 0, shiftY = 0, shiftZ = 0;

		for (std::pair<MetaV, MetaV> originalEdge : splitNeighbors)
		{

			MetaV connectionToSplit;
			//if the srcEdge node0 is equal to either of the nodes on the splitting edge, then 
			//we will break the connection between srcEdge node0 and node1, and create a new edge
			//between node0 and the new target node (determined by mapping of the duplicated nodes
			//to the old target node of the srcEdge
			if (originalEdge.first == splitEdge.first || originalEdge.first == splitEdge.second)
			{
				connectionToSplit = originalEdge.first;
			}
			else
			{
				connectionToSplit = originalEdge.second;
			}

			MetaE neighbor;
			boost::tie(neighbor, exists) = boost::edge(originalEdge.first, originalEdge.second, *this);
			BMetaEdge nedge = operator[](neighbor);
			SkelVert shiftTowards;
			if (vertNodeMap[nedge.start()] == connectionToSplit)
			{
				shiftTowards = nedge.mVertices[1];
			}
			else
			{
				shiftTowards = nedge.mVertices[nedge.mVertices.size() - 2];
			}

			float* goTowards = mSkeleton[shiftTowards].p;
			float* from = operator[](connectionToSplit).p;

			float deltaX = goTowards[0] - from[0];
			float deltaY = goTowards[1] - from[1];
			float deltaZ = goTowards[2] - from[2];
			shiftX += deltaX *0.75;
			shiftY += deltaY * 0.75;
			shiftZ += deltaZ * 0.75;
		}

		Point3d shift = Point3d(shiftX, shiftY, shiftZ);
		std::vector<SkelVert> splitVerts = toSplit.mVertices;
		std::vector<RootAttributes> splitEdges = toSplit.mEdges;
		std::map<MetaV, MetaV> originalToDuplicateMap = std::map<MetaV, MetaV>();

		std::vector<SkelVert> newVertices = {};

		for (int i = 0; i < splitVerts.size(); ++i)
		{
			Point3d src = mSkeleton[splitVerts[i]];
			Point3d dst = src + shift;
			newVertices.push_back(mSkeleton.addVertex(dst));
		}

		for (int i = 0; i < newVertices.size() - 1; ++i)
		{
			SkelEdge added = mSkeleton.addEdge(newVertices[i], newVertices[i + 1], splitEdges[i], edgeVertices);
		}

		MetaV v0 = addNode(newVertices[0], &mSkeleton);
		MetaV v1 = addNode(newVertices.back(), &mSkeleton);

		originalToDuplicateMap[vertNodeMap[splitVerts[0]]] = v0;
		originalToDuplicateMap[vertNodeMap[splitVerts.back()]] = v1;

		MetaE addedEdge = addEdge(v0, v1, newVertices, &mSkeleton, false);

		for (std::pair<MetaV, MetaV> originalEdge : splitNeighbors)
		{

			MetaV splitNode, remainingNode;
			//if the srcEdge node0 is equal to either of the nodes on the splitting edge, then 
			//we will break the connection between srcEdge node0 and node1, and create a new edge
			//between node0 and the new target node (determined by mapping of the duplicated nodes
			//to the old target node of the srcEdge

			
			MetaE e;
			bool exists = false;
			boost::tie(e, exists) = boost::edge(originalEdge.first, originalEdge.second, *this);

			std::vector<SkelVert> originalVerts = {};
			for (SkelVert v : operator[](e).mVertices)
			{
				originalVerts.push_back(v);
			}
			std::vector<RootAttributes> originalEdges = operator[](e).mEdges;

			if (originalEdge.first == splitEdge.first || originalEdge.first == splitEdge.second)
			{
				splitNode = originalEdge.first;
				remainingNode = originalEdge.second;
			}
			else
			{
				splitNode = originalEdge.second;
				remainingNode = originalEdge.first;
			}

			if (vertNodeMap[originalVerts[0]] == splitNode)
			{
				originalVerts[0] = operator[](originalToDuplicateMap[splitNode]).mSrcVert;
				mSkeleton.addEdge(originalVerts[0], originalVerts[1], originalEdges[0], edgeVertices);
			}
			else
			{
				originalVerts[originalVerts.size() - 1] = operator[](originalToDuplicateMap[splitNode]).mSrcVert;
				mSkeleton.addEdge(originalVerts[originalVerts.size() - 1], originalVerts[originalVerts.size() - 2], originalEdges.back(), edgeVertices);
			}
			

			removeEdge(e);
			addEdge(remainingNode, originalToDuplicateMap[splitNode], originalVerts, &mSkeleton, false, false);
		}


		//duplicateEdge(srcE, &mSkeleton, dupV0, dupV1, dupE);
		//for (std::pair<MetaV, MetaV> originalEdge : splitNeighbors)
		//{
		//	int srcNode, oldTarget, newTarget;
		//	//if the srcEdge node0 is equal to either of the nodes on the splitting edge, then 
		//	//we will break the connection between srcEdge node0 and node1, and create a new edge
		//	//between node0 and the new target node (determined by mapping of the duplicated nodes
		//	//to the old target node of the srcEdge
		//	if (originalEdge.first == splitEdge.first || originalEdge.first == splitEdge.second)
		//	{
		//		srcNode = originalEdge.second;
		//		oldTarget = originalEdge.first;
		//		newTarget = mDuplicateNodeMap[dupV0] == oldTarget ? dupV0 : dupV1;
		//	}
		//	else if (originalEdge.second == splitEdge.first || originalEdge.second == splitEdge.second)
		//	{
		//		srcNode = originalEdge.first;
		//		oldTarget = originalEdge.second;
		//		newTarget = mDuplicateNodeMap[dupV0] == oldTarget ? dupV0 : dupV1;
		//	}

		//	MetaE removeE;
		//	boost::tie(removeE, exists) = boost::edge(srcNode, oldTarget, *this);

		//	if (exists)
		//	{
		//		std::vector<SkelVert> verts = operator[](removeE).mVertices;
		//		removeEdge(removeE);
		//		addEdge(srcNode, newTarget, verts, &mSkeleton, false);
		//	}
		//	BMetaNode &newNode = operator[](newTarget);
		//	BMetaNode &goTowardsNode = operator[](srcNode);
		//	float deltaX = goTowardsNode.x() - newNode.x();
		//	float deltaY = goTowardsNode.y() - newNode.y();
		//	float deltaZ = goTowardsNode.z() - newNode.z();
		//	float deltaNorm = sqrt(deltaX*deltaX + deltaY * deltaY + deltaZ * deltaZ);
		//	float xProp = deltaX / deltaNorm;
		//	float yProp = deltaY / deltaNorm;
		//	float zProp = deltaZ / deltaNorm;

		//	float xDist = std::min(abs(xProp *toSplit.averageThickness * 2.0), abs(deltaX / 3.0));
		//	float yDist = std::min(abs(yProp *toSplit.averageThickness * 2.0), abs(deltaY / 3.0));
		//	float zDist = std::min(abs(zProp *toSplit.averageThickness * 2.0), abs(deltaZ / 3.0));

		//	newNode[0] += xProp < 0 ? -xDist : xDist;
		//	newNode[1] += yProp < 0 ? -yDist : yDist;
		//	newNode[2] += zProp < 0 ? -zDist : zDist;
		//			
		//}

		//attempt to bridge the nodes at either end of both sides of the split
		//eg. if both the duplicated edge and the src edge will now have a single edge
		//connected at 1 end, then join those two edges together

		std::vector<MetaV> orderedNodes = { node0, node1, v0, v1 };
		std::sort(orderedNodes.begin(), orderedNodes.end());

		/*BridgeNode(orderedNodes[3]);
		BridgeNode(orderedNodes[2]);

		BridgeNode(orderedNodes[1]);
		BridgeNode(orderedNodes[0]);*/

		
		findAndLabelConnectedComponents();

		for (int i = orderedNodes.size() - 1; i >= 0; --i)
		{
			BridgeNode(orderedNodes[i]);
		}


		return;
	}


	void BMetaGraph::PromoteOperation(SkelVert toPromote)
	{
		std::cout << "Promote operation triggered " << std::endl;
		metaEdgeIter mei = boost::edges(*this);
		bool edgeFound = false;
		for (; mei.first != mei.second; ++mei)
		{
			BMetaEdge *edge = &operator[](*mei.first);

			for (SkelVert v : edge->mVertices)
			{
				if (v == toPromote)
				{
					edgeFound = true;
					break;
				}
			}
			if (edgeFound)
			{
				MetaV leftNode = vertNodeMap[edge->start()];
				MetaV rightNode = vertNodeMap[edge->end()];
				MetaV myNode = addNode(toPromote, &mSkeleton);
				std::pair<BMetaEdge, BMetaEdge> split = edge->split(toPromote, &mSkeleton);
				MetaE leftEdge, rightEdge;
				bool success = false;
				boost::tie(leftEdge, success) = add_edge(leftNode, myNode, *this);
				if (success)
				{
					operator[](leftEdge) = split.first;
				}
				boost::tie(rightEdge, success) = add_edge(myNode, rightNode, *this);
				if (success)
				{
					operator[](rightEdge) = split.second;
				}

				updateNodeDegrees();
				removeEdgeNoBridging(*mei.first);

				break;
			}
		}
	}

	bool BMetaGraph::BridgeNode(MetaV nodeToBridge)
	{
		std::cout << "Beginning to bridge " << std::endl;
		if (nodeToBridge >= m_vertices.size())
		{
			std::cout << "Invalid index " << std::endl;
			return false;
		}
		//return;
		if (boost::degree(nodeToBridge, *this) == 2)
		{
			BMetaGraph::adjacency_iterator adjIt, adjEnd;
			boost::tie(adjIt, adjEnd) = boost::adjacent_vertices(nodeToBridge, *this);
			MetaV v1 = *adjIt;
			MetaV v2 = *(adjIt + 1);
			MetaE e1, e2;
			bool exists = false;
			boost::tie(e1, exists) = boost::edge(nodeToBridge, v1, *this);
			boost::tie(e2, exists) = boost::edge(nodeToBridge, v2, *this);
			if (v1 == v2 || v1 == nodeToBridge || v2 == nodeToBridge)
			{
				return false;
			}
			BMetaEdge joinedEdge = operator[](e1).join(operator[](e2), &mSkeleton);
			MetaE addedEdge = addEdge(v1, v2, joinedEdge.mVertices, &mSkeleton, false);

			boost::tie(adjIt, adjEnd) = boost::adjacent_vertices(nodeToBridge, *this);
			std::vector<MetaE> toRemove = {};
			for (; adjIt != adjEnd; ++adjIt)
			{
				MetaE edge;
				boost::tie(edge, exists) = boost::edge(nodeToBridge, *adjIt, *this);
				
				if (exists)
				{
					toRemove.push_back(edge);
				}
				
			}
			for (MetaE edge : toRemove)
			{
				removeEdge(edge);
			}
			removeNode(nodeToBridge);
			findAndLabelConnectedComponents();
			operator[](addedEdge).updateColors(edgeOptions, edgeColors);

			return true;
		}
		return false;
	}

	/*void BMetaGraph::BridgeNodeSweeper()
	{
		return;
		bool bridged = true;

		while (bridged)
		{
			bridged = false;

			metaVertIter mvi = boost::vertices(*this);
			--mvi;
			--mvi.first;
			for (; mvi.first != mvi.second; --mvi)
			{
				bridged = BridgeNode(*mvi.second);
				if (bridged)
				{
					break;
				}
			}
		}
	}*/

	int BMetaGraph::GetEdgeComponent(MetaE edge)
	{
		BMetaEdge src = operator[](edge);
		//ensure that both endpoint vertices are mapped as metanodes

		if (vertNodeMap.count(src.start()) == 0 ||
			vertNodeMap.count(src.end()) == 0)
		{
			//std::cout << "The endpoints of this metaedge are not mapped as metanodes" << std::endl;
			return -1;
		}

		MetaV node1 = vertNodeMap[src.start()];
		MetaV node2 = vertNodeMap[src.end()];

		int component1 = operator[](node1).connectedComponent;
		int component2 = operator[](node2).connectedComponent;

		if (component1 != component2)
		{
			ConnectComponents(component1, component2);
		}

		return operator[](node1).connectedComponent;
	}


	void BMetaGraph::bridgeUtil(int u, bool visited[], int disc[],
		int low[], int parent[], std::vector<std::pair<MetaV, MetaV>> &bridgeNodes)
	{
		// A static variable is used for simplicity, we can 
		// avoid use of static variable by passing a pointer.
		static int time = 0;

		// Mark the current node as visited
		visited[u] = true;

		// Initialize discovery time and low value
		disc[u] = low[u] = ++time;

		// Go through all vertices aadjacent to this
		BMetaGraph::adjacency_iterator adjIt, adjEnd;
		boost::tie(adjIt, adjEnd) = boost::adjacent_vertices(u, *this);
		
		for (; adjIt != adjEnd; ++adjIt)
		{
			size_t v = *adjIt;

			if (!visited[v])
			{
				parent[v] = u;
				bridgeUtil(v, visited, disc, low, parent, bridgeNodes);

				low[u] = std::min(low[u], low[v]);

				if (low[v] > disc[u])
				{
					bridgeNodes.push_back(std::make_pair(u, v));
					//std::cout << u << " " << v << std::endl;
				}
			}
			else if (v != parent[u])
			{
				low[u] = std::min(low[u], disc[v]);
			}
		}
	}

	// DFS based function to find all bridges. It uses recursive 
	// function bridgeUtil()
	void BMetaGraph::findBridges()
	{
		// Mark all the vertices as not visited
		bool *visited = new bool[m_vertices.size()];
		int *disc = new int[m_vertices.size()];
		int *low = new int[m_vertices.size()];
		int *parent = new int[m_vertices.size()];



		// Initialize parent and visited arrays
		for (int i = 0; i < m_vertices.size(); i++)
		{
			parent[i] = -1;
			visited[i] = false;
		}

		metaEdgeIter mei = boost::edges(*this);
		for (; mei.first != mei.second; ++mei)
		{
			operator[](*mei.first).isBridge = false;
		}


		std::vector<std::pair<MetaV, MetaV>> bridgeNodes = {};

		std::cout << "==================Bridges exist between the following vertices==================" << std::endl;
		// Call the recursive helper function to find Bridges
		// in DFS tree rooted with vertex 'i'
		for (int i = 0; i < m_vertices.size(); i++)
			if (visited[i] == false)
				bridgeUtil(i, visited, disc, low, parent, bridgeNodes);

		
		for each(auto nodepair in bridgeNodes)
		{
			MetaE e;
			bool exists;
			boost::tie(e, exists) = boost::edge(nodepair.first, nodepair.second, *this);

			if (exists)
			{
				operator[](e).isBridge = true;
			}
		}

		//nonBridgeEdgeIndices = {};
		//bridgeEdgeIndices = {};

		//mei = boost::edges(*this);
		//for (; mei.first != mei.second; ++mei)
		//{
		//	if (operator[](*mei.first).isBridge)
		//	{
		//		operator[](*mei.first).addToIndicesList(bridgeEdgeIndices);
		//	}
		//	else
		//	{
		//		operator[](*mei.first).addToIndicesList(nonBridgeEdgeIndices);
		//	}
		//}



		delete[] visited;
		delete[] disc;
		delete[] low;
		delete[] parent;
		std::cout << "==================End Bridges==================" << std::endl;
	}

	
}

//
//PyMetaGraph::PyMetaGraph(std::string filename)
//{
//	mGraph = Roots::BMetaGraph();
//	mGraph.loadFromFile(filename);
//	mSkeleton = PySkeleton(&mGraph.mSkeleton);
//	mGraph.findAndLabelConnectedComponents();
//	reload();
//}
//
//void PyMetaGraph::initializeFromSkeleton()
//{
//	mGraph.initializeFromSkeleton();
//	mGraph.findAndLabelConnectedComponents();
//	reload();
//}
//
//void PyMetaGraph::labelComponents()
//{
//	mGraph.findAndLabelConnectedComponents();
//	reload();
//}
//
//void PyMetaGraph::joinOperation(int v0, int v1)
//{
//	//std::cout << "Attempting to join vertices " << v0 << " " << v1 << std::endl;
//	Point3d p1 = mGraph.mSkeleton[v0];
//	Point3d p2 = mGraph.mSkeleton[v1];
//	mGraph.JoinOperation(p1, p2);
//	//mSkeleton = PySkeleton(&mGraph.mSkeleton);
//	reload();
//
//}
//
//void PyMetaGraph::breakOperation(Roots::MetaEdge3d edge)
//{
//	//std::cout << "Attempting to break edge between " << edge.node0 << " " << edge.node1 << std::endl;
//	std::cout << "Breaking edge between " << edge.node0 << " and " << edge.node1 << ". This is the " << edge.order << "th edge." << std::endl;
//	mGraph.BreakOperation(edge);
//	mGraph.findAndLabelConnectedComponents();
//	reload();
//}
//
//void PyMetaGraph::splitOperation(Roots::MetaEdge3d toSplit, boost::python::list connectedEdgeSet)
//{
//	std::vector<Roots::MetaEdge3d> edgeSet = toStdVector<Roots::MetaEdge3d>(connectedEdgeSet);
//	mGraph.SplitOperation(toSplit, edgeSet);
//	this->mSkeleton = PySkeleton(&mGraph.mSkeleton);
//	mGraph.findAndLabelConnectedComponents();
//	reload();
//}
//
//void PyMetaGraph::reload()
//{
//	if (!mGraph.pythonUpToDate)
//	{
//		//mGraph.findAndLabelConnectedComponents();
//		//std::cout << "Entering boost python update loop " << std::endl;
//		mGraph.FindMinimumSpanningTree();
//		//std::cout << "Num edges to break " << mGraph.GetNumEdgesToBreak() << std::endl;
//		numEdgesToBreak = mGraph.GetNumEdgesToBreak();
//		mGraph.findBridges();
//		mMetaNodeLocations = boost::python::list();
//		mMetaEdgeConnections = boost::python::list();
//		std::vector<Roots::MetaNode3d> metaNodeLocs = {};
//		std::vector<Roots::MetaEdge3d> metaEdgeCons = {};
//		componentNodeMap = boost::python::dict();
//		componentEdgeMap = boost::python::dict();
//		componentNodesOfInterestMap = boost::python::dict();
//
//		//stand in std format maps to edit in place before saving to python equivalents
//		std::map<int, boost::python::list> tempComponentNodeMap;
//		std::map<int, boost::python::list> tempComponentNodesOfInterestMap;
//		std::map<int, boost::python::list> tempComponentEdgeMap;
//		
//
//		metaVertIter mvi = boost::vertices(mGraph);
//		int vertI = 0;
//		for (; mvi.first != mvi.second; ++mvi, ++vertI)
//		{
//			//std::cout << mGraph.m_vertices.size() << std::endl;
//			//std::cout << "Trying to find vertex to add : ";
//			//std::cout << *mvi.first << std::endl;
//
//			//std::cout << "Getting meta node " << std::endl;
//			//Roots::BMetaNode curNode = mGraph[*mvi.first];
//			//std::cout << curNode.mSrcVert << std::endl;
//			//std::cout << "Successfully accessed meta node " << std::endl;
//
//			int nodeDegree = boost::degree(*mvi.first, mGraph);
//
//			Roots::MetaNode3d toAdd = Roots::MetaNode3d(mGraph[*mvi.first], &mGraph.mSkeleton, vertI, nodeDegree);
//			//std::cout << toAdd << std::endl;
//			metaNodeLocs.push_back(toAdd);
//
//			if (tempComponentNodeMap.count(toAdd.connectedComponent) == 0)
//			{
//				boost::python::list toStart = boost::python::list();
//				tempComponentNodeMap[toAdd.connectedComponent] = toStart;
//				tempComponentNodesOfInterestMap[toAdd.connectedComponent] = boost::python::list();
//			}
//			//std::cout << "Finding node degree for " << *mvi.first << "th node.  It is : ";
//			
//			//std::cout << nodeDegree << std::endl;
//			if (nodeDegree < 2)
//			{
//				tempComponentNodesOfInterestMap[toAdd.connectedComponent].append<int>(vertI);
//			}
//			//std::cout << "Adding node to tempConnectedComponentMap " << std::endl;
//			tempComponentNodeMap[toAdd.connectedComponent].append<int>(vertI);
//		}
//		//std::cout << "completed mapping nodes and components " << std::endl;
//
//		//std::cout << "Beginning to map edges and components" << std::endl;
//		metaEdgeIter mei = boost::edges(mGraph);
//		int edgeI = 0;
//		for (; mei.first != mei.second; ++mei, ++edgeI)
//		{
//			MetaE edge = *mei.first;
//			Roots::BMetaEdge srcEdge = mGraph[edge];
//			Roots::MetaEdge3d toAdd = Roots::MetaEdge3d(edge.m_source, edge.m_target,
//				srcEdge.getAvgThickness(), srcEdge.averageWidth, mGraph.GetEdgeComponent(edge), edgeI, srcEdge.mEdges, srcEdge.isBridge);
//			metaEdgeCons.push_back(toAdd);
//			//mMetaEdgeConnections.append<MetaEdge3d>(toAdd);
//
//			if (tempComponentEdgeMap.count(toAdd.connectedComponent) == 0)
//			{
//				boost::python::list toStart = boost::python::list();
//				tempComponentEdgeMap[toAdd.connectedComponent] = toStart;
//			}
//			tempComponentEdgeMap[toAdd.connectedComponent].append<int>(edgeI);
//		}
//
//		for each(std::pair<int, boost::python::list> nodeMap in tempComponentNodeMap)
//		{
//			componentNodeMap[nodeMap.first] = nodeMap.second;
//		}
//
//		for each(std::pair<int, boost::python::list> nodeMap in tempComponentNodesOfInterestMap)
//		{
//			componentNodesOfInterestMap[nodeMap.first] = nodeMap.second;
//		}
//
//		for each(std::pair<int, boost::python::list> edgeMap in tempComponentEdgeMap)
//		{
//			componentEdgeMap[edgeMap.first] = edgeMap.second;
//		}
//
//		mMetaNodeLocations = toPythonList<Roots::MetaNode3d>(metaNodeLocs);
//
//		mMetaEdgeConnections = toPythonList<Roots::MetaEdge3d>(metaEdgeCons);
//
//		mGraph.pythonUpToDate = true;
//
//
//	}
//	//std::cout << "completed boost python update loop " << std::endl;
//}