#include "BoostMetaGraph.h"

//namespace
//{
//	typedef boost::graph_traits<Roots::BMetaGraph>::vertex_iterator vertIter;
//	typedef boost::graph_traits<Roots::BMetaGraph>::edge_iterator edgeIter;
//	using namespace boost;
//	struct metaVertIter : std::pair<vertIter, vertIter>
//	{
//		metaVertIter(const std::pair<vertIter, vertIter> &other)
//			:std::pair<vertIter, vertIter>(other)
//		{
//		}
//		metaVertIter operator++()
//		{
//			++this->first;
//			return *this;
//		}
//		metaVertIter operator--()
//		{
//			--this->second;
//			return *this;
//		}
//	};
//
//	struct metaEdgeIter : std::pair<edgeIter, edgeIter>
//	{
//		metaEdgeIter(const std::pair<edgeIter, edgeIter> &other)
//			:std::pair<edgeIter, edgeIter>(other)
//		{
//		}
//		metaEdgeIter operator++()
//		{
//			++this->first;
//			return *this;
//		}
//
//		metaEdgeIter operator--()
//		{
//			--this->second;
//			return *this;
//		}
//	};
//
//	template <class T>
//	boost::python::list toPythonList(std::vector<T> vector) {
//		typename std::vector<T>::iterator iter;
//		boost::python::list list;
//		for (iter = vector.begin(); iter != vector.end(); ++iter) {
//			list.append(*iter);
//		}
//		return list;
//	}
//
//	template <class T>
//	std::vector<T> toStdVector(boost::python::list list)
//	{
//		std::vector<T> result = {};
//		for (int i = 0; i < boost::python::len(list); ++i)
//		{
//			result.push_back(boost::python::extract<T>(list[i]));
//		}
//		return result;
//	}
//
//	GLfloat defaultColor[3] = { 1.0, 0.0, 0.0 };
//}


namespace Roots
{

	//////////////////////////////////////////////EDGES////////////////////////////////////////////
	EdgeVisualizationOptions::EdgeVisualizationOptions()
	{

		show = true;
		scale = 1.0;
		magnifyNonBridges = false;
		showOnlyNonBridges = false;

		minColorCutoff = 0.0;
		maxColorCutoff = 1.0;
		colorization = ColorizationOptions::ByThickness;
		flatSelectionColor[0] = 1.0;
		flatSelectionColor[1] = 1.0;
		flatSelectionColor[2] = 1.0;
		flatSelectionColor[3] = 1.0;
		heatmap = { { 1.0, 0.0, 0.0, 1.0 },{ 1.0, 0.0, 0.0, 1.0 } };

	}

	void BMetaGraph::showEdges(bool doShow)
	{
		edgeOptions.show = doShow;
		std::cout << "Do show edges " << doShow << std::endl;
	}

	void BMetaGraph::setEdgeScale(float scale)
	{
		edgeOptions.scale = scale;
		std::cout << "Set edge scale " << scale << std::endl;
	}

	void BMetaGraph::magnifyNonBridges(bool doMagnify)
	{
		edgeOptions.magnifyNonBridges = doMagnify;
		std::cout << "Magnify non-bridges " << std::to_string(doMagnify) << std::endl;
	}

	void BMetaGraph::showOnlyNonBridges(bool showOnly)
	{
		edgeOptions.showOnlyNonBridges = showOnly;
		std::cout << "Show only non-bridges " << std::to_string(showOnly) << std::endl;
	}

	void BMetaGraph::assignEdgeHeatMap(boost::python::list heatmap)
	{
		edgeOptions.heatmap = {};
		for (int i = 0; i < boost::python::len(heatmap); ++i)
		{
			std::vector<GLfloat> vectorColor = {};
			boost::python::list pythonColor = boost::python::extract<boost::python::list>(heatmap[i]);
			for (int c = 0; c < 3; ++c)
			{
				vectorColor.push_back(boost::python::extract<float>(pythonColor[c]));
			}
			edgeOptions.heatmap.push_back(vectorColor);
		}
		metaEdgeIter mei = boost::edges(*this);
		for (; mei.first != mei.second; ++mei)
		{
			operator[](*mei.first).updateColors(edgeOptions, edgeColors);
		}

	}

	void BMetaGraph::colorizeEdgesByThickness()
	{
		edgeOptions.colorization = ColorizationOptions::ByThickness;
		std::cout << "Colorize edges by thickness " << std::endl;
	}
	void BMetaGraph::colorizeEdgesByWidth()
	{
		edgeOptions.colorization = ColorizationOptions::ByWidth;
		std::cout << "Colorize edges by width " << std::endl;
	}
	void BMetaGraph::colorizeEdgesByRatio()
	{
		edgeOptions.colorization = ColorizationOptions::ByRatio;
		std::cout << "Colorize edges by ratio " << std::endl;
	}
	void BMetaGraph::colorizeEdgesByComponent()
	{
		edgeOptions.colorization = ColorizationOptions::ByComponent;
		std::cout << "Colorize edges by component " << std::endl;
	}
	
	void BMetaGraph::setEdgeColorFloor(float val)
	{
		edgeOptions.minColorCutoff = val;

		metaEdgeIter mei = boost::edges(*this);
		for (; mei.first != mei.second; ++mei)
		{
			operator[](*mei.first).updateColors(edgeOptions, edgeColors);
		}
	}
	void BMetaGraph::setEdgeColorCeiling(float val)
	{
		edgeOptions.maxColorCutoff = val;
		if (edgeOptions.maxColorCutoff == edgeOptions.minColorCutoff)
		{
			edgeOptions.maxColorCutoff += 0.001;
		}

		metaEdgeIter mei = boost::edges(*this);
		for (; mei.first != mei.second; ++mei)
		{
			operator[](*mei.first).updateColors(edgeOptions, edgeColors);
		}
	}

	void BMetaGraph::setEdgeSelectionColor(float r, float g, float b)
	{
		edgeOptions.flatSelectionColor[0] = r;
		edgeOptions.flatSelectionColor[1] = g;
		edgeOptions.flatSelectionColor[2] = b;
		return;
	}


	/////////////////////////////////////////////NODES/////////////////////////////////////////////
	NodeVisualizationOptions::NodeVisualizationOptions()
	{

		showEndpoints = true;
		endpointScale = 0.5f;
		showJunctions = true;
		junctionScale = 0.5f;

		useConstantNodecolor = false;
		constantNodeColor[0] = 1.0;
		constantNodeColor[1] = 0.5;
		constantNodeColor[2] = 0.0;
		constantNodeColor[3] = 1.0;

		minColorCutoff = 0.0f;
		maxColorCutoff = 1.0;
		colorization = ColorizationOptions::ByDegree;
		heatmap = { {1.0, 0.0, 0.0, 1.0},  {1.0, 0.0, 0.0, 1.0} };
	}

	void BMetaGraph::showEndpoints(bool doShow)
	{
		nodeOptions.showEndpoints = doShow;
		std::cout << "Do show endpoints " << doShow << std::endl;
	}

	void BMetaGraph::showJunctions(bool doShow)
	{
		nodeOptions.showJunctions = doShow;
		std::cout << "Do show junctions " << doShow << std::endl;
	}

	void BMetaGraph::setEndpointScale(float scale)
	{
		nodeOptions.endpointScale = scale;
	}

	void BMetaGraph::setJunctionScale(float scale)
	{
		nodeOptions.junctionScale = scale;
	}

	void BMetaGraph::setConstantNodeColor(float r, float g, float b)
	{
		nodeOptions.constantNodeColor[0] = r;
		nodeOptions.constantNodeColor[1] = g;
		nodeOptions.constantNodeColor[2] = b;
		nodeOptions.constantNodeColor[3] = 1.0;
	}

	void BMetaGraph::assignNodeHeatMap(boost::python::list heatmap)
	{
		//nodeOptions.heatmap = heatmap;
		nodeOptions.heatmap = {};

		for (int i = 0; i < boost::python::len(heatmap); ++i)
		{
			std::vector<GLfloat> vectorColor = {};
			boost::python::list pythonColor = boost::python::extract<boost::python::list>(heatmap[i]);

			for (int c = 0; c < 3; ++c)
			{
				vectorColor.push_back(boost::python::extract<float>(pythonColor[c]));
			}
			
			nodeOptions.heatmap.push_back(vectorColor);
		}

		metaVertIter mvi = boost::vertices(*this);
		for (; mvi.first != mvi.second; ++mvi)
		{
			MetaV  node = *mvi.first;
			operator[](node).updateColors(nodeOptions);
		}
	}

	void BMetaGraph::colorizeNodesByThickness()
	{
		nodeOptions.colorization = ColorizationOptions::ByThickness;
		metaVertIter mvi = boost::vertices(*this);
		for (; mvi.first != mvi.second; ++mvi)
		{
			MetaV  node = *mvi.first;
			operator[](node).currentColor = operator[](node).glThicknessColor;
		}
		std::cout << "Colorize nodes by thickness " << std::endl;
	}
	void BMetaGraph::colorizeNodesByWidth()
	{
		nodeOptions.colorization = ColorizationOptions::ByWidth;
		metaVertIter mvi = boost::vertices(*this);
		for (; mvi.first != mvi.second; ++mvi)
		{
			MetaV  node = *mvi.first;
			operator[](node).currentColor = operator[](node).glWidthColor;
		}
		std::cout << "Colorize ndes by width " << std::endl;
	}
	void BMetaGraph::colorizeNodesByDegree()
	{
		nodeOptions.colorization = ColorizationOptions::ByDegree;
		metaVertIter mvi = boost::vertices(*this);
		for (; mvi.first != mvi.second; ++mvi)
		{
			MetaV  node = *mvi.first;
			operator[](node).currentColor = operator[](node).glDegreeColor;
		}
		std::cout << "Colorize nodes by degree " << std::endl;
	}
	void BMetaGraph::colorizeNodesByComponent()
	{
		nodeOptions.colorization = ColorizationOptions::ByComponent;
		metaVertIter mvi = boost::vertices(*this);
		for (; mvi.first != mvi.second; ++mvi)
		{
			MetaV  node = *mvi.first;
			operator[](node).currentColor = operator[](node).glComponentColor;
		}
		std::cout << "Colorize nodes by component " << std::endl;
	}

	void BMetaGraph::colorizeNodesByConstantColor()
	{
		nodeOptions.colorization = ColorizationOptions::Other;
		metaVertIter mvi = boost::vertices(*this);
		for (; mvi.first != mvi.second; ++mvi)
		{
			MetaV  node = *mvi.first;
			operator[](node).currentColor = nodeOptions.constantNodeColor;
		}
		std::cout << "Colorize nodes by constant color" << std::endl;
	}

	void BMetaGraph::setNodeColorFloor(float val)
	{
		nodeOptions.minColorCutoff = val;
		metaVertIter mvi = boost::vertices(*this);
		for (; mvi.first != mvi.second; ++mvi)
		{
			MetaV  node = *mvi.first;
			operator[](node).updateColors(nodeOptions);
		}
	}
	void BMetaGraph::setNodeColorCeiling(float val)
	{
		nodeOptions.maxColorCutoff = val;
		if (nodeOptions.maxColorCutoff == nodeOptions.minColorCutoff)
		{
			nodeOptions.maxColorCutoff += 0.001;
		}
		metaVertIter mvi = boost::vertices(*this);
		for (; mvi.first != mvi.second; ++mvi)
		{
			MetaV  node = *mvi.first;
			operator[](node).updateColors(nodeOptions);
		}
	}

	//////////////////////////////////////////COMPONENT OPTIONS////////////////////////////////////

	void BMetaGraph::setDisplayOnlySelectedComponents(bool displayOnlySelectedComponents)
	{
		onlyDisplaySelectedComponents = displayOnlySelectedComponents;
	}

	void BMetaGraph::setComponent1(int component)
	{
		selectedComponent1 = component;
		std::cout << "Component 1 set to " << component << std::endl;
	}

	void BMetaGraph::setComponent2(int component)
	{
		selectedComponent2 = component;
		std::cout << "Component 2 set to " << component << std::endl;
	}

	void BMetaGraph::setShowBoundingBoxes(bool doShow)
	{
		showBoundingBoxes = doShow;
		std::cout << "Show bounding boxes set to " << doShow << std::endl;
	}

	/////////////////////////////////////////DRAWING///////////////////////////////////////////////

	void BMetaGraph::drawEdges()
	{
		if (useArcball)
		{
			glTranslatef(0, 0, -mSkeleton.mRadius * 2);
			arcball_rotate();
		}
		
		if (!edgeOptions.show || !isLoaded)
		{
			return;
		}
		GLfloat *colorArray;
		switch (edgeOptions.colorization)
		{
		case ColorizationOptions::ByThickness:
			colorArray = &edgeColors[0][0];
			break;
		case ColorizationOptions::ByWidth:
			colorArray = &edgeColors[1][0];
			break;
		case ColorizationOptions::ByRatio:
			colorArray = &edgeColors[2][0];
			break;
		case ColorizationOptions::ByComponent:
			colorArray = &edgeColors[3][0];
			break;
		default:
			colorArray = &edgeColors[0][0];
			break;
		}

		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_COLOR_ARRAY);
		glVertexPointer(3, GL_FLOAT, 0, &edgeVertices[0]);
		glColorPointer(3, GL_FLOAT, 0, colorArray);

		metaEdgeIter mei = boost::edges(*this);

		for (; mei.first != mei.second; ++mei)
		{
			BMetaEdge *edge = &operator[](*mei.first);
			bool isDrawn = false;
			if (onlyDisplaySelectedComponents)
			{
				if (edge->connectedComponent != selectedComponent1 && edge->connectedComponent != selectedComponent2)
				{
					continue;
				}
			}
			if (edgeOptions.showOnlyNonBridges)
			{
				if (edge->isBridge)
				{
					continue;
				}
			}
			if (edge->isSelected)
			{
				glLineWidth(edgeOptions.scale * edgeSelectionScaling);
				glDrawElements(GL_LINES, edge->indicesList.size(), GL_UNSIGNED_INT, &edge->indicesList[0]);
				isDrawn = true;
			}
			else if (edgeOptions.magnifyNonBridges)
			{
				if (!edge->isBridge)
				{
					glLineWidth(edgeOptions.scale *nonBridgeScaling);
					glDrawElements(GL_LINES, edge->indicesList.size(), GL_UNSIGNED_INT, &edge->indicesList[0]);
					isDrawn = true;
				}
			}
			if (!isDrawn)
			{
				glLineWidth(edgeOptions.scale);
				glDrawElements(GL_LINES, edge->indicesList.size(), GL_UNSIGNED_INT, &edge->indicesList[0]);
			}
		}

		//if (edgeOptions.magnifyNonBridges)
		//{
		//	glLineWidth(edgeOptions.scale * nonBridgeScaling);
		//	glDrawElements(GL_LINES, nonBridgeEdgeIndices.size(), GL_UNSIGNED_INT, &nonBridgeEdgeIndices[0]);

		//	glLineWidth(edgeOptions.scale);
		//	glDrawElements(GL_LINES, bridgeEdgeIndices.size(), GL_UNSIGNED_INT, &bridgeEdgeIndices[0]);

		//	glLineWidth(edgeOptions.scale * selectionScaling);
		//	glDrawElements(GL_LINES, selectionEdgeIndices.size(), GL_UNSIGNED_INT, &selectionEdgeIndices[0]);
		//}
		//else if(onlyDisplaySelectedComponents)
		//{
		//	glLineWidth(edgeOptions.scale);
		//	glDrawElements(GL_LINES, componentIndices[selectedComponent1].size(), GL_UNSIGNED_INT, &componentIndices[selectedComponent1][0]);
		//	glDrawElements(GL_LINES, componentIndices[selectedComponent2].size(), GL_UNSIGNED_INT, &componentIndices[selectedComponent2][0]);
		//	
		//	glLineWidth(edgeOptions.scale*selectionScaling);
		//	glDrawElements(GL_LINES, selectionEdgeIndices.size(), GL_UNSIGNED_INT, &selectionEdgeIndices[0]);

		//}
		//else
		//{
		//	glLineWidth(edgeOptions.scale);
		//	glDrawElements(GL_LINES, baseEdgeIndices.size(), GL_UNSIGNED_INT, &baseEdgeIndices[0]);
		//}

		glDisableClientState(GL_VERTEX_ARRAY);
		glDisableClientState(GL_COLOR_ARRAY);
	}

	void BMetaGraph::drawNodes()
	{
		
		metaVertIter mvi = boost::vertices(*this);
		
		GLfloat *nodeColor;
		for (; mvi.first != mvi.second; ++mvi)
		{
			BMetaNode *node = &operator[](*mvi.first);
			int degree = boost::degree(*mvi.first, *this);
			//int degree = node->degree;

			if ((degree == 1 && !nodeOptions.showEndpoints) || (degree > 2 && !nodeOptions.showJunctions) || degree == 2)
			{
				continue;
			}

			

			bool isSelected = false;

			float scale = 0.0;

			if (onlyDisplaySelectedComponents)
			{
				if (node->connectedComponent != selectedComponent1 && node->connectedComponent != selectedComponent2)
				{
					continue;
				}
			}
			//if (degree == 0)
			//{
			//	node->updateColors(nodeOptions);
			//}

			if (nodeOptions.colorization == ColorizationOptions::ByThickness)
			{
				nodeColor = node->glThicknessColor;
			}
			else if (nodeOptions.colorization == ColorizationOptions::ByWidth)
			{
				nodeColor = node->glWidthColor;
			}
			else if (nodeOptions.colorization == ColorizationOptions::ByDegree)
			{
				nodeColor = node->glDegreeColor;
			}
			else if (nodeOptions.colorization == ColorizationOptions::ByComponent)
			{
				nodeColor = node->glComponentColor;
			}
			else if (nodeOptions.colorization == ColorizationOptions::Other)
			{
				nodeColor = node->currentColor;
			}
			if ((*mvi.first == selectNode1 && selectNode1Valid) || (*mvi.first == selectNode2 && selectNode2Valid))
			{
				isSelected = true;
			}

			if (degree <= 1)
			{
				scale = nodeOptions.endpointScale;
			}
			else if (degree > 2)
			{
				scale = nodeOptions.junctionScale;
			}

			if (isSelected)
			{
				scale *= nodeSelectionScaling;
			}

			if (degree == 0)
			{
				drawCube.fancierDraw(nodeColor, node->x(), node->y(), node->z(), scale);
			}

			else
			{
				drawSphere.fancierDraw(nodeColor, node->x(), node->y(), node->z(), scale);
			}
			

		}

		//if (selectVertValid)
		//{
		//	nodeColor = selectionColor;
		//	Point3d *p = &mSkeleton[selectVert];
		//	drawSphere.fancierDraw(nodeColor, p->x(), p->y(), p->z(), nodeOptions.endpointScale * nodeSelectionScaling);
		//}
	}

	void BMetaGraph::drawBoxes()
	{
		if (showBoundingBoxes)
		{
			if (onlyDisplaySelectedComponents)
			{
				BoundingBox b1 = componentBounds[selectedComponent1];
				BoundingBox b2 = componentBounds[selectedComponent2];

				std::vector<GLfloat> c1 = ColorTable::getComponentColor(selectedComponent1);
				std::vector<GLfloat> c2 = ColorTable::getComponentColor(selectedComponent2);

				b1.draw(c1, edgeOptions.scale * 2);
				b2.draw(c2, edgeOptions.scale * 2);
			}
			else
			{
				for (int i = 0; i < componentBounds.size(); ++i)
				{
					std::vector<GLfloat> c = ColorTable::getComponentColor(i);

					componentBounds[i].draw(c, edgeOptions.scale * 2);
				}
			}
		}

		//GLfloat boxColor[4] = { 1.0, 0.0, 1.0, 1.0 };
		//for (float x = -10; x < 11; x += 8)
		//{
		//	for (float y = -10; y < 11; y += 5)
		//	{
		//		for (float z = -10; z < 11; z += 3)
		//		{
		//			drawCube.fancierDraw(boxColor, x, y, z, 1.0);
		//		}
		//	}
		//}

	}

	void BMetaGraph::draw()
	{
		drawNodes();
		drawEdges();
		drawBoxes();
	}
}