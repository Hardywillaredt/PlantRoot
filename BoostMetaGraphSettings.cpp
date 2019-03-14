#include "BoostMetaGraph.h"

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
		traitSelectionColor[0] = 1.0;
		traitSelectionColor[1] = 1.0;
		traitSelectionColor[2] = 1.0;
		traitSelectionColor[3] = 1.0;
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

	// show loops
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
			operator[](*mei.first).updateColors(edgeOptions, vertexColors, &mSkeleton);
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
			operator[](*mei.first).updateColors(edgeOptions, vertexColors, &mSkeleton);
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
			operator[](*mei.first).updateColors(edgeOptions, vertexColors, &mSkeleton);
		}
	}

	void BMetaGraph::setEdgeSelectionColor(float r, float g, float b)
	{
		edgeOptions.flatSelectionColor[0] = r;
		edgeOptions.flatSelectionColor[1] = g;
		edgeOptions.flatSelectionColor[2] = b;
		return;
	}

	void BMetaGraph::setCurrentPrimaryNodeSelectionColor(float r, float g, float b)
	{
		edgeOptions.traitSelectionColor[0] = r;
		edgeOptions.traitSelectionColor[1] = g;
		edgeOptions.traitSelectionColor[2] = b;
		return;
	}

	void BMetaGraph::showMesh(bool doShow)
	{
		std::cout << "display mesh set to " << doShow << std::endl;
		displayMesh = doShow;
	}

	void BMetaGraph::setMeshAlpha(float alpha)
	{
		std::cout << "Mesh alpha set to " << alpha << std::endl;
		alphaMesh.setAlpha(alpha);
	}
	
	void BMetaGraph::setMeshColor(float red, float green, float blue)
	{
		std::cout << "Mesh color (rgb) set to " << red << " " << green << " " << blue << std::endl;
		alphaMesh.setColor(red, green, blue);
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

	//////////////////////////////////// EXTRA OPTIONS ////////////////////////////////////
	
	/* component options */
	void BMetaGraph::setDisplayOnlySelectedComponents(bool displayOnlySelectedComponents)
	{
		onlyDisplaySelectedComponents = displayOnlySelectedComponents;
		buildEdgeVBOs();
	}

	void BMetaGraph::setComponent1(int component)
	{
		selectedComponent1 = component;
		std::cout << "Component 1 set to " << component << std::endl;
		if (onlyDisplaySelectedComponents)
		{
			buildEdgeVBOs();
		}
	}

	void BMetaGraph::setComponent2(int component)
	{
		selectedComponent2 = component;
		if (selectedComponent2 > numComponents)
		{
			selectedComponent2 = numComponents - 1;
		}
		std::cout << "Component 2 set to " << component << std::endl;
		if (onlyDisplaySelectedComponents)
		{
			buildEdgeVBOs();
		}
	}

	void BMetaGraph::setShowBoundingBoxes(bool doShow)
	{
		showBoundingBoxes = doShow;
		std::cout << "Show bounding boxes set to " << doShow << std::endl;
	}

	/* stem options */
	void BMetaGraph::setDisplayStem(bool doShowStem)
	{
		showStem = doShowStem;
		buildEdgeVBOs();
	}
	
	void BMetaGraph::setDisplayPrimaryNodes(bool doShowPriamryNodes)
	{
		showPrimaryNodes = doShowPriamryNodes;
		buildEdgeVBOs();
	}

	void BMetaGraph::setCurrentPrimaryNode(int node)
	{
		std::cout << "CurrentPrimaryNode " << CurrentPrimaryNode << std::endl;
		if (CurrentPrimaryNode != node)
		{
			CurrentPrimaryNode = node;
			std::cout << "Primary node set to " << node << std::endl;
			unselectAll();
			buildEdgeVBOs();
			CurrentPredecessors = NodesPredecessors[CurrentPrimaryNode].predecessors;
		}
	}

	void BMetaGraph::setDisplayConfirmedPrimaryBranches(bool doShowConfirmedBranches)
	{
		showConfirmedPrimaryBranches = doShowConfirmedBranches;
		buildEdgeVBOs();
	}

	void BMetaGraph::setRandomColorizePrimaryNodes(bool doShow)
	{
		isRandomColorizePrimaryNodes = doShow;
		buildEdgeVBOs();
	}

	void BMetaGraph::setDisplayOnlyBranchesOfCurrentPrimaryNode(bool doShow)
	{
		showOnlyBranchesOfCurrentPrimaryNode = doShow;
		buildEdgeVBOs();
	}

	void BMetaGraph::setDisplayTraitsOnly(bool doShow)
	{
		showTraitsOnly = doShow;
		buildEdgeVBOs();
	}


	void BMetaGraph::setDisplaySelectedSegment(bool doShow)
	{
		showSelectedSegment = doShow;
		buildEdgeVBOs();
	}

	void BMetaGraph::setSegmentHorizontalSliderRadius(float val)
	{
		SegmentHorizontalRadius = val;
		/*
		metaEdgeIter mei = boost::edges(*this);
		for (; mei.first != mei.second; ++mei)
		{
			operator[](*mei.first).updateColors(edgeOptions, vertexColors, &mSkeleton);
		}*/
		SetSelectSegmentPointOperation();
	}
	
	/////////////////////////////////////////DRAWING///////////////////////////////////////////////

	void BMetaGraph::drawEdges()
	{
		if (!edgeOptions.show || !isLoaded)
		{
			return;
		}
		GLfloat *colorArray;
		switch (edgeOptions.colorization)
		{
		case ColorizationOptions::ByThickness:
			colorArray = &vertexColors[0][0];
			break;
		case ColorizationOptions::ByWidth:
			colorArray = &vertexColors[1][0];
			break;
		case ColorizationOptions::ByRatio:
			colorArray = &vertexColors[2][0];
			break;
		case ColorizationOptions::ByComponent:
			colorArray = &vertexColors[3][0];
			break;
		default:
			colorArray = &vertexColors[0][0];
			break;
		}

		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_COLOR_ARRAY);
		glVertexPointer(3, GL_FLOAT, 0, &mSkeleton.glVertices[0]);
		glColorPointer(3, GL_FLOAT, 0, colorArray);

		if (selectionVBO.size() > 0)
		{
			glLineWidth(edgeOptions.scale * edgeSelectionScaling);
			glDrawElements(GL_LINES, selectionVBO.size(), GL_UNSIGNED_INT, &selectionVBO[0]);
		}

		// edges on loop: more thick
		if (edgeOptions.magnifyNonBridges)
		{
			glLineWidth(edgeOptions.scale * nonBridgeScaling);
		}
		else
		{
			glLineWidth(edgeOptions.scale);
		}

		// show loop
		if (nonBridgeVBO.size() > 0 && !showTraitsOnly && !showSelectedSegment)
		{
			glDrawElements(GL_LINES, nonBridgeVBO.size(), GL_UNSIGNED_INT, &nonBridgeVBO[0]);
		}

		// show non-loop
		if (!edgeOptions.showOnlyNonBridges && bridgeVBO.size() > 0 && !showTraitsOnly && !showSelectedSegment)
		{
			glLineWidth(edgeOptions.scale);
			glDrawElements(GL_LINES, bridgeVBO.size(), GL_UNSIGNED_INT, &bridgeVBO[0]);
		}

		if (showTraitsOnly)
		{
			glDrawElements(GL_LINES, stemVBO.size(), GL_UNSIGNED_INT, &stemVBO[0]);
		}

		if (showStem && stemSelected)
		{
			for (std::vector<MetaE>::reverse_iterator riter = StemPath.rbegin(); riter != StemPath.rend(); ++riter)
			{
				MetaV u_tmp = boost::source(*riter, *this);
				MetaV v_tmp = boost::target(*riter, *this);
				MetaE e_tmp = boost::edge(u_tmp, v_tmp, *this).first;
				highLightEdge(e_tmp);
			}
		}

		// only draw selected segment
		if (selectedSegmentVBO.size() > 0 && showSelectedSegment)
		{
			glDrawElements(GL_LINES, selectedSegmentVBO.size(), GL_UNSIGNED_INT, &selectedSegmentVBO[0]);
		}

		// only draw primary branches with color setting
		if (PrimaryBranchesObj.size() > 0)
		{
			if (showTraitsOnly && primaryBranchesVBO.size() > 0)
			{
				glDrawElements(GL_LINES, primaryBranchesVBO.size(), GL_UNSIGNED_INT, &primaryBranchesVBO[0]);
			}

			if (showConfirmedPrimaryBranches && !showOnlyBranchesOfCurrentPrimaryNode)
			{
				for (auto vectorit = PrimaryBranchesObj.begin(); vectorit != PrimaryBranchesObj.end(); ++vectorit)
				{
					//std::cout << vectorit->primaryNode << std::endl;
					//glDrawElements(GL_LINES, vectorit->metaEdges.size(), GL_UNSIGNED_INT, &vectorit->metaEdges[0]);
					for (MetaE edge : vectorit->metaEdges)
					{
						// set confirmed branches the same color as the connected node
						if (isRandomColorizePrimaryNodes)
						{
							int index = vectorit->primaryNodeIndex;
							index = index % 10;
							GLfloat edgeColor[4];
							for (int i = 0; i < 4; ++i)
							{
								edgeColor[i] = randomColorLoopUpTable[index][i];
							}
							colorEdge(edge, edgeColor);
						}
						else
						{
							highLightEdge(edge);
						}
					}
				}
			}
			if (showOnlyBranchesOfCurrentPrimaryNode && !showConfirmedPrimaryBranches)
			{
				int index = 0;
				for (auto vectorit = PrimaryBranchesObj.begin(); vectorit != PrimaryBranchesObj.end(); ++vectorit)
				{
					//glDrawElements(GL_LINES, vectorit->metaEdges.size(), GL_UNSIGNED_INT, &vectorit->metaEdges[0]);
					//std::cout << vectorit->primaryNode << std::endl;
					if (vectorit->primaryNodeIndex == CurrentPrimaryNode)
					{
						GLfloat edgeColor[4];
						for (int i = 0; i < 4; ++i)
						{
							edgeColor[i] = randomColorLoopUpTable[index][i];
						}
						for (MetaE edge : vectorit->metaEdges)
						{
							colorEdge(edge, edgeColor);
						}
						++index;
						index = index % 10;
					}
				}
			}
			if (showOnlyBranchesOfCurrentPrimaryNode && showConfirmedPrimaryBranches)
			{
				int index = 0;
				for (auto vectorit = PrimaryBranchesObj.begin(); vectorit != PrimaryBranchesObj.end(); ++vectorit)
				{
					//glDrawElements(GL_LINES, vectorit->metaEdges.size(), GL_UNSIGNED_INT, &vectorit->metaEdges[0]);
					// show branches of current node using random color
					if (vectorit->primaryNodeIndex == CurrentPrimaryNode)
					{
						GLfloat edgeColor[4];
						for (int i = 0; i < 4; ++i)
						{
							edgeColor[i] = randomColorLoopUpTable[index][i];
						}
						for (MetaE edge : vectorit->metaEdges)
						{
							colorEdge(edge, edgeColor);
						}
						++index;
						index = index % 10;
					}
					// show branches of other nodes using the same color as their connected nodes
					else
					{
						for (MetaE edge : vectorit->metaEdges)
						{
							if (isRandomColorizePrimaryNodes)
							{
								int index2 = vectorit->primaryNodeIndex;
								index2 = index2 % 10;
								GLfloat edgeColor[4];
								for (int i = 0; i < 4; ++i)
								{
									edgeColor[i] = randomColorLoopUpTable[index2][i];
								}
								colorEdge(edge, edgeColor);
							}
							else
							{
								highLightEdge(edge);
							}
						}
					}
				}
			}
		}

		glDisableClientState(GL_VERTEX_ARRAY);
		glDisableClientState(GL_COLOR_ARRAY);
	}

	void BMetaGraph::edgePickRender()
	{
		glDisable(GL_LIGHTING);
		glClearColor(1.0, 1.0, 1.0, 1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
		glLoadMatrixf(projectionTransform);

		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glLoadMatrixf(modelViewTransform);

		if (!edgeOptions.show || !isLoaded)
		{
			return;
		}
		GLubyte edgeIdColor[3];

		glEnableClientState(GL_VERTEX_ARRAY);
		glVertexPointer(3, GL_FLOAT, 0, &mSkeleton.glVertices[0]);


		metaEdgeIter mei = boost::edges(*this);

		int id = 0;
		for (; mei.first != mei.second; ++mei, ++id)
		{
			BMetaEdge *e = &operator[](*mei.first);
			e->instanceId = id;
			int updatedId = operator[](*mei.first).instanceId;
			int desiredId = id;
			if (updatedId != desiredId)
			{
				std::cout << "Updated and desired Id's do not match " << std::endl;
				std::cout << "Desired : " << desiredId << " updated : " << updatedId << std::endl;
			}
			

			edgeIdColor[0] = (e->instanceId & 0x000000FF) >> 0;
			edgeIdColor[1] = (e->instanceId & 0x0000FF00) >> 8;
			edgeIdColor[2] = (e->instanceId & 0x00FF0000) >> 16;
			glColor3ub(edgeIdColor[0], edgeIdColor[1], edgeIdColor[2]);

			if (edgeOptions.showOnlyNonBridges && e->isBridge)
			{
				continue;
			}

			if (e->isSelected)
			{
				glLineWidth(edgeOptions.scale * edgeSelectionScaling);
			}
			else if (edgeOptions.magnifyNonBridges && !e->isBridge)
			{
				glLineWidth(edgeOptions.scale * nonBridgeScaling);
			}
			else
			{
				glLineWidth(edgeOptions.scale);
			}
			

			glDrawElements(GL_LINES, e->indicesList.size(), GL_UNSIGNED_INT, &e->indicesList[0]);
		}

		glDisableClientState(GL_VERTEX_ARRAY);


		glPopMatrix();
		glMatrixMode(GL_PROJECTION);
		glPopMatrix();

	}

	void BMetaGraph::drawNodes()
	{
		
		metaVertIter mvi = boost::vertices(*this);
		
		GLfloat *nodeColor;
		for (; mvi.first != mvi.second; ++mvi)
		{
			BMetaNode *node = &operator[](*mvi.first);
			int degree = boost::degree(*mvi.first, *this);

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
			
			//if (!PrimaryBranchesNode.empty() && (std::find(PrimaryBranchesNode.begin(), PrimaryBranchesNode.end(), *mvi.first) != PrimaryBranchesNode.end()))
			if (PrimaryBranchSelectionValid && (*mvi.first == PrimaryBranchSelection))
			{
				isSelected = true;
				//nodeColor = selectedPrimaryNodeColor;
			}

			if ((*mvi.first == selectNode1 && selectNode1Valid) || (*mvi.first == selectNode2 && selectNode2Valid))
			{
				isSelected = true;
			}

			if ((*mvi.first == selectStemStart && selectStemStartValid) || (*mvi.first == selectStemEnd && selectStemEndValid))
			{
				isSelected = true;
			}

			// primary nodes setting
			MetaV temp = *mvi.first;
			auto it = std::find_if(PrimaryNodes.begin(), PrimaryNodes.end(), [&temp](const std::pair<int, MetaV>& element) { return element.second == temp; });
			if ( it != PrimaryNodes.end() )
			{
				//isSelected = true; // all primary nodes : show bigger
				if (isRandomColorizePrimaryNodes)
				{
					int index = std::distance(PrimaryNodes.begin(), it);
					index = index % 10;
					GLfloat temp[4];
					for (int i = 0; i < 4; ++i)
					{
						temp[i] = randomColorLoopUpTable[index][i];
					}
					nodeColor = temp;
				}
				//if (selectPrimaryNodesValid)//&& (temp == PrimaryNodes.back().second))
				if (showPrimaryNodes)
				{
					isSelected = true;
				}
			}

			if ((CurrentPrimaryNode != -1) && (it != PrimaryNodes.end()))
			{
				if (*mvi.first == PrimaryNodes[CurrentPrimaryNode].second)
				{
					isSelected = true;
				}
			}

			if ((*mvi.first == selectSegmentPoint1 && selectSegmentPoint1Valid) || (*mvi.first == selectSegmentPoint2 && selectSegmentPoint2Valid))
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
				// std::cout << "drawNode degree = 0" << std::endl;
			}
			else
			{
				drawSphere.fancierDraw(nodeColor, node->x(), node->y(), node->z(), scale);
			}
		}
	}

	void BMetaGraph::nodePickRender()
	{
		glDisable(GL_LIGHTING);
		glClearColor(1.0, 1.0, 1.0, 1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
		glLoadMatrixf(projectionTransform);

		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glLoadMatrixf(modelViewTransform);

		metaVertIter mvi = boost::vertices(*this);
		std::cout << "Node id's " << std::endl;
		GLubyte nodeIdColor[3];
		for (; mvi.first != mvi.second; ++mvi)
		{
			BMetaNode *node = &operator[](*mvi.first);
			int degree = boost::degree(*mvi.first, *this);

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

			int num = (*mvi.first);
			nodeIdColor[0] = (num & 0x000000FF) >> 0;
			nodeIdColor[1] = (num & 0x0000FF00) >> 8;
			nodeIdColor[2] = (num & 0x00FF0000) >> 16;

			if ((*mvi.first == selectNode1 && selectNode1Valid) || (*mvi.first == selectNode2 && selectNode2Valid))
			{
				isSelected = true;
			}
			
			if ((*mvi.first == selectStemStart && selectStemStartValid) || (*mvi.first == selectStemEnd && selectStemEndValid))
			{
				isSelected = true;
			}

			MetaV temp = *mvi.first;
			auto it = std::find_if(PrimaryNodes.begin(), PrimaryNodes.end(), [&temp](const std::pair<int, MetaV>& element) { return element.second == temp; });
			if (it != PrimaryNodes.end())
			{
				isSelected = true;
				//std::cout << "nodepickrender pri node" << std::endl;
			}

			//if (!PrimaryBranchesNode.empty() && (std::find(PrimaryBranchesNode.begin(), PrimaryBranchesNode.end(), *mvi.first) != PrimaryBranchesNode.end()))
			if (PrimaryBranchSelectionValid && (*mvi.first == PrimaryBranchSelection))
			{
				isSelected = true;
			}

			if ((*mvi.first == selectSegmentPoint1 && selectSegmentPoint1Valid) || (*mvi.first == selectSegmentPoint2 && selectSegmentPoint2Valid))
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
				if (!PrimaryBranchSelectionValid)
				{
					drawCube.pickDraw(nodeIdColor, node->x(), node->y(), node->z(), scale);
					std::cout << "nodePickRender degree = 0" << std::endl;
				}
			}
			else
			{
				drawSphere.pickDraw(nodeIdColor, node->x(), node->y(), node->z(), scale);
				//std::cout << "test draw node 2333 " << std::endl;
			}
		}

		glPopMatrix();
		glMatrixMode(GL_PROJECTION);
		glPopMatrix();
	}


	void BMetaGraph::vertPickRender()
	{
		glDisable(GL_LIGHTING);
		glClearColor(1.0, 1.0, 1.0, 1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
		glLoadMatrixf(projectionTransform);

		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glLoadMatrixf(modelViewTransform);


		metaEdgeIter mei = boost::edges(*this);

		GLubyte vertIdColor[3];
		
		float scale = 0.5 * (nodeOptions.endpointScale + nodeOptions.junctionScale);
		for (; mei.first != mei.second; ++mei)
		{
			BMetaEdge *e = &operator[](*mei.first);
			if (onlyDisplaySelectedComponents)
			{
				if (e->connectedComponent != selectedComponent1 && e->connectedComponent != selectedComponent2)
				{
					continue;
				}
			}

			//we only need to render the middle vertices (not endpoints eg. metanodes) because we have already passed over those in the picking loop
			for (int i = 1; i < e->mVertices.size() - 1; ++i)
			{
				SkelVert v = e->mVertices[i];
				vertIdColor[0] = (v & 0x000000FF) >> 0;
				vertIdColor[1] = (v & 0x0000FF00) >> 8;
				vertIdColor[2] = (v & 0x00FF0000) >> 16;
				Point3d *p = &mSkeleton[v];
				drawSphere.pickDraw(vertIdColor, p->x(), p->y(), p->z(), scale);
			}
		}

		glPopMatrix();
		glMatrixMode(GL_PROJECTION);
		glPopMatrix();
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

	}

	

	void BMetaGraph::draw()
	{
		if (useArcball)
		{
			//glMatrixMode(GL_MODELVIEW);
			glTranslatef(eyeShiftX, eyeShiftY, -mSkeleton.mRadius * 2);
			arcball_rotate();
			glTranslatef(-viewCenter.x(), -viewCenter.y(), -viewCenter.z());
		}

		/*
		if (isDisplaySelectedSegment())
		{
			//getClipPlane();
			glEnable(GL_CLIP_PLANE0);
			glEnable(GL_CLIP_PLANE1);
			GLdouble *c = new GLdouble[4];
			for (int i = 0; i < 4; ++i)
			{
				c[i] = rootClipPlaneNormal[0][i];
			}
			glClipPlane(GL_CLIP_PLANE0, c);
			for (int i = 0; i < 4; ++i)
			{
				c[i] = rootClipPlaneNormal[1][i];
			}
			glClipPlane(GL_CLIP_PLANE1, c);
		}*/

		glEnable(GL_LIGHTING);
		drawNodes();
		glDisable(GL_LIGHTING);
		drawEdges();
		drawBoxes();
		glGetFloatv(GL_MODELVIEW_MATRIX, modelViewTransform);
		glGetFloatv(GL_PROJECTION_MATRIX, projectionTransform);
		/*
		if (isDisplaySelectedSegment())
		{
			glDisable(GL_CLIP_PLANE0);
			glDisable(GL_CLIP_PLANE1);
		}*/

		if (displayMesh)
		{
			alphaMesh.render();
		}
		glEnable(GL_LIGHTING);
	}
}