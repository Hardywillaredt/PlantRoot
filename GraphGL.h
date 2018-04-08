#pragma once
#include "BoostMetaGraph.h"


#include <GLFW/glfw3.h>
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include "boost/python.hpp"

using namespace Roots;



struct Color
{
	double r, g, b, a;
	Color()
		:r(0.0), g(0.0), b(0.0), a(0.0) {}
	Color(double red, double green, double blue, double alpha = 1.0)
		:r(red), g(green), b(blue), a(alpha) {}
};





struct GraphVisualizationOptions
{
	EdgeVisualizationOptions edgeOptions;
	NodeVisualizationOptions nodeOptions;
	Color selectionColor;
};

struct VBOEdges
{
	EdgeVisualizationOptions options;
	std::vector<GLfloat> vertices;
	std::map<ColorizationOptions, std::vector<GLfloat>> colorizations;
	std::vector<GLuint> highlightIndices;
	std::vector<GLuint> baseIndices;


	BMetaGraph &metaGraph;

	VBOEdges(BMetaGraph &graph, );
	void draw();

};

struct VBONodes
{
	NodeVisualizationOptions options;
	
	BMetaGraph &metaGraph;
};

struct MetaGraphModelViewer
{

};