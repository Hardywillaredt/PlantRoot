#include "PythonModule.h"
#include "BoostMetaGraph.h"
#include <GLFW/glfw3.h>
#include <stdlib.h>
#include <stdio.h>
#include "Sphere.h"

static void error_callback(int error, const char* description)
{
	fputs(description, stderr);
}
static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
		glfwSetWindowShouldClose(window, GL_TRUE);
}
int main(int argc, char **argv)
{
	if (argc == 2)
	{
		char * filename = argv[1];
		std::string fileString = std::string(filename);

		Roots::BMetaGraph mgraph = Roots::BMetaGraph();
		mgraph.loadFromFile(filename);
		metaEdgeIter mei = boost::edges(mgraph);
		for (; mei.first != mei.second; ++mei)
		{
			MetaE e = *mei.first;
			if (boost::degree(e.m_source, mgraph) >= 3 && boost::degree(e.m_target, mgraph) >= 3)
			{
				MetaV v0 = e.m_source;
				MetaV v1 = e.m_target;
				BMetaGraph::adjacency_iterator adjIt, adjEnd;
				boost::tie(adjIt, adjEnd) = boost::adjacent_vertices(v0, mgraph);
				MetaE e0;
				bool exists = false;
				for (; adjIt != adjEnd; ++adjIt)
				{
					boost::tie(e0, exists) = boost::edge(*adjIt, v0, mgraph);
					if (e0 != e)
					{
						break;
					}
				}
				MetaE e1;
				boost::tie(adjIt, adjEnd) = boost::adjacent_vertices(v1, mgraph);
				for (; adjIt != adjEnd; ++adjIt)
				{
					boost::tie(e1, exists) = boost::edge(*adjIt, v1, mgraph);
					if (e1 != e)
					{
						break;
					}
				}

				mgraph.splitEdge = e;
				mgraph.splitNeighbors = { e0, e1 };
				mgraph.splitEdgeValid = true;
				mgraph.SplitOperation();
			}
		}

		mei = boost::edges(mgraph);
		MetaE e = *mei.first;
		mgraph.removeComponentEdge = e;
		mgraph.removeComponentEdgeValid = true;
		mgraph.RemoveComponentOperation();
	}
	

}