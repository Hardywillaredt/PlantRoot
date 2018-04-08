#include "Skeleton.h"
#include "HierarchicalRoot.h"
#include "PythonModule.h"
#include "BoostMetaGraph.h"

//int main(int argc, char **argv)
//{
//	//assume if 1 argument that its the filename of the data
//	if (argc == 2)
//	{
//		char * filename = argv[1];
//		std::string fileString = std::string(filename);
//
//		Roots::BMetaGraph mgraph = Roots::BMetaGraph(fileString);
//		Roots::MetaEdge3d edge1 = Roots::MetaEdge3d();
//		std::vector<Roots::MetaEdge3d> edges = std::vector<Roots::MetaEdge3d>(2);
//		edge1.node0 = 0;
//		edge1.node1 = 1;
//
//		edges[0].node0 = 0;
//		edges[0].node1 = 2;
//		edges[1].node0 = 1;
//		edges[1].node1 = 5;
//
//		mgraph.SplitOperation(edge1, edges);
//		//Roots::Skeleton skel = Roots::Skeleton(fileString);
//		
//		//std::cout << std::endl << std::endl;
//		//std::cout << filename << std::endl;
//
//
//		//Roots::Skeleton skeleton;
//		//skeleton.LoadFromTextFile(filename);
//
//		//std::cout << "trying my thing" << std::endl;
//
//		//std::cout << skeleton;
//
//		//Roots::vertList verts = skeleton.getVertices();
//		//
//
//		//Roots::edgePtrList temp = { skeleton.mEdgePtrs[0] };
//		//std::cout << "Printing edges " << std::endl;
//		//std::cout << skeleton.mEdgePtrs[0]->v0 << " " << skeleton.mEdgePtrs[0]->v1 << std::endl;
//		//std::cout << "Done printing edges " << std::endl;
//
//		//Roots::HierarchicalRoot hierarchical = Roots::HierarchicalRoot(temp);
//		//skeleton.RemoveEdge(skeleton.mEdgePtrs[0][0]);
//		
//		int test = 5;
//		return 0;
//	}
//}

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
int main(void)
{
	GLFWwindow* window;
	glfwSetErrorCallback(error_callback);
	if (!glfwInit())
		exit(EXIT_FAILURE);
	window = glfwCreateWindow(640, 480, "Simple example", NULL, NULL);
	if (!window)
	{
		glfwTerminate();
		exit(EXIT_FAILURE);
	}
	glfwMakeContextCurrent(window);
	glfwSwapInterval(1);
	glfwSetKeyCallback(window, key_callback);
	drawing::VBOSphere sphere = drawing::VBOSphere(1.0, 0.0, 1.0, 1.0, 0.5, 3);
	GLfloat color[4] = { 1.0, 0.0, 1.0, 1.0 };
	GLfloat lightPos[3] = { 10, 10, 10 };
	while (!glfwWindowShouldClose(window))
	{
		float ratio;
		int width, height;
		glfwGetFramebufferSize(window, &width, &height);
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);

		glLightfv(GL_LIGHT0, GL_POSITION, lightPos);

		ratio = width / (float)height;
		glViewport(0, 0, width, height);
		glClear(GL_COLOR_BUFFER_BIT);
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glOrtho(-ratio, ratio, -1.f, 1.f, 1.f, 10000.f);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, color);
		glColor4f(color[0], color[1], color[2], color[3]);
		sphere.draw();
		//glRotatef((float)glfwGetTime() * 50.f, 0.f, 0.f, 1.f);
		//glBegin(GL_TRIANGLES);
		//glColor3f(1.f, 0.f, 0.f);
		//glVertex3f(-0.6f, -0.4f, 0.f);
		//glColor3f(0.f, 1.f, 0.f);
		//glVertex3f(0.6f, -0.4f, 0.f);
		//glColor3f(0.f, 0.f, 1.f);
		//glVertex3f(0.f, 0.6f, 0.f);
		//glEnd();
		glfwSwapBuffers(window);
		glfwPollEvents();
	}
	glfwDestroyWindow(window);
	glfwTerminate();
	exit(EXIT_SUCCESS);
}