#pragma once
#include <GLFW/glfw3.h>
#include "Geometry.h"
#include <vector>
#include <iostream>
namespace Roots
{
	struct Face
	{
		std::vector<GLuint> vertices;
		Point3d center;

		Face(std::vector<GLuint> verts, float *vertexData);
		Face();
		void fancierDraw(GLfloat *color, bool ShowGeom);
		void pickDraw(GLubyte *color, bool showGeom);
		void findCenter(float *vertexData);
		friend std::ostream& operator<<(std::ostream &os, const Face &f);
	};
}