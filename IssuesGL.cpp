#include "IssuesGL.h"

#include <GLFW/glfw3.h>
#include <stdlib.h>
#include <stdio.h>


IssuesGL::IssuesGL()
{
	//do nothing
}


void IssuesGL::issueGL()
{
	glBegin(GL_TRIANGLES);
	glColor3f(1.f, 0.f, 0.f);
	glVertex3f(-0.6f, -0.4f, 0.f);
	glColor3f(0.f, 1.f, 0.f);
	glVertex3f(0.6f, -0.4f, 0.f);
	glColor3f(0.f, 0.f, 1.f);
	glVertex3f(0.f, 0.6f, 0.f);
	glEnd();
}