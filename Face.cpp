#include "Face.h"
#include <set>
namespace Roots
{

	Face::Face(std::vector<GLuint> verts, float *vertexData)
	{
		vertices = verts;
		findCenter(vertexData);
	}
	
	Face::Face()
	{
		vertices = {};
		center = Point3d();
	}

	void Face::fancierDraw(GLfloat *color, bool ShowGeom)
	{
		return;
	}
	void Face::pickDraw(GLubyte *color, bool showGeom)
	{
		return;
	}
	void Face::findCenter(float *vertexData)
	{
		center = Point3d();

		std::set<GLuint> uniqueVerts = std::set<GLuint>();

		for (GLuint v : vertices)
		{
			uniqueVerts.insert(v);
		}

		for (std::set<GLuint>::iterator iter = uniqueVerts.begin(); iter != uniqueVerts.end(); ++iter)
		{
			for (int dim = 0; dim < 3; ++dim)
			{
				center[dim] += vertexData[(*iter) * 3 + dim];
			}
		}

		for (int dim = 0; dim < 3; ++dim)
		{
			center[dim] /= uniqueVerts.size();
		}
	}

	std::ostream& operator<<(std::ostream& os, const Face &f)
	{
		return os << 3 << " " << f.vertices[0] << " " << f.vertices[1] << " " << f.vertices[2];
	}

}