#include "Skeleton.h"
#include "HierarchicalRoot.h"

int main(int argc, char **argv)
{
	//assume if 1 argument that its the filename of the data
	if (argc == 2)
	{
		char * filename = argv[1];
		filename = "D:/2017_Fall/RootsProject/Data/simple.txt";
		
		std::cout << std::endl << std::endl;
		std::cout << filename << std::endl;


		Roots::Skeleton skeleton;
		skeleton.LoadFromTextFile(filename);

		std::cout << "trying my thing" << std::endl;

		std::cout << skeleton;

		//Roots::vertList verts = skeleton.getVertices();
		//

		//Roots::edgePtrList temp = { skeleton.mEdgePtrs[0] };
		//std::cout << "Printing edges " << std::endl;
		//std::cout << skeleton.mEdgePtrs[0]->v0 << " " << skeleton.mEdgePtrs[0]->v1 << std::endl;
		//std::cout << "Done printing edges " << std::endl;

		//Roots::HierarchicalRoot hierarchical = Roots::HierarchicalRoot(temp);
		//skeleton.RemoveEdge(skeleton.mEdgePtrs[0][0]);
		
		int test = 5;

	}
}