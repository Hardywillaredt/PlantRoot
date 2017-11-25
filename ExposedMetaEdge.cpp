#include "MetaGraph.h"
#include "boost/python.hpp"
//#include "boost\python\suite\indexing\map_indexing_suite.hpp"
//#include "boost\python\suite\indexing\vector_indexing_suite.hpp"

/////////////////////Vector Typedefs//////////////////////
typedef std::vector<double> dVec;
typedef std::vector<std::string> sVec;
typedef std::vector<int> iVec;


/////////////////////EXPOSURE CLASSES//////////////////////
//
//struct TestClass
//{
//	dVec GetTestVecD()
//	{
//		return myDVec;
//	}
//	void SetTestVecD(const dVec& dVec)
//	{
//		myDVec = dVec;
//	}
//
//private:
//	dVec myDVec;
//};
//
///////////////////////WRAPPER CODE//////////////////////////
//
//using namespace boost::python;
//
//BOOST_PYTHON_MODULE(Roots)
//{
//	class_<dVec>("dVec")
//		.def(vector_indexing_suite<dVec>());
//	class_<sVec>("sVec")
//		.def(vector_indexing_suite<sVec>());
//	class_<iVec>("iVec")
//		.def(vector_indexing_suite<iVec>());
//
//	class_<TestClass>("TestClass")
//		.def("GetTestVecD", &TestClass::GetTestVecD)
//		.def("SetTestVecD", &TestClass::SetTestVecD)
//		;
//};
//
//
//struct ExposedMetaGraph
//{
//
//};


//struct ExposedMetaGraph
//{
//
//};

