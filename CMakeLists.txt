cmake_minimum_required (VERSION 3.0)


project(RootsToolC)
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../python CACHE PATH "python runtime directory")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

file(GLOB Roots_src
"*.h"
"*.cpp"
)


# set(BOOST_INCLUDEDIR D:/programs/boost_1_65_1/boost_1_65_1/boost_1_65_1/boost CACHE PATH "Boost header directory" FORCE)
# set(BOOST_ROOT D:/programs/boost_1_65_1/boost_1_65_1/boost_1_65_1 CACHE PATH "Boost root directory" FORCE)

# set(GLFW_INCLUDEDIR D:/dev/programs/glfw-3.2.1/glfw-3.2.1/include CACHE PATH "GLFW include directory")
# set(GLFW_LIB D:/dev/programs/glfw-3.2.1/glfw-3.2.1/glfw-build/src/Debug/glfw3.lib CACHE STRING "glfw compiled library")
# set(GLAD_INCLUDEDIR D:/dev/programs/glad/include CACHE PATH "GLAD include directory")

# set(BoostPythonLib D:/dev/school/rootsProject/SupportLibs/boostPython.lib CACHE STRING "boost python library")

# set(OpenCV_DIR D:/dev/programs/opencv/build CACHE PATH "Opencv Include")
# set(OpenCV_VERSION "3.4.2" CACHE STRING "OpenCV installation version")
# option(OPENCV_ENABLE "description" ON)

set(BOOST_INCLUDEDIR E:/will/boost_1_65_1/boost CACHE PATH "Boost header directory" FORCE)
set(BOOST_ROOT E:/will/boost_1_65_1 CACHE PATH "Boost root directory" FORCE)

set(GLFW_INCLUDEDIR E:/will/glfw-3.2.1.bin.WIN64/include CACHE PATH "GLFW include directory")
set(GLFW_LIB E:/will/glfw-3.2.1.bin.WIN64/lib-vc2015/glfw3.lib CACHE STRING "glfw compiled library")
set(GLAD_INCLUDEDIR E:/will/glad/include CACHE PATH "GLAD include directory")

set(BoostPythonLib E:/will/boost_1_65_1/stage/lib/libboost_python3-vc141-mt-gd-1_65_1.lib CACHE STRING "boost python library")

# set(OpenCV_DIR E:/will/opencv-3.4.2/build CACHE PATH "Opencv Include")
# set(OpenCV_VERSION "3.4.2" CACHE STRING "OpenCV installation version")
option(OPENCV_ENABLE "description" OFF)

if(${OPENCV_ENABLE})
	find_package(OpenCV REQUIRED)
endif()

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc")

find_package(Boost REQUIRED)

include_directories(E:/will/boost_1_65_1)
if (NOT "$ENV{BOOST_INCLUDE_DIRS}" STREQUAL "")
    set(BOOST_INCLUDE_DIRS "$ENV{BOOST_INCLUDE_DIRS}" CACHE STRING "Boost include directory, eg. boost_v_xx_x/boost" FORCE)
endif()
if ("$ENV{BOOST_INCLUDE_DIRS}" STREQUAL "")
	set(BOOST_INCLUDE_DIRS "E:/will/boost_1_65_1" CACHE STRING "Boost include directory, eg. boost_v_xx_x/boost")
endif()

if (NOT "$ENV{BOOST_LIBRARY_DIRS}" STREQUAL "")
    set(BOOST_LIBRARY_DIRS "$ENV{BOOST_LIBRARY_DIRS}" CACHE STRING "Boost library directory, eg. boost_v_xx_x/lib64-msvc-14.1" FORCE)
endif()
if ("$ENV{BOOST_LIBRARY_DIRS}" STREQUAL "")
	set(BOOST_LIBRARY_DIRS "E:/will/boost_1_65_1/stage/lib" CACHE STRING "Boost library directory, eg. boost_v_xx_x/lib64-msvc-14.1")
endif()

set(OUTPUT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../python" CACHE PATH "Output directory for binaries/dlls, typically project python directory")

foreach( OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES} )
    string( TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG )
    set( CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${OUTPUT_DIR} )
    set( CMAKE_LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${OUTPUT_DIR} )
    set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${OUTPUT_DIR} )
endforeach( OUTPUTCONFIG CMAKE_CONFIGURATION_TYPES )

#add_subdirectory(SupportLibs)


#set(JSON_CPP_LIB ${CMAKE_CURRENT_BINARY_DIR}/SupportLibs/src/lib_json/Debug/jsoncpp.lib CACHE STRING "JSON CPP lib" FORCE)
#set(JSON_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/SupportLibs/jsoncpp/include;${CMAKE_CURRENT_SOURCE_DIR}/SupportLibs/jsoncpp/src CACHE STRING "JSON cpp src dirs" FORCE)
#message(${JSON_CPP_LIB})
#link_libraries(${JSON_CPP_LIB})


file(GLOB BOOST_LIBS "${BOOST_LIBRARY_DIRS}/*.lib")
link_directories(${BOOST_LIBRARY_DIRS})
# link_directories(${OpenCV_DIR}/lib/Debug)

set(PYTHON_DIR "C:/Python35" CACHE PATH "Python includes directory.  Eg. user/AppData/Local/Programs/Python/Python35/include")

add_library(RootsTool SHARED ${Roots_src})
#target_include_directories(RootsTool PUBLIC ${JSON_INCLUDES})
target_include_directories(RootsTool PUBLIC ${BOOST_INCLUDE_DIRS})
target_include_directories(RootsTool PUBLIC ${BOOST_INCLUDE_DIRS}/boost)
target_include_directories(RootsTool PUBLIC ${PYTHON_DIR})
target_include_directories(RootsTool PUBLIC ${PYTHON_DIR}/include ${PYTHON_DIR}/libs ${PYTHON_DIR}/DLLs)
target_include_directories(RootsTool PUBLIC ${GLFW_INCLUDEDIR})
target_include_directories(RootsTool PUBLIC ${GLAD_INCLUDEDIR})

#comment this line for debug version
target_compile_definitions(RootsTool PRIVATE BOOST_PYTHON_STATIC_LIB)
target_link_libraries(RootsTool ${BoostPythonLib})

if(${OPENCV_ENABLE})
	# include_directories(E:/will/opencv-3.4.2/include)
	# target_include_directories(RootsTool PUBLIC ${OpenCV_DIR}/lib/Debug ${OpenCV_DIR}/install/include)
	target_include_directories(RootsTool PUBLIC ${OpenCV_DIR}/lib/Debug ${OpenCV_DIR}/include)
	target_include_directories(RootsTool PUBLIC E:/will/opencv/build/include)
	# target_link_libraries(RootsTool opencv_calib3d342d.lib opencv_core342d.lib opencv_dnn342d.lib opencv_features2d342d.lib opencv_flann342d.lib opencv_highgui342d.lib opencv_imgcodecs342d.lib opencv_imgproc342d.lib opencv_ml342d.lib opencv_objdetect342d.lib opencv_photo342d.lib opencv_shape342d.lib opencv_stitching342d.lib opencv_superres342d.lib opencv_ts342d.lib opencv_video342d.lib opencv_videoio342d.lib opencv_videostab342d.lib)
	target_link_libraries(RootsTool opencv_calib3d342d.lib opencv_core342d.lib opencv_features2d342d.lib opencv_flann342d.lib opencv_highgui342d.lib opencv_imgcodecs342d.lib opencv_imgproc342d.lib opencv_ts342d.lib opencv_videoio342d.lib)
	target_compile_definitions(RootsTool PRIVATE WITH_OPENCV=TRUE)
endif()
target_link_libraries(RootsTool ${PYTHON_DIR}/libs/python35.lib ${GLFW_LIB} opengl32.lib glu32.lib)





add_executable(RootsEXE ${Roots_src})
#target_include_directories(RootsEXE PUBLIC ${JSON_INCLUDES})
target_include_directories(RootsEXE PUBLIC ${BOOST_INCLUDE_DIRS})
target_include_directories(RootsEXE PUBLIC ${BOOST_INCLUDE_DIRS}/boost)
target_include_directories(RootsEXE PUBLIC ${PYTHON_DIR})
target_include_directories(RootsEXE PUBLIC ${PYTHON_DIR}/include ${PYTHON_DIR}/libs ${PYTHON_DIR}/DLLs)
target_include_directories(RootsEXE PUBLIC ${GLFW_INCLUDEDIR})
target_include_directories(RootsEXE PUBLIC ${GLAD_INCLUDEDIR})

#comment this line for debug version
target_compile_definitions(RootsEXE PRIVATE BOOST_PYTHON_STATIC_LIB)
target_link_libraries(RootsEXE ${BoostPythonLib})


target_link_libraries(RootsEXE ${PYTHON_DIR}/libs/python35.lib ${GLFW_LIB} opengl32.lib glu32.lib)

if(${OPENCV_ENABLE})
	target_include_directories(RootsEXE PUBLIC ${OpenCV_DIR}/lib/Debug ${OpenCV_DIR}/include)
	target_include_directories(RootsTool PUBLIC E:/will/opencv/build/include)
	# target_link_libraries(RootsEXE opencv_calib3d342d.lib opencv_core342d.lib opencv_dnn342d.lib opencv_features2d342d.lib opencv_flann342d.lib opencv_highgui342d.lib opencv_imgcodecs342d.lib opencv_imgproc342d.lib opencv_ml342d.lib opencv_objdetect342d.lib opencv_photo342d.lib opencv_shape342d.lib opencv_stitching342d.lib opencv_superres342d.lib opencv_ts342d.lib opencv_video342d.lib opencv_videoio342d.lib opencv_videostab342d.lib)
	target_link_libraries(RootsEXE opencv_calib3d342d.lib opencv_core342d.lib opencv_features2d342d.lib opencv_flann342d.lib opencv_highgui342d.lib opencv_imgcodecs342d.lib opencv_imgproc342d.lib opencv_ts342d.lib opencv_videoio342d.lib)
	target_compile_definitions(RootsEXE PRIVATE WITH_OPENCV=TRUE)
endif()

set_target_properties(RootsEXE PROPERTIES RUNTIME_OUTPUT_DIRECTORY "E:/PlantRoot-master/x64/Debug")


target_compile_options(
	RootsTool PRIVATE
	"/EHsc")

#add_executable(RootsEXE ${Roots_src})
#target_include_directories(RootsEXE PUBLIC ${JSON_INCLUDES})
#target_include_directories(RootsEXE PUBLIC ${BOOST_INCLUDE_DIRS})
#target_include_directories(RootsEXE PUBLIC ${BOOST_INCLUDE_DIRS}/boost)
#target_include_directories(RootsEXE PUBLIC ${PYTHON_DIR})
#target_include_directories(RootsTool PUBLIC ${PYTHON_DIR}/include ${PYTHON_DIR}/libs ${PYTHON_DIR}/DLLs)
#target_link_libraries(RootsEXE ${PYTHON_DIR}/libs/python35.lib)

#target_compile_options(
#	RootsEXE PRIVATE
#	"/EHsc")

# target_link_libraries(RootsTool ${JSON_CPP_LIB})
# target_link_libraries(RootsTool ${BOOST_LIBS})




set_target_properties(RootsTool PROPERTIES SUFFIX ".pyd")
#								ARCHIVE_OUTPUT_DIRECTORY "${OUTPUT_DIR}"
#								LIBRARY_OUTPUT_DIRECTORY "${OUTPUT_DIR}"
#								RUNTIME_OUTPUT_DIRECTORY "${OUTPUT_DIR}"
#								)
# add_executable(RootsEXE ${Roots_src})
# target_include_directories(RootsEXE PUBLIC ${JSON_INCLUDES})
# target_include_directories(RootsEXE PUBLIC ${Boost_INCLUDE_DIRS})
# target_link_libraries(RootsEXE ${JSON_CPP_LIB})
# target_link_libraries(RootsEXE ${BOOST_LIBS})