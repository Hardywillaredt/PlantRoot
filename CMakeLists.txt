cmake_minimum_required (VERSION 3.0)


project(RootsToolC)


file(GLOB Roots_src
"*.h"
"*.cpp"
)


set(BOOST_INCLUDEDIR D:/programs/boost_1_65_1/boost_1_65_1/boost_1_65_1/boost CACHE PATH "Boost header directory" FORCE)
set(BOOST_ROOT D:/programs/boost_1_65_1/boost_1_65_1/boost_1_65_1 CACHE PATH "Boost root directory" FORCE)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc")

find_package(Boost REQUIRED)


if (NOT "$ENV{BOOST_INCLUDE_DIRS}" STREQUAL "")
    set(BOOST_INCLUDE_DIRS "$ENV{BOOST_INCLUDE_DIRS}" CACHE STRING "Boost include directory, eg. boost_v_xx_x/boost" FORCE)
endif()
if ("$ENV{BOOST_INCLUDE_DIRS}" STREQUAL "")
	set(BOOST_INCLUDE_DIRS "C:/Programs/boost_1_65_1/boost" CACHE STRING "Boost include directory, eg. boost_v_xx_x/boost")
endif()

if (NOT "$ENV{BOOST_LIBRARY_DIRS}" STREQUAL "")
    set(BOOST_LIBRARY_DIRS "$ENV{BOOST_LIBRARY_DIRS}" CACHE STRING "Boost library directory, eg. boost_v_xx_x/lib64-msvc-14.1" FORCE)
endif()
if ("$ENV{BOOST_LIBRARY_DIRS}" STREQUAL "")
	set(BOOST_LIBRARY_DIRS "C:/Programs/boost_1_65_1/lib64-msvc-14.1" CACHE STRING "Boost library directory, eg. boost_v_xx_x/lib64-msvc-14.1")
endif()


add_subdirectory(SupportLibs)


set(JSON_CPP_LIB ${CMAKE_CURRENT_BINARY_DIR}/SupportLibs/src/lib_json/Debug/jsoncpp.lib CACHE STRING "JSON CPP lib" FORCE)
set(JSON_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/SupportLibs/jsoncpp/include;${CMAKE_CURRENT_SOURCE_DIR}/SupportLibs/jsoncpp/src CACHE STRING "JSON cpp src dirs" FORCE)
message(${JSON_CPP_LIB})
link_libraries(${JSON_CPP_LIB})


file(GLOB BOOST_LIBS "${BOOST_LIBRARY_DIRS}/*.lib")


set(PYTHON_DIR "C:/Users/Will/AppData/Local/Programs/Python/Python35" CACHE PATH "Python includes directory.  Eg. user/AppData/Local/Programs/Python/Python35/include")


add_library(RootsTool SHARED ${Roots_src})
target_include_directories(RootsTool PUBLIC ${JSON_INCLUDES})
target_include_directories(RootsTool PUBLIC ${BOOST_INCLUDE_DIRS})
target_include_directories(RootsTool PUBLIC ${BOOST_INCLUDE_DIRS}/boost)
target_include_directories(RootsTool PUBLIC ${PYTHON_DIR})
target_include_directories(RootsTool PUBLIC ${PYTHON_DIR}/include ${PYTHON_DIR}/libs ${PYTHON_DIR}/DLLs)



# target_link_libraries(RootsTool ${JSON_CPP_LIB})
# target_link_libraries(RootsTool ${BOOST_LIBS})
target_link_libraries(RootsTool ${PYTHON_DIR}/libs/python35.lib)

target_compile_options(
	RootsTool PRIVATE
	"/EHsc")


# add_executable(RootsEXE ${Roots_src})
# target_include_directories(RootsEXE PUBLIC ${JSON_INCLUDES})
# target_include_directories(RootsEXE PUBLIC ${Boost_INCLUDE_DIRS})
# target_link_libraries(RootsEXE ${JSON_CPP_LIB})
# target_link_libraries(RootsEXE ${BOOST_LIBS})