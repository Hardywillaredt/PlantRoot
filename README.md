# PlantRoot
### Aim ###
Discover biologically significant traits form CT scan of plant root.
### How ###
1. Using [persistent homology](https://www.cse.wustl.edu/~taoju/research/homology_CCCG.pdf) to deal with topological noise.
2. Build graph structure for skeleton.
3. Build root heirarchy to understand biologically meaningful traits.
### Skeleton Data Structure ###
1. Lower level
   - node - defined by 3D location (x,y,z)
   - edge - defined by 2 nodes
2. Higher level - visible in GUI
   - metaNode - the intersection of three or more edges
   - metaEdge - consist of a set of edges between two metaNode
   - metaFace


## Installation Guide ##
 - Environment: Windows 10 64bit, Python 3.7.3 64bit, Visual Studio 2015 or 2017
 - Download python from: https://www.python.org/ftp/python/3.7.3/python-3.7.3-amd64.exe
 - After installation, if you cannot run python in command line, right click on "This PC", select "Properties -- Advanced System Settings -- Environment Variables.." and add your directory which contains python.exe, then you should be able to run python and pip3.
 - Note that Visual Studio 2015 does NOT install the languages by default. If you are using Visual Studio 2015, you need to manually choose to install the languages in the customized section during setup.
#### Install Dependencies ####
1. Download Boost 1.68.0 from: https://dl.bintray.com/boostorg/release/1.68.0/source/boost_1_68_0.zip and extract the contents to a folder of your choice.
2. In command line, cd to its directory and run ".\bootstrap.bat". You should see an executable file "b2" is generated.
3. In command line, run ".\b2 --address-mode=64 --with-python" to install libboost python library. After a successful installation, the .lib file should be in boost_1_68_0/stage/lib.
4. Download glfw-3.3 64bit from: https://github.com/glfw/glfw/releases/download/3.3/glfw-3.3.bin.WIN64.zip and extract the contents to a folder of your choice.
5. Generate glad from: http://glad.dav1d.de/#profile=compatibility&specification=gl&api=gl%3D4.6&api=gles1%3D1.0&api=gles2%3D3.2&api=glsc2%3D2.0&language=c&loader=on Language: C/C++, Specification: OpenGL, Profile: Compatibility, API: Select all latest versions (gl: 4.6, gles1: 1.0, gles2: 3.2, glsc2: 2.0), Extensions: Leave blank, Options: Generate a loader. Click generate and extract the contents to a folder of your choice.
6. Download OpenCV 4.1.0 from: https://sourceforge.net/projects/opencvlibrary/files/4.1.0/opencv-4.1.0-vc14_vc15.exe/download and install it to a folder of your choice.
Having the same version as listed above guarentees a working environment. Having other versions might work but is not tested.
#### Known conflicts and bugs ####
1. Boost 1.65.0 doesn't work with python 3.7.
2. Boost 1.67.0 has a bug in python/details/config.hpp that makes the program unable to find libboost_python library file.
3. Boost 1.70 doesn't respond when running bootstrap.
4. Python 32bit version doesn't work with libboost_python 64bit version. Make sure you download both as 64bit when you have a 64bit windows.
#### Modify CMakeLists.txt #### 
1. Navigate to your project directory and find CMakeLists.txt.
2. Modify the following lines:
   - set(BOOST_INCLUDEDIR ... <- This is your boost header (.../boost_1_68_0/boost)
   - set(BOOST_ROOT ... <- This is your boost directory (.../boost_1_68_0)
   - set(GLFW_INCLUDEDIR ... <- This is your glfw include directory (.../glfw-3.3.bin.WIN64/include)
   - set(GLFW_LIB ... <- This is your glfw library file (.../glfw-3.3.bin.WIN64/lib-vc2015/glfw3.lib)
   - set(GLAD_INCLUDEDIR ... <- This is your glad include directory (.../glad/include)
   - set(BoostPythonLib ... <- This is your boost python library (.../boost_1_68_0/stage/lib/libboost_python37-vc140-mt-gd-x64-1_68.lib) Make sure you use "mt-gd-x64" if there are other ones.
   - set(OpenCV_DIR ... <- This is your opencv build directory (.../opencv/build)
   - set(BOOST_INCLUDE_DIRS ... <- This is the same as your boost directory (.../boost_1_68_0)
   - set(BOOST_LIBRARY_DIRS ... <- This is your boost python libary directory (.../boost_1_68_0/stage/lib)
   - set(PYTHON_DIR ... <- This is the directory that has your python.exe (.../Python/Python37)
   - set_target_properties(RootsEXE PROPERTIES RUNTIME_OUTPUT_DIRECTORY ... <- This is the debug directory (.../PlantRoot/x64/Debug)
3. Save the file and in command line, cd to your project directory and run 'cmake -G "Visual Studio Win64"'.
4. Now you should see "RootsToolC.sln" solution generated in your project directory.

#### Compile #### 
1. Open "RootsToolC.sln" in Visual Studio. If it is unable to open the solution and gives an error like "Could not open 'xxx.VC.VC.opendb': The process cannot access the file because it is being used by another process.", in Team Explorer, click Home|PlantRoot and go to "Settings -- Repository Settings". Under "Ignore & Attributes Files", generate a new .gitignore file if there isn't one and click edit. Add "db.lock" to the list under "#Build Results". Save it and re-open Visual Studio.
2. In Class View, Right click on "RootsTool" and select "Set as Startup Project"
3. Select Build->Build Solution
   - If it reports error message like "LINK : fatal error LNK1104: cannot open file ...", check if your dependencies are correctly installed and CMakeLists.txt is correctly modified. If it still doesn't work, go to "Solution Explorer -- right click RootsTool -- Properities" and manually add these directories:
     - VC++ Directories -> Include Directories
     - VC++ Directories -> Library Directories
     - C/C++ -> General -> Additional Include directories
     - Linker -> General -> Additional Library Directories
     - Linker -> Input -> Additional dependency
   - If there is error message relate to "conflict target machine", check Properities -> Configuration Properties -> Linker -> Advanced -> Target Machine. Make sure it's "MachineX64"
4. A successful build should give only warnings and have "3 succeeded, 0 failed, 1 skipped". Check your output folder (parallel to your project folders and named "python"), there should be files that looks like "RootsEXE.xxx" and "RootsTool.xxx".

#### Finally ####
Congrats! Now you are finished on the C++ part. Go to https://github.com/chunyuan1/PlantRootsRelease for further instruction on the python part.