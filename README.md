# PlantRoot

## Compile Instruction ##

#### Dependency ####
	1. Compile Boost_1_65_1
	2. Download glfw-3.2.1 64bit version
	3. Download glad
	4. Install python3.5 64bit version

#### Modify CMakeLists.txt #### 
	1. Substitute include path to corresponding path in your computer
	2. Substitute library path to corresponding path in your computer

#### Compile #### 
	1. Navigate to project folder, enter command:
```
mkdir build
cd build
cmake -G "Visual Studio 15 2017 Win64" ..
```

	2. Open sln file in vs2017

	Right click 'RootsTool' select “Set as Startup Project”

Select Build->Build Solution
- If it reports error message like "LINK : fatal error LNK1104: cannot open file ...", check dependency (not necessary if directories are correcyly defined in CMakeLists.txt) in Solution Explorer -> right click RootsTool -> properities
  - VC++ Directories -> Include Directories
  - VC++ Directories -> Library Directories
  - C/C++ -> General -> Additional Include directories
  - Linker -> General -> Additional Library Directories
  - Linker -> Input -> Additional dependency

