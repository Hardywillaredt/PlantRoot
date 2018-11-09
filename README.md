# PlantRoot

### Compile Instruction ###
Environment: Windows, Python 3.5 64bit, Visual Studio 2017
#### Dependency ####
1. Compile Boost_1_65_1
   
   "libboost_python" library is needed. 
   ```
   .\b2 --address-mode=64 --with-python
   ```
2. Download glfw-3.2.1 64bit version
3. Download glad

It's not necessary to have exactly the same version as listed above. But we didn't test on other versions.
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
2. Open .sln file by double click
3. Right click 'RootsTool' select “Set as Startup Project”
4. Select Build->Build Solution
   - If it reports error message like "LINK : fatal error LNK1104: cannot open file ...", check dependency (not necessary if directories are correcyly defined in CMakeLists.txt) in Solution Explorer -> right click RootsTool -> Properities. You probabily need to manually add these directories
     - VC++ Directories -> Include Directories
     - VC++ Directories -> Library Directories
     - C/C++ -> General -> Additional Include directories
     - Linker -> General -> Additional Library Directories
     - Linker -> Input -> Additional dependency
   - If there is error message relate to "conflict target machine", check Properities -> Configuration Properties -> Linker -> Advanced -> Target Machine. Make sure it's "MachineX64"
