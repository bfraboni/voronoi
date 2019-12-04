solution "Voronoi"
	configurations { "release", "debug" }
	platforms { "x64" }
	includedirs { "src/", "src/gKit" }
	
	configuration "debug"
		targetdir "bin/debug"
		defines { "DEBUG" }
		if _PREMAKE_VERSION >="5.0" then
			symbols "on"
		else
			flags { "Symbols" }
		end

	configuration "release"
		targetdir "bin/release"
		if _PREMAKE_VERSION >="5.0" then
			optimize "speed"
		else
			flags { "OptimizeSpeed" }
		end

	configuration "linux"
		buildoptions { "-mtune=native -march=native" }
		buildoptions { "-std=c++11" }
		buildoptions { "-W -Wall -Wextra -Wsign-compare -Wno-unused-parameter -Wno-unused-function -Wno-unused-variable", "-pipe" }
		buildoptions { "-flto"}
		linkoptions { "-flto"}
		buildoptions { "-fopenmp" }
		linkoptions { "-fopenmp" }
		links { "SDL2", "SDL2_image" }

	configuration { "linux", "debug" }
		buildoptions { "-g"}
		linkoptions { "-g" } -- bugfix : premake4 ne genere pas le flag debug pour le linker... 

	configuration { "windows" }
		defines { "WIN32", "NVWIDGETS_EXPORTS", "_USE_MATH_DEFINES", "_CRT_SECURE_NO_WARNINGS" }
		defines { "NOMINMAX" } -- allow std::min() and std::max()

		-- extern libs
		includedirs { "extern/visual/include" }
		libdirs { "extern/visual/lib" }
		links { "SDL2", "SDL2main", "SDL2_image" }
		
		buildoptions { "/openmp" } 
		
		if _PREMAKE_VERSION >="5.0" then
			system "Windows"
			architecture "x64"
			flags { "MultiProcessorCompile", "NoMinimalRebuild" }
			disablewarnings { "4244", "4305" }
		end
		
	configuration { "windows", "release" }
		if _PREMAKE_VERSION >="5.0" then
			flags { "LinkTimeOptimization" }
		end
		
	configuration "macosx"
		frameworks= "-F /Library/Frameworks/"
		buildoptions { "-std=c++11" }
		defines { "GK_MACOS" }
		buildoptions { frameworks }
		linkoptions { frameworks .. " -framework SDL2 -framework SDL2_image " }

-- common files
gkit_files = { "src/gKit/*.cpp", "src/gKit/*.h" }

project("voronoization")
	language "C++"
	kind "ConsoleApp"
	targetdir "bin"
	includedirs { "src/", "src/extern/" }
	files { gkit_files , "src/main.cpp" }
