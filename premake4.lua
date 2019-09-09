solution "voronoi"
	configurations { "release", "debug" }
	-- switch gcc / clang
    -- premake.gcc.cc = "gcc-8"
    -- premake.gcc.cxx = "g++-8"
    premake.gcc.cc = "clang"
    premake.gcc.cxx = "clang++"

	includedirs { ".", "src/gKit" }

	configuration "release"
		targetdir "bin/release"
		defines { "GK_RELEASE" }
		flags { "OptimizeSpeed" }
		buildoptions{ "-O3" }
		linkoptions{ "-O3" }
	
	configuration "debug"
		targetdir "bin/debug"
		defines { "DEBUG" }
		flags { "Symbols" }

	configuration "linux"
		buildoptions { "-mtune=native -march=native -std=c++11 -W -Wall -Wextra -Wsign-compare -Wno-unused-parameter -Wno-unused-function -Wno-unused-variable", "-pipe" }
		buildoptions { "-flto"}
		linkoptions { "-flto"}
		buildoptions { "-fopenmp" }
		linkoptions { "-fopenmp" }
		links { "SDL2", "SDL2_image" }

	configuration { "linux", "debug" }
		buildoptions { "-g"}
		linkoptions { "-g"}
    
project("voronoi")
	language "C++"
	kind "ConsoleApp"
	targetdir "bin"
	files { "src/gKit/*", "src/main.cpp" }

project("test-kdtree")
	language "C++"
	kind "ConsoleApp"
	targetdir "bin"
	files { "src/gKit/*", "src/test/kdtree.cpp" }
