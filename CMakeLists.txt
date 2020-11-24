

cmake_minimum_required(VERSION 2.8)
project (ConstrainedDelaunay)

include_directories("./source")



set(CTL_HEADERS
	source/BlockList.h
	source/CGrid.h
	source/ConstraintMap.h
	source/ctl.h
	source/DelaunayTriangulation.h
	source/Edge.h
	source/ID.h
	source/LocationResult.h
	source/QTriangulate.h
	source/QuadEdge.h
	source/Subdivision.h
	source/TIN.h
	source/Util.h
	source/Vector.h
	source/Vertex.h
	)

set(CTL_SOURCES 
	source/CGrid.cpp
	source/ConstraintMap.cpp
	source/DelaunayTriangulation.cpp
	source/Edge.cpp
	source/ID.cpp
	source/LocationResult.cpp
	source/QTriangulate.cpp
	source/QuadEdge.cpp
	source/Subdivision.cpp
	source/TIN.cpp
	source/Util.cpp
	source/Vector.cpp
	source/Vertex.cpp
	)

add_executable(simple_example ${CTL_SOURCES} ${CTL_HEADERS} "examples/simple/delaunay_simple.cpp")

add_library(ctltriangulator SHARED ${CTL_SOURCES} ${CTL_HEADERS} "source/DelaunayTriangulationWrapper.cpp")
