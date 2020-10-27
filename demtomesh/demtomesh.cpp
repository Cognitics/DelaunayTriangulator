// demtomesh.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "cdb_tile/Tile.h"
#include "ccl/FileInfo.h"
#include "ctl/Vector.h"
#include "ccl/gdal.h"
#include "elev/SimpleDEMReader.h"
#include "ctl/DelaunayTriangulation.h"
#include "ctl/TIN.h"
#include <scenegraph/Scene.h>
#include <scenegraphflt/scenegraphflt.h>
#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>
#undef min



using namespace cognitics::cdb;

/*! \brief Saves triangulation results to a file

	First the caller will add all the vertices using AddVertex(). Each order gets an index

	Derived classes implement the AddPolygon(), AddVertex() and AddTriangle() methods
	AddPolygon() is optional
*/

class MeshWriter
{
public:
	virtual ~MeshWriter() {}
	virtual bool AddVertex(const ctl::Point& p) = 0;
	virtual bool AddFacet(int index1, int index2, int index3) = 0;
	virtual bool AddPolygon(const ctl::PointList& polygon) { return true; }
	virtual bool Close() { return true; }
};

//! \brief Writer that generates a simple .obj mesh file
class ObjMeshWriter : public MeshWriter
{
	std::string _fileName;
	std::ofstream _fd;
	int _vertexCount = 0;

public:
	ObjMeshWriter(const std::string& fileName) : _fileName(fileName)
	{}
	~ObjMeshWriter()
	{
		Discard();
	}

//! \brief Closes the file and deletes it so that it is not left around
	void Discard()
	{
		Close();
	}
	bool Open()
	{
		if (!_fd.is_open())
		{
			_fd.open(_fileName, std::ios_base::binary | std::ios_base::out | std::ios_base::trunc);
			_fd << std::setprecision(3);
		}
		return _fd.is_open();
	}
	bool Close()
	{
		if (_fd.is_open())
			_fd.close();
		_vertexCount = 0;
		return true;
	}
	bool AddVertex(const ctl::Point& p)
	{
		if (!Open())
			return false;

		// We save the geometric coordinates only
		char str[200];
		double r = 0, g = 0, b = 0;
		if (p.z <= 0)
			b = 1.0;
		else
			r = p.z/300;
		sprintf_s(str, "v %.3f %.3f %.3f %.3f %.3f %.3f\n", p.x, p.y, p.z, r, g, b);
	
//        _fd << "v " << p.x << " " << p.y << " " << p.z "\n";
		_fd << str;
		++_vertexCount;
		return true;
	}

	virtual bool AddFacet(int index1, int index2, int index3)
	{
		if (!Open())
			return false;
		if (index1 < 0 || index1 >= _vertexCount ||
			index2 < 0 || index2 >= _vertexCount ||
			index3 < 0 || index3 >= _vertexCount)
				return false;
		// Note that OBJ indexes start at 1
		_fd << "f " << index1+1 << " " << index2+1 << " " << index3+1 << "\n";
		return true;
	}

}; // ObjMeshWriter


 //! \brief The quantized elevations are given in cm, and the low 8 bits 
 //         contain the flags POINT_xxx
int MeshElevationToQuantized(double elevation)
{
	return ((int)(elevation * 100)) << 8;
}

int MeshQuantizedZ(int quantized)
{
	return quantized >> 8;
}

double MeshQuantizedToElevation(int quantized)
{
	return (double)MeshQuantizedZ(quantized)/100.0;
}


//! \brief Quantized elevation flags
//      To facilitate the traversal of the elevation array this flags indicate 
//      how each point should be handled
enum {
	MESHPOINT_BOUNDARY = 1,  // The point is the bounday of the target area
	MESHPOINT_REMOVED = 2,   // The point was discarded, considered not relevant
	MESHPOINT_EDGE = 4,      // The point was identified as an edge between removed areas, so it should be preserved
	MESHPOINT_PENDING = 8    // The point has been already pushed in the analysis stack
};



/*! \brief Generic interface to tweak the selected points in a line of quantized elevations
	The line is defined by its length and the increment between values, so it can
	handle lines which are not contiguous in the array

Derived classes implement the Process() method that carry out the processing. The intention is
that the should only remove some dots from the triangulation, or preserve them.

\code
		class Simplifier : public DEMesher::LineSelector 
		{
		void Process();
		};
		Simplifier s;

		// Lines are 65 values long, the distance between them is 513 values
		// Force to preserve the first two values and the last two values
		int *firstPoint = GetFirstPoint();
		firstPoint[0*513] |= MESHPOINT_EDGE;
		firstPoint[1*513] |= MESHPOINT_EDGE;
		firstPoint[63*513] |= MESHPOINT_EDGE;
		firstPoint[64*513] |= MESHPOINT_EDGE;
		s.Select(firstPoint, 65, 513);
\endcode
*/
class LineSelector
{
protected:
	int *_first = nullptr; // Pointer to the first quantized point
	int _n = 0;    // Number of points
	int _step = 1; // Step between points

	virtual ~LineSelector() {}

/*! \brief Implements simplification of the line

	The fields _first, _n and _step contain the location of tha target line
*/
	virtual void Process() = 0;

	// A few simple macros to manipulate the quantized elevations
	// Reference to an entry
	int& Item(int x) { return _first[x*_step]; }

	// Return the elevation of a line point, without the flags
	int Z(int x) { return MeshQuantizedZ(Item(x)); }

	void Remove(int x)    { Item(x) |= MESHPOINT_REMOVED; };   // Remove from triangulation
	void Preserve(int x) { Item(x) &= (~MESHPOINT_REMOVED); }; // Marks the value to be used
	void KeepEdge(int x) { Item(x) |= MESHPOINT_EDGE; Item(x) &= (~MESHPOINT_REMOVED); }; // Edge to be forced into the triangulation
	bool IsRemoved(int x) { return (Item(x) & MESHPOINT_REMOVED) != 0; }
	bool IsEdge(int x)    { return (Item(x) & MESHPOINT_EDGE) != 0; }

public:
//! \brief Sends a line (defined by the start and the step interval) to a line selector
	void Select(int* line, int n, int step)
	{
		if (!line || n <= 0)
			return;
		_first = line;
		_n = n;
		_step = step;
		Process();
	}
};

/*! \brief Generic interface that selects which vertices are to be preserved in the mesh

    Derived classes implement the Process() method. When called the quantized elevation grid
	is in the _data member
*/
class MeshSelector
{
protected:
	// Points to the quantized elevation grid
	int* _data = nullptr;
	int _qSize = 0;

	int* ptr(int y)
	{
		if (y <= 0)
			return _data;
		if (y >= _qSize) 
			y = _qSize - 1;
		return _data + y * _qSize;
	}
	
	virtual ~MeshSelector() {}
	virtual void Process() = 0;
public:
	void Select(int* grid, int gridSize)
	{
		if (grid && gridSize >= 2)
		{
			_data = grid;
			_qSize = gridSize;
			Process();
		}
	}
}; // MeshSelector


/*!
	\brief Generates the mesh for a CDB tile.

	This implementation assumes that the DEM files are aligned with the tile.
*/
class DEMesher
{
	// Location of the geoTile for this tile
	Tile _tile;

	// The base for the CDB
	std::string _root;

	// Quantized altitude vector. The altitudes are defined at the corners of the
	// grid, starting at the NorthWest corner.
	// Note that the first entry of the table is the altitude of the NW corner
	// of the target area. It is not shifted as in the CDB specification. So the
	// array contains one more row and one more column than the tile size
	int _qSize = 0;
	std::vector<int> _quantized;

	int* ptr(int y) 
	{ 
		assert(y >= 0 && y < _qSize);
		if (y < 0) y = 0;
		if (y >= _qSize) y = _qSize - 1;
		return _quantized.data() + _qSize * y; 
	}

	public: 

 //! \brief Creates a mesher for the tile at a certain latitude and longitude
	DEMesher(const std::string &cdbRoot, const Tile& t) :
		  _root(cdbRoot),
		  _tile(t) {}

	DEMesher(const std::string &cdbRoot, Latitude south, Longitude west, int LOD) :
		_root(cdbRoot)
	{
		TileLatitude geoLat(south);
		int tileWidth = get_tile_width(geoLat);
		TileLongitude geoLng(geoLat, west);
		if (LOD < 0)
			LOD = 0;
		int scale = (1 << LOD);
		int up = static_cast<int>((south.value() - geoLat.value()) * scale);
		int right = static_cast<int>((west.value() - geoLng.value()) * scale)/tileWidth;

		// Use DS001, S001, T001, tiles for the terrain elevation
		_tile.setDataset(Dataset(Dataset::Elevation));
		_tile.setLod(LOD);
		_tile.setCs1(1);
		_tile.setCs2(1);
		_tile.setUref(up);
		_tile.setRref(right);

		// The coordinates of the tile
		CoordinatesRange range(geoLng.value() + right * tileWidth / (double)scale,
			geoLng.value() + (right + 1) * tileWidth / (double)scale,
			geoLat.value() + up / (double)scale,
			geoLat.value() + (up + 1)/(double) scale);
		_tile.setCoordinates(range);
	}

	Tile& TargetTile() { return _tile; }

	bool CheckForCDB()
	{
		_root = ccl::standardizeSlashes(_root);
		if(!_root.size() || _root.back()!='/')
			_root += "/";
		ccl::FileInfo fi;
		return ccl::directoryExists(_root);
	}

	void Generate(const Tile &north, const Tile &east, const Tile &northEast);

	void Triangulate(scenegraph::Scene *saved);

private:
	void SelectBoundary(LineSelector& processor);
	bool CopyGrid(elev::SimpleDEMReader &reader, int x, int y, int subsample=2);
	void SaveQuantized(const char* fileName);
};

//! \brief Writes a JSON with the quantized values, for analysis outside
void DEMesher::SaveQuantized(const char* fileName)
{
	std::ofstream fd(fileName, std::ios_base::out | std::ios_base::binary);
	fd << "{ \"data\": [\n";
	for(int y = 0; y < _qSize; y++)
	{
		int* src = ptr(y);
		if(y > 0)
			fd << ",\n";
		fd << "[ ";
		for(int x = 0; x < _qSize; x++)
		{
			if(x > 0)
			{
				fd << ",";
				if((x % 100) == 0)
					fd << "\n  ";
			}
			fd << src[x];
		}
		fd << " ]";
	}
	fd << " \n] }\n";
}


//! \brief Copies the grid from a reader to the quantized grid
//! \param x,y are the indexes of the top-left corner of the incoming data
//       relative to the tile grid. They can be negative
bool DEMesher::CopyGrid(elev::SimpleDEMReader &reader, int x, int y, int subsample)
{
	std::vector<double> points;
	int width = reader.getWidth();
	int height = reader.getHeight();
	bool valid = reader.getGrid(points);
	if(!valid)
		return false;
	int srcx = 0, srcy = 0;
	if(x < 0)
	{
		srcx = -x;
		x = 0;
	}
	if(y < 0)
	{
		srcy = -y;
		y = 0;
	}

	// Copy the sections that overlap in the two arrays to the destination. Note that
	// we are subsampling the input.
	int copyWidth = std::min(_qSize - x, width / subsample - srcx);
	for(;y<_qSize && srcy<height/subsample; y++, srcy++)
	{ 
		int* dst = ptr(y) + x;
		double* srcRow = points.data() + width * subsample *srcy + subsample * srcx;
		for(int i = 0; i < copyWidth; i++)
			dst[i] = MeshElevationToQuantized(srcRow[i * subsample]);
	}
	return true;
}

/*! \brief Removes triangulation points from flat areas

	The points at the boundary must be marked with MESHPOINT_BOUNDARY so that the algorithm
	does not exceed the limits of the tile.

	The function checks areas of the image where the point altitudes are within a range.
	All the points inside are removed.

	Elevation zero is handled differently because it corresponds to sea level. Only 
	posts with elevation zero are combined.

	The algorithm is connected components:
	   -for each point, put in stack each neighbor which is within range AND is not
		an edge, boundary, or discarded
	   -if all four neighbors are within range, mark the point as discarded. Otherwise
		mark as edge
	   -pop next point from stack
*/
class SimplifyFlat : public MeshSelector
{
public:
	SimplifyFlat(double toleranceMeters) :
		_tolerance(toleranceMeters) {}

private:
	double _tolerance = 0;
	std::vector<int> _stack; // The stack contains the positions of the next point pending to check. It contains indexes

	void Process()
	{
		int tolerance = 0; // Quantized tolerance. It is zero when selecting points at sea level

		int limit = _qSize * (_qSize-1);

		// The altitude range of the current region. We need to make sure the range does 
		// not exceed the tolerance
		int low, high;

		auto CheckNeighbor = [&](int index)->bool
		{
			assert(index >= 0 && index < _qSize* _qSize);
			int altitude = MeshQuantizedZ(_data[index]);
			if (altitude == 0 && tolerance != 0)
				return false;
			bool inRange = altitude >= high - tolerance && altitude <= low + tolerance;
			if (inRange)
			{
				if (low > altitude)
					low = altitude;
				else if (high < altitude)
					high = altitude;
				if ((_data[index] & (MESHPOINT_BOUNDARY | MESHPOINT_REMOVED | MESHPOINT_EDGE | MESHPOINT_PENDING)) == 0)
				{
					_data[index] |= MESHPOINT_PENDING;
					_stack.push_back(index);
				}
			}
			else if (altitude != 0 && tolerance != 0 && _data[index] & (MESHPOINT_BOUNDARY | MESHPOINT_EDGE))
			{
				// If the neighbor is an edge and the altitude is close to the limit, let's
				// discard it anyway
				int highTolerance = (tolerance * 12) / 10;
				inRange = altitude >= high - highTolerance && altitude <= low + highTolerance;
			}

			return inRange;
		};

		int removed = 0;
		int edges = 0;

		// Check all the point, skipping the first and last rows because they are boundaries
		for (int testIndex = _qSize; testIndex < limit; testIndex++)
		{
			int v = _data[testIndex];
			if (v & (MESHPOINT_BOUNDARY | MESHPOINT_REMOVED | MESHPOINT_EDGE))
				continue;
			low = high = MeshQuantizedZ(v);

			// For a point at sea level the tolerance is zero so that the area only contains points at sea level
			tolerance = (low == 0) ? 0 : MeshQuantizedZ(MeshElevationToQuantized(_tolerance));
			_stack.push_back(testIndex);
			while (_stack.size())
			{
				int idx = _stack.back();
				_stack.pop_back();
				_data[idx] &= ~MESHPOINT_PENDING; // Clear pending status
				if (_data[idx] & (MESHPOINT_BOUNDARY | MESHPOINT_REMOVED | MESHPOINT_EDGE))
					continue;
				int inRange = 0;
				if (CheckNeighbor(idx + 1))  // Right neighbor
					inRange++;
				if (CheckNeighbor(idx - 1))  // Left neighbor
					inRange++;
				if (CheckNeighbor(idx + _qSize)) // Below neighbor
					inRange++;
				if (CheckNeighbor(idx - _qSize)) // Above neighbor
					inRange++;
				if (inRange == 4)
				{
					_data[idx] |= MESHPOINT_REMOVED;
					removed++;
				}
				else
				{
					_data[idx] |= MESHPOINT_EDGE;
					edges++;
				}
			}
		}
		std::cout << removed << " removed points, " << edges << " edge points" << std::endl;
	}
}; // SimplifyFlat

/*! \brief Removes points which are approximately colinear with their neighbors
*/
class SimplifySlope : public MeshSelector
{
public:
	SimplifySlope(double toleranceRatio) :
		_tolerance(toleranceRatio) {}

private:
	// Margin of error, as a fraction of the range between maximum and minimum
	double _tolerance = 0;
	void Process()
	{
		double tolerance = _tolerance;
		auto IsInline = [tolerance](int e1, int e2, int e3)
		{
			if ((e1 | e2 | e3) & MESHPOINT_REMOVED)
				return false;
			// Preserve shoreline edges. If a point is at sea level, it can only
			// be removed is the neighbors are at sea level
			e1 = MeshQuantizedZ(e1);
			e2 = MeshQuantizedZ(e2);
			e3 = MeshQuantizedZ(e3);
			if (e2 == 0)
				return e1 == 0 && e3 == 0;
			int error = e2 - (e1 + e3) / 2;
			int margin = (int)(tolerance * (e3 - e1));
			if (margin < 0)
				margin = -margin;
			if (error < 0)
				error = -error;
			return error <= 100+margin;
		};
		int removed = 0;
		int accepted = 0;
		for (int y = 1; y < _qSize-1; y++)
		{
			int* row = ptr(y);
			int* above = ptr(y - 1);
			int* below = ptr(y + 1);
			for (int x = 1; x < _qSize - 1; x++)
			{
				// There are 8 lines that combine the center pixel with the neighbor
				// Check if the center is aligned with any of them. There must be a better way
				if (IsInline(row[x - 1], row[x], row[x + 1]) ||
					IsInline(above[x], row[x], below[x]) ||
					IsInline(above[x - 1], row[x], below[x + 1]) ||
					IsInline(above[x + 1], row[x], below[x - 1]) ||

					IsInline(above[x - 1], row[x], row[x + 1]) ||
					IsInline(below[x - 1], row[x], row[x + 1]) ||
					IsInline(row[x - 1], row[x], above[x + 1]) ||
					IsInline(row[x - 1], row[x], below[x + 1])
					)
				{
					removed++;
					row[x] |= MESHPOINT_REMOVED;
				}
				else if ((row[x] & MESHPOINT_REMOVED) == 0)
					accepted++;
			}
		}
		std::cout << removed << " colinear points, " << accepted << " remaining points" << std::endl;
	}
}; // SimplifySlope


/*! \brief Removes points from a line keeping the error below a limit
 
	The function checks how much error would be generated if the point was removed
	and the previous and next segments were joined. If the error is less than the
	objective the point is removed and the two segments joined.

	Point N considers two segments [M,N] and [N,N+1]. The algorithm evaluates the
	error in the range [M, N+1]

	Points which are marked as MESHPOING_EDGE are always preserved, points 
	marked as MESHPOINT_REMOVE are always removed.
*/
class LineSimplifyByError : public LineSelector
{
	// Error allowed in the approximation
	double _toleranceMeters = 0;

	// Maximum spacing between preserved points. We want to keep a few points in the
	// line to make triangulation simpler
	static int const MAX_DISTANCE = 32;

public:
	LineSimplifyByError(double toleranceMeters) :
		_toleranceMeters(toleranceMeters)
	{}

	void Process()
	{
		int xprev = 0, xnext = 0; // index of the previous point preserved
		int x;
		int tolerance = MeshQuantizedZ(MeshElevationToQuantized(_toleranceMeters));

		// Calculates the error at point X, given the current segment
		auto ZError = [&](int x)->int
		{
			int zprev = Z(xprev);
			int znext = Z(xnext);
			// Note that the multiplication will NOT overflow because the range [xprev,xmax] is less than 32
			int z = Z(x);
			int error = z - (zprev + ((x - xprev) * (znext - zprev)) / (xnext - xprev));
			return error >= 0 ? error : -error;
		};

		Preserve(0);
		Preserve(_n-1);

		for (x = 1; x < _n-1; x++)
		{
			if (IsRemoved(x))
				continue;
			if (x >= xprev + MAX_DISTANCE || IsEdge(x))
			{
				Preserve(x);
				xprev = x;
				continue;
			}
			xnext = x + 1;
			// Calculate the error for all the points in the range
			int error = 0;
			for (int i = x; i > xprev && error <= tolerance; i--)
				error = ZError(x);
			if (error <= tolerance)
				Remove(x);
			else
			{
				Preserve(x);
				xprev = x;
			}
		}
	}
}; // LineSimplifyByError


/*! \brief Preserve sea level boundaries

	 This class inserts grid points where there is a transition from sea to land
*/
class KeepSeaLevelInLine : public LineSelector
{
	// Maximum and minimum sea level
	double _maxLevel, _minLevel;
	void Process()
	{
		int low = MeshQuantizedZ(MeshElevationToQuantized(_minLevel));
		int high = MeshQuantizedZ(MeshElevationToQuantized(_maxLevel));
		auto InRange = [&](int x)->bool
		{
			int z = Z(x);
			return (z >= low && z <= high);
		};

		for (int x = 1; x < _n-1; x++)
		{
			if (InRange(x))
			{
				if (!InRange(x+1))
				{
					Preserve(x);
					Preserve(x+1);

				}
				if (!InRange(x-1))
				{
					Preserve(x);
					Preserve(x-1);
				}
			}
		}
	}
public: 
	KeepSeaLevelInLine(double maxLevelMeters, double minLevelMeters) :
		_maxLevel(maxLevelMeters),
		_minLevel(minLevelMeters) {}
};

//! \brief Marks the points in a line as boundary, the corners as edges
class BoundaryMarker : public LineSelector
{
	void Process()
	{
		for (int x = 0; x < _n; x++)
			Item(x) |= MESHPOINT_BOUNDARY;
		Item(0) |= MESHPOINT_EDGE;
		Item(_n-1) |= MESHPOINT_EDGE;
	}
};


class TestLineSimple
{
public:
	void DoTest(std::vector<double> values, std::vector<int> constraints, int step)
	{
		std::vector<int> elev(values.size()*step, -10000);
		for (int i = 0; i < values.size(); i++)
			elev[i*step] = MeshElevationToQuantized(values[i]) | MESHPOINT_BOUNDARY;
		for (int c : constraints)
			if (c >= 0 && c < values.size())
				elev[c*step] |= MESHPOINT_EDGE;
		LineSimplifyByError obj1(10);
		obj1.Select(elev.data(), (int)elev.size(), step);
	}

	TestLineSimple()
	{
		std::vector<int> empty;
		std::vector<int> keep1 = { 1 };
		std::vector<int> keep2 = { 1, 3 };
		std::vector<double>v0;
		std::vector<double>v1 = { 1, 2, 3, 4, 5 };
		std::vector<double>v2 = { 1, 2, 20, 4, 5 };
		std::vector<double>v3 = { 1, 2, -20, 4, 5 };
		DoTest(v0, empty, 1);
		DoTest(v0, keep1, 1);
		DoTest(v1, empty, 1);
		DoTest(v1, keep1, 1);
		DoTest(v1, keep2, 1);
		DoTest(v2, empty, 1);
		DoTest(v2, keep1, 1);
		DoTest(v2, keep2, 1);
		DoTest(v3, empty, 1);
		DoTest(v3, keep1, 1);
		DoTest(v3, keep2, 1);
	}
};

//TestLineSimple t;

void DEMesher::SelectBoundary(LineSelector& processor)
{
	// Run the four corners in the standard order (minx, miny) (maxx, miny) (maxx,maxy) (minx,maxy)
	processor.Select(ptr(0),                 _qSize, 1);
	processor.Select(ptr(0)+_qSize-1,        _qSize, _qSize);
	processor.Select(ptr(_qSize-1)+_qSize-1, _qSize, -1);
	processor.Select(ptr(_qSize - 1),        _qSize, -_qSize);
}


void DEMesher::Generate(const Tile &north, const Tile &east, const Tile &northEast)

{
	if (!CheckForCDB())
		return;

	auto fileName = _root+_tile.getFilename(".tif");
	OGRSpatialReference oSRS; // GDAL reference system
	oSRS.SetWellKnownGeogCS("WGS84");
	int width, height;
	int subsample = 2;

	{
		std::cout << "Center: " << fileName << std::endl;
		elev::SimpleDEMReader reader(fileName, oSRS);
		if (!reader.Open())
			return;
		width = reader.getWidth()/subsample;
		height = reader.getHeight()/subsample;

		// Generate the grid at half the source resolution, one additional as guard
		_qSize = std::min(width, height) + 1;
		_quantized.resize(_qSize * (_qSize+1), MeshElevationToQuantized(-10000));
		// Center reader skips the NORTH edge
		CopyGrid(reader,0,1,subsample);
		// Copy the missing top and left edges
		for (int i = 0; i < _qSize; i++)
		{
			ptr(0)[i] = ptr(1)[i];
			ptr(i)[_qSize - 1] = ptr(i)[_qSize - 2];
		}
		ptr(0)[_qSize - 1] = ptr(0)[_qSize - 2];
	}
	{
		fileName = _root + north.getFilename(".tif");
		std::cout << "North: " << fileName << std::endl;
		elev::SimpleDEMReader reader(fileName, oSRS);
		if (reader.Open())
			// Copy the top-most row
			CopyGrid(reader,0,-height+1,subsample);
	}
	{
		fileName = _root + east.getFilename(".tif");
		std::cout << "East: " << fileName << std::endl;
		elev::SimpleDEMReader reader(fileName, oSRS);
		if (reader.Open())
			// Copy the right-most column
			CopyGrid(reader,width,1,subsample);
	}
	{
		fileName = _root + northEast.getFilename(".tif");
		std::cout << "NorthEast: " << fileName << std::endl;
		elev::SimpleDEMReader reader(fileName, oSRS);
		if (reader.Open())
			// Copy the top-right corner
			CopyGrid(reader,width,-height+1,subsample);
	}
	
	BoundaryMarker marker;
	SelectBoundary(marker);

	double tolerance = 20;
	double seaLevel = 0;
	SimplifyFlat mesher(tolerance);
	mesher.Select(_quantized.data(), _qSize);
	SimplifySlope flattener(0.2);
	flattener.Select(_quantized.data(), _qSize);

	// Run the mesh simplification and the sea level detail on the boundary
	LineSimplifyByError simplify(tolerance / 2);
	KeepSeaLevelInLine forceSeaLevel(seaLevel, seaLevel);
	SelectBoundary(simplify);
	SelectBoundary(forceSeaLevel);


	SaveQuantized("c:\\build\\temp\\quantized.json");
}

/*! \brief Saves the triangulated mesh into the scene

	The points in the DEM elevation grid that are preserved make a triangle mesh. The boundary
	of the grid is a constraint so that all the edges are included.

	The mesh becomes a Face in the scene. The (X,Y) coordinates are Y pointing up, the expectation 
	is that the scene matrix will set the right orientation.

	REVISIT:: Need to figure out the texture and the UV coordinates.
*/
void DEMesher::Triangulate(scenegraph::Scene* scene)
{
	// The boundary rectangle is counterclockwise, but on a 'normal' (Y up) axis
	// the best way to build it is: (xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax)
	int* top = ptr(0);
	int* bottom = ptr(_qSize - 1);
	ctl::PointList boundary{ 
		ctl::Point(       -1,        -1, MeshQuantizedToElevation(top[0])),
		ctl::Point(_qSize,        -1, MeshQuantizedToElevation(top[_qSize-1])), 
		ctl::Point(_qSize, _qSize, MeshQuantizedToElevation(bottom[_qSize-1])), 
		ctl::Point(-1,        _qSize, MeshQuantizedToElevation(bottom[0])) 
	};

	ctl::DelaunayTriangulation dt(boundary);

	// All the points that have not been removed are constraints

	// First, insert the constrainted outside contour
	std::vector<ctl::Point> contour;

	// The coordinate system for the display has the origin on the bottom-right, but the images have the
	// origin at the top-left.
	auto ScenePoint = [&](int x, int y, int altitude)->ctl::Point
	{
		return ctl::Point(x, _qSize-y-1, MeshQuantizedToElevation(altitude));
	};

	auto AddToContour = [&](int x, int y)
	{
		int altitude = ptr(y)[x];
		if ((altitude & MESHPOINT_REMOVED) == 0)
			contour.push_back(ScenePoint(x,y,altitude));
	};
	int x=0, y=0;
	for (; x < _qSize - 1; x++)
		AddToContour(x, y);
	for (; y < _qSize-1; y++)
		AddToContour(x, y);
	for (; x>0; x--)
		AddToContour(x, y);
	// Note that the first point is repeated at the end
	for (; y>=0; y--)
		AddToContour(x, y);

	int failedPoints = 0;
	if (!dt.InsertConstrainedLineString(contour))
		failedPoints = (int)contour.size();
	int insertedPoints = (int)contour.size();
	for (y=1; y<_qSize-1; y++)
		for (x = 1; x < _qSize-1; x++)
		{
			int altitude = ptr(y)[x];
			if ((altitude & MESHPOINT_REMOVED) == 0)
			{
				insertedPoints++;
				if (dt.InsertConstrainedPoint(ScenePoint(x,y,altitude)) == 0)
					failedPoints++;
			}
		}
	std::vector<ctl::Edge*> edges;
	edges = dt.GatherTriangles(ctl::PointList());
	ctl::TIN container(&dt, edges);
	if (scene)
	{
		scene->hasVertexNormals = true;
		scene->faces.emplace_back();
		scenegraph::Face& face = scene->faces.back();
		face.primaryColor = scenegraph::Color(0.5f, 1.0f, 0.5f);
		face.alternateColor = scenegraph::Color(0.5f, 1.0f, 0.5f);

		bool success = true;
		auto &verts = container.verts;
		auto &normals = container.normals;

		// Indexes of the corners. The last entry always contains the last enumerated corner.
		int corners[4] = { -1,-1,-1,-1 };
		int numPoints = 0;
		int numCorners = 0;
		for (int i=0; i<verts.size(); i++)
		{
			ctl::Point v = verts[i];
			if (v.x < 0 || v.x >= _qSize)
			{
				if (numCorners<4)
					corners[numCorners] = corners[3] = numPoints;
				numCorners++;
			}
			face.addVert(sfa::Point(v.x, v.y, v.z));
			face.vertexNormals.push_back(sfa::Point(normals[i].x, normals[i].y, normals[i].x));
			numPoints++;
		}

		int numTriangles = 0;
		auto triangles = container.triangles;
		auto IsCorner = [&](int index)->bool { return index<=corners[3] && (index==corners[3] || index == corners[0] || index == corners[1] || index == corners[2]); };
		for (size_t i = 0; i+2 < triangles.size(); i += 3)
		{
			if (!IsCorner(triangles[i]) && !IsCorner(triangles[i+1]) && !IsCorner(triangles[i+2]))
			{
				success = success && face.AddFacet(triangles[i+0], triangles[i+1], triangles[i+2]);
				++numTriangles;
			}
		}

		
		// Add a new face with blue color that marks the sea level
		if (numCorners >= 4)
		{
			scene->faces.emplace_back();
			scenegraph::Face& seaSurface = scene->faces.back();
			seaSurface.primaryColor = scenegraph::Color(0.3, 0.3f, 1.0f);
			seaSurface.alternateColor = scenegraph::Color(0.3f, 0.3f, 1.0f);
			for (int c = 0; c < 4; c++)
			{
				ctl::Point v = verts[corners[c]];
				// Corner altitude just above zero so that the surface hides the mesh below
				v.z = 0.001;
				if (v.x < 0) v.x = 0;
				if (v.x >= _qSize) v.x = _qSize - 1;
				if (v.y < 0) v.y = 0;
				if (v.y >= _qSize) v.y = _qSize - 1;
				seaSurface.addVert(sfa::Point(v.x, v.y, v.z));
				seaSurface.vertexNormals.push_back(sfa::Point(0, 0, 1));
			}
			success = success && seaSurface.AddFacet(0, 1, 2);
			success = success && seaSurface.AddFacet(0, 2, 3);
		}


		if (!success)
			std::cout << "Cannot save the mesh" << std::endl;
		else
			std::cout << "Generated mesh: " << numPoints << " points and " << numTriangles << " facets" << std::endl;
	}

}



// Returns the tile at a certain coordinate, with the tiles North and West
void FindTiles(std::vector<Tile>& tiles, const Tile& target, Tile& north, Tile& east, Tile &northEast)
{
	Latitude southLat = target.getCoordinates().low().latitude();
	Longitude westLng = target.getCoordinates().low().longitude();
	Latitude northLat = target.getCoordinates().high().latitude();
	Longitude eastLng = target.getCoordinates().high().longitude();
	for (auto& t : tiles)
	{
		auto low = t.getCoordinates().low();
		auto high = t.getCoordinates().high();
		if (low.latitude() == northLat && low.longitude()==westLng)
			north = t;
		else if (low.latitude()==southLat && low.longitude() == eastLng)
			east = t;
		else if (low.latitude()==northLat && low.longitude() == eastLng)
			northEast = t;
	}

}

int main(int argc, char **argv)
{
	cognitics::gdal::init(argv[0]);
	char const* cdb = "C:/ocb/CDB_Yemen_4.0.0";
	int LOD = 4;
	double increment = 1.0 / (2 << LOD);
	//Coordinates target(12.75, 45);
	//Coordinates target(12.75, 44.97);
	Coordinates target(12.78, 45.043); LOD = 6;
	CoordinatesRange corner(target.longitude().value()-increment, target.longitude().value()+increment*2,
		target.latitude().value()-increment, target.latitude().value()+increment*2);
	auto tiles = generate_tiles(corner,Dataset::Elevation,LOD);
	Tile center;
	for (auto& t : tiles)
		if (t.getCoordinates().low().latitude().value() <= target.latitude().value() &&
			t.getCoordinates().low().longitude().value() <= target.longitude().value() &&
			t.getCoordinates().high().latitude().value() > target.latitude().value() &&
			t.getCoordinates().high().longitude().value() > target.longitude().value())
			center = t;

	Tile north, east, northEast;
	FindTiles(tiles,center,north,east,northEast);
	DEMesher test(cdb,center);
	test.Generate(north,east,northEast);
	std::unique_ptr<scenegraph::Scene> scene(new scenegraph::Scene());

	test.Triangulate(scene.get());
	scenegraph::buildOpenFlightFromScene("c:/build/temp/terrain.flt", scene.get());


	for (double s =89; s<=90; s+=0.1)
		for (double e = 44; e <= 48.1; e += 0.1)
		{
			DEMesher  mesher(cdb,s,e,3);
			std::cout << "A: "<< mesher.TargetTile().getFilename("tif") << std::endl;
			CoordinatesRange corner(e,e+0.0001,s,s+0.0001);
			auto list = generate_tiles(corner, mesher.TargetTile().getDataset(), 3);
			std::cout << "B: "<< list[0].getFilename() << std::endl;

		}
}
