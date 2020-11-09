// demtomesh.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "cdb_tile/Tile.h"
#include "ccl/FileInfo.h"
#include "ctl/Vector.h"
#include "ccl/gdal.h"
#include "elev/SimpleDEMReader.h"
#include "ctl/DelaunayTriangulation.h"
#include "ctl/TIN.h"
#include "ctl/Util.h"
#include <scenegraph/Scene.h>
#include <scenegraphflt/scenegraphflt.h>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#undef min
#undef max



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


/*! \brief Represents a vector or array of elevation elements
    The structure points to a table of elevation values organized by rows. 
    If the structure represents a vector, only the first element of each row is used

*/
class QuantizedElevation
{
protected:
    int *_data = nullptr; // Pointer to the first quantized array
    int _n = 0;    // Number of rows
    int _step = 1; // Step between rows, could be negative. If it is an array, this is the number of columns

public:
    int rows() { return _n; }
    int cols() { return _step >= 0 ? _step : -_step; }

    int* ptr(int y, int x=0) 
    {
        if (y < 0 || y >= rows() || x<0 || x >= cols())
            return nullptr;
        return _data + y * _step + x;
    }

    int* top() { return ptr(0); }
    int* bot() { return ptr(_n-1); }

    // A few simple macros to manipulate the quantized elevations
    // Reference to an entry
    int& at(int y, int x=0) { return ptr(y,x)[0]; }


    // Return the elevation of a point, without the flags
    int Z(int y, int x=0) { return MeshQuantizedZ(at(y,x)); }

    // Return the elevation of a point, in meters
    double Elevation(int y, int x=0) { return MeshQuantizedToElevation(at(y,x)); }

    void Remove(int y, int x=0)    { at(y,x) |= MESHPOINT_REMOVED; };   // Remove from triangulation
    void Preserve(int y, int x=0) { at(y,x) &= (~MESHPOINT_REMOVED); }; // Marks the value to be used
    void KeepEdge(int y, int x=0) { at(y,x) |= MESHPOINT_EDGE; at(y,x) &= (~MESHPOINT_REMOVED); }; // Edge to be forced into the triangulation
    bool IsRemoved(int y, int x=0) { return (at(y,x) & MESHPOINT_REMOVED) != 0; }
    bool IsEdge(int y, int x=0)    { return (at(y,x) & MESHPOINT_EDGE) != 0; }
    bool IsBoundary(int y, int x=0)    { return (at(y,x) & MESHPOINT_BOUNDARY) != 0; }

    void FromVector(std::vector<int> &values, int rowSize)
    {
        rowSize = std::max(rowSize, 0);
        _data = values.data();
        _n = (int)values.size() / rowSize;
        _step = rowSize;
    }
};



/*! \brief Generic interface to tweak the selected points in quantized elevations
    
Derived classes implement the Process() method that carry out the processing. The intention is
that the should only remove some elements from the triangulation, or preserve them.

\code
        class Simplifier : public DEMesher::GridSelector 
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
class GridSelector : public QuantizedElevation
{
protected:
    virtual ~GridSelector() {}

/*! \brief Implements simplification of the line

    The fields _first, _n and _step contain the location of tha target line
*/
    virtual void Process() = 0;

public:
//! \brief Sends a line (defined by the start and the step interval) to a line selector
    void Select(int* line, int n, int step)
    {
        if (!line || n <= 0)
            return;
        _n = n;
        _data = line;
        _step = step;
        Process();
    }
};

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
    QuantizedElevation _elev;

    // Some configuration parameters: the color for the faces
    int _faceRGB[3] = { 255, 255, 255 };

    // Add an additional face with the sea level
    bool _addSeaLevel = false;

    // Number of grid units per tile. It depends on the LOD
    int _nominalTileSize = 1024;

    // Position of the SW corner of this tile relative to the origin. This is
    // given in tile units, before scaling
    ctl::Point _origin = ctl::Point(0,0,0);

    // Scaling for the three dimensions to convert from tile units to meters. 
    // It depends on the flat projection, so it varies per tile.
    ctl::Point _scale = ctl::Point(1.0, 1.0, 1.0);

    std::vector<ctl::PointList*> _polygons;

    bool IsGenerated() const { return _quantized.size() > 0; }

public: 

 //! \brief Creates a mesher for the tile at a certain latitude and longitude
    DEMesher(const std::string& cdbRoot, const Tile& t, const ctl::Point& origin, const ctl::Point& scale);
    DEMesher(const DEMesher &) = delete;
    DEMesher& operator =(const DEMesher&) = delete;

    virtual ~DEMesher();


    bool AddPolygon(ctl::PointList& p);


    bool Generate();

    void Triangulate(scenegraph::Scene* saved);

private:
    bool CheckForCDB();
    void SelectBoundary(GridSelector& processor);
    bool CopyGrid(const Tile &source, int x, int y);
    void SaveQuantized(const char* fileName);

};

//! \brief Calculates distance from a geographical coordinate to the SW corner of the tile, in tile units
ctl::Point ConvertGeoToTileUnits(LOD lod, const Coordinates &tileOrigin, Coordinates geo, double elevation = 0)
{
    int tileGridSize = lod.dimensions();
    double tileStep =  get_tile_width(TileLatitude(tileOrigin.latitude()));

    double xUnits = (geo.longitude().value() - tileOrigin.longitude().value()) / tileStep * lod.cols() * tileGridSize;
    double yUnits = (geo.latitude().value() - tileOrigin.latitude().value()) * lod.rows() * tileGridSize;
    return ctl::Point(xUnits, yUnits, elevation);
}


//! \brief Sets the name and location of the scene that represents the tile
void SetSceneAttributes(scenegraph::Scene* scene, const std::string &name, const CoordinatesRange &location);


/*! \brief Cache of tile data.
    
    To complete the tile each mesher needs to read the three surrounding tiles. To avoid reading 
    the DEM files over and over, this cache keeps the data as it has been read. 
*/
class DEMCache
{
public:
    static int _tileSize;

    // To save a lot of memory the cache entry discards the grid when it is no
    // longer needed and leaves only four lines with the boundaries:
    //   -top, bottom, left, right

    static std::vector<int> *GetGrid(const Tile &tile, const std::string& cdbRoot);
    static void FreeGrid(const Tile& tile);

private:
    static std::vector<int> *ReplaceGrid(const Tile &tile, const std::string& cdbRoot);
    
};

int DEMCache::_tileSize = 512;
std::map<std::string, std::vector<int> *> elevationCache;

void DEMCache::FreeGrid(const Tile &tile)
{
    auto range = tile.getFilename();
    auto found = elevationCache.find(range);
    if (found != elevationCache.end())
    {
        delete found->second;
        elevationCache.erase(range);
    }
}
//! \brief Generates a replacement for the elevation grid from a tile with lower LOD
std::vector<int>* DEMCache::ReplaceGrid(const Tile& tile, const std::string& cdbRoot)
{
    // Finding the parent tile is not a trivial affair
    LOD newLOD(tile.getLod() - 1);
    if (newLOD.dimensions() <= 4)
        return nullptr;

    auto higherTile = generate_tiles(tile.getCoordinates(), tile.getDataset(), newLOD);
    assert(higherTile.size() == 1);
    auto highGrid = GetGrid(higherTile[0], cdbRoot);
    if (!highGrid)
        return nullptr;
    // Now, take one of the four quadrants, and resample by interpolation
    int xoffset = 0, yoffset = 0;
    if (higherTile[0].getCoordinates().low().latitude() == tile.getCoordinates().low().latitude())
        yoffset = _tileSize / 2; // Quadrant at the bottom
    if (higherTile[0].getCoordinates().low().longitude() < tile.getCoordinates().low().longitude())
        xoffset = _tileSize / 2; // Quadrant to the right
    std::cout << "Replace with: " << higherTile[0].getFilename(".tif")<<"offsets: " << xoffset<< " "<< yoffset << std::endl;

    std::vector<int>* grid = new std::vector<int>(_tileSize * _tileSize);

    for (int x = 0; x < _tileSize; x += 2)
        for (int y = 0; y < _tileSize; y += 2)
        {
            int srcIndex = ((y/2)+yoffset) * _tileSize + x/2 + xoffset;
            int left = (*highGrid)[srcIndex];
            int right = left;
            int bottom = left;
            if (x/2+xoffset+1 < _tileSize)
                right = (*highGrid)[srcIndex + 1];
            if (y/2+yoffset+1 < _tileSize)
                bottom = (*highGrid)[srcIndex + _tileSize];
            if (left == 0 && right > 100000)
                std::cout << "Check out\n";
            int dstIndex = y * _tileSize + x;
            (*grid)[dstIndex] = left;
            (*grid)[dstIndex+1] = (left+right) / 2;
            (*grid)[dstIndex+_tileSize] = (left+bottom) / 2;
            (*grid)[dstIndex+_tileSize+1] = (right+bottom) / 2;
        }
    return grid;
}


std::vector<int> * DEMCache::GetGrid(const Tile& tile, const std::string& cdbRoot)
{
    auto range = tile.getFilename(".tif");
    auto found = elevationCache.find(range);
    if (found != elevationCache.end())
    {
        std::vector<int>* grid = found->second;
        if (grid->size() == _tileSize * _tileSize)
            return grid;
        FreeGrid(tile);
    }
    std::string fileName = cdbRoot+range;
    if (!ccl::fileExists(fileName))
    {
       auto newGrid =  ReplaceGrid(tile, cdbRoot);
       if (newGrid)
           elevationCache[range] = newGrid;
       return newGrid;
    }

    OGRSpatialReference oSRS; // GDAL reference system
    oSRS.SetWellKnownGeogCS("WGS84");
    elev::SimpleDEMReader reader(fileName, oSRS);
    if (!reader.Open())
        return nullptr;

    int width = reader.getWidth();
    int height = reader.getHeight();
    if (width != height || width<_tileSize || height < _tileSize)
        return nullptr;
    double scaleFactor = (float)_tileSize / (float)width;
    reader.setScaleFactor((float)_tileSize/(float)width);

    width = reader.getScaledWidth();
    height = reader.getScaledHeight();
    //REVISIT: Can we have rounding errors
    if (width < _tileSize || height < _tileSize)
        return nullptr;


    std::vector<double> originalGrid;
    reader.getGrid(originalGrid);
    std::vector<int>* grid = new std::vector<int>(_tileSize * _tileSize);
    for (int y = 0; y < _tileSize; y++)
    {
        for (int x = 0; x < _tileSize; x++)
            (*grid)[y*_tileSize + x] = MeshElevationToQuantized(originalGrid[y*width + x]);
    }
    elevationCache[range] = grid;
    return grid;
}

DEMesher::DEMesher(const std::string& cdbRoot, const Tile& t, const ctl::Point &origin, const ctl::Point& scale) :
        _root(cdbRoot),
        _tile(t),
        _origin(origin),
        _scale(scale),
        _nominalTileSize(LOD(t.getLod()).dimensions()) {}

DEMesher::~DEMesher()
{
    for (ctl::PointList* p : _polygons)
        delete p;
}

bool DEMesher::CheckForCDB()
{
    _root = ccl::standardizeSlashes(_root);
    if(!_root.size() || _root.back()!='/')
        _root += "/";
    ccl::FileInfo fi;
    return ccl::directoryExists(_root);
}



//! \brief Includes a polygon in the mesh6
// The polygon is given in geographic coordinates. For now the altitude is ignored.
bool DEMesher::AddPolygon(ctl::PointList& p)
{
    // Convert the polygon to tile units
    ctl::PointList *converted = new ctl::PointList(p.size());
    for (int i = 0; i < p.size(); i++)
    {
        Coordinates geo(p[i].y, p[i].x);
        converted[0][i] = ConvertGeoToTileUnits(_tile.getLod(), _tile.getCoordinates().low(), geo);
        converted[0][i].y = 1024 - converted[0][i].y;
    }
    ctl::PointList tightBoundary =
    {
        ctl::Point(0, 0, 0), ctl::Point(_nominalTileSize, 0, 0), ctl::Point(_nominalTileSize, _nominalTileSize, 0),  
                             ctl::Point(0, _nominalTileSize, 0), ctl::Point(0, 0, 0)
    };
    ctl::PointList clipped = ctl::ClipToPolygon(*converted, tightBoundary, 1.0/(_nominalTileSize*2));
    if (clipped.size() <= 2)
    {
        delete converted;
        return false;
    }
    converted->resize(clipped.size());
    std::copy(clipped.begin(), clipped.end(), converted->begin());
    _polygons.push_back(converted);
    return true;
}


//! \brief Writes a JSON with the quantized values, for analysis outside
void DEMesher::SaveQuantized(const char* fileName)
{
    std::ofstream fd(fileName, std::ios_base::out | std::ios_base::binary);
    fd << "{ \"data\": [\n";
    for(int y = 0; y < _qSize; y++)
    {
        int* src = _elev.ptr(y);
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
bool DEMesher::CopyGrid(const Tile &source, int x, int y)
{
    if (source.getDataset().code() == Dataset::Invalid)
        return false;
    auto grid = DEMCache::GetGrid(source, _root);
    if (!grid)
    {
        std::cout << "Missing: " << source.getFilename(".tif") << std::endl;
        return false;
    }
    std::cout << "File: " << source.getFilename(".tif") << std::endl;

    int gridSize = _qSize - 1;
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
    int copyWidth = std::min(_qSize - x, gridSize - srcx);
    for(; y<_qSize && srcy<gridSize; y++, srcy++)
    { 
        int* dst = _elev.ptr(y, x);
        int* srcRow = (*grid).data() + gridSize *srcy + srcx;
        for(int i = 0; i < copyWidth; i++)
            dst[i] = srcRow[i];
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
class SimplifyFlat : public GridSelector
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
        if (_step <=1 || _n<=1)
            return;

        int limit = rows()*cols()-1;

        // The altitude range of the current region. We need to make sure the range does 
        // not exceed the tolerance
        int low, high;

        auto CheckNeighbor = [&](int index)->bool
        {
            assert(index >= 0 && index < rows() * cols());
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
                // If the neighbor is edge discard the point anyway
                inRange = true;
            }

            return inRange;
        };

        int removed = 0;
        int edges = 0;

        // Check all the point, skipping the first and last rows because they are boundaries
        for (int testIndex = cols(); testIndex < limit; testIndex++)
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
                if (CheckNeighbor(idx + _step)) // Below neighbor
                    inRange++;
                if (CheckNeighbor(idx - _step)) // Above neighbor
                    inRange++;
                if (inRange == 4 || (inRange==3))// && MeshQuantizedZ(_data[idx]!=0)))
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
class SimplifySlope : public GridSelector
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
        if (_n <= 1 || _step <= 1)
            return;
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
        for (int y = 1; y < rows()-1; y++)
        {
            int* row = ptr(y);
            int* above = ptr(y - 1);
            int* below = ptr(y + 1);
            for (int x = 1; x < cols() - 1; x++)
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
                {
                    accepted++;
                    // Accepted pixels at sea level are marked as edges, so that they are not 
                    // simplified out
                    if (MeshQuantizedZ(row[x]) == 0)
                        row[x] |= MESHPOINT_EDGE;
                }
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
class LineSimplifyByError : public GridSelector
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
            int z = Z(x);
            int z1 = Z(x-1);
            int z2 = Z(x+1);
            // Preserve peaks and valleys
            if ((z > z1 && z > z2) || (z < z1 && z < z2))
                KeepEdge(x);

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


/*! \brief Preserve sea level boundaries in a line

     This class keeps the points in the line where there is a transition from sea to land
*/
class KeepSeaLevelInLine : public GridSelector
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
class BoundaryMarker : public GridSelector
{
    void Process()
    {
        for (int x = 0; x < _n; x++)
            at(x) |= MESHPOINT_BOUNDARY;
        at(0) |= MESHPOINT_EDGE;
        at(_n-1) |= MESHPOINT_EDGE;
    }
};

void DEMesher::SelectBoundary(GridSelector& processor)
{
    processor.Select(_elev.ptr(0),                 _qSize, 1);
    processor.Select(_elev.ptr(0),                 _qSize, _qSize);
    processor.Select(_elev.ptr(_qSize-1),          _qSize, 1);
    processor.Select(_elev.ptr(0, _qSize-1),       _qSize, _qSize);
}

bool DEMesher::Generate()
{
    if (!CheckForCDB())
        return false;

    auto fileName = _root+_tile.getFilename(".tif");
    OGRSpatialReference oSRS; // GDAL reference system
    oSRS.SetWellKnownGeogCS("WGS84");
    Coordinates sw = _tile.getCoordinates().low();
    Coordinates ne = _tile.getCoordinates().high();
    Coordinates nw (ne.latitude().value(), sw.longitude().value());
    Coordinates se (sw.latitude().value(), ne.longitude().value());
    
    {
        // Generate the grid at half the source resolution, one additional as guard
        _qSize = DEMCache::_tileSize+1;
        _quantized.resize(_qSize*_qSize+1, MeshElevationToQuantized(-10000));
        _elev.FromVector(_quantized, _qSize);
        // Center reader skips the NORTH edge
        if (!CopyGrid(_tile, 0, 1))
        {
            _quantized.resize(0);
            _elev.FromVector(_quantized, _qSize);
            return false;
        }
;
        // Copy the missing top and left edges
        for (int i = 0; i < _qSize; i++)
        {
            _elev.ptr(0)[i] = _elev.ptr(1)[i];
            _elev.ptr(i)[_qSize - 1] = _elev.ptr(i)[_qSize - 2];
        }
        _elev.ptr(0)[_qSize - 1] = _elev.ptr(0)[_qSize - 2];

        // Copy the top-most row
        auto northTiles = generate_tiles(CoordinatesRange(nw,nw), _tile.getDataset(), _tile.getLod());
        CopyGrid(northTiles[0], 0, -(_qSize-2));

        // Copy the right-most column
        auto eastTiles = generate_tiles(CoordinatesRange(se,se), _tile.getDataset(), _tile.getLod());
        CopyGrid(eastTiles[0], _qSize-1, 1);

        // Copy the top-right corner
        auto northEastTiles = generate_tiles(CoordinatesRange(ne,ne), _tile.getDataset(), _tile.getLod());
        CopyGrid(northEastTiles[0], _qSize-1, -(_qSize-2));
    }
    
    BoundaryMarker marker;
    SelectBoundary(marker);

    int minHeight = _elev.Z(0);
    int maxHeight = minHeight+1;
    for (int y=1; y < _qSize-1; y++)
        for (int x=1; x<_qSize-1; x++)
        {
            int z = _elev.Z(y, x);
            maxHeight = std::max(z, maxHeight);
            minHeight = std::min(z, minHeight);
        }

    // If the terrain is very shallow force a smaller tolerance so that the
    // relief is at least a little visible
    double minimumHeightTolerance = (maxHeight-minHeight)/400.0;

    // The process to simplify the mesh uses a height tolerance in meters.
    // It is set at 2% the vertical dimension of the tile. 
    double heightTolerance = _scale.y*_nominalTileSize/50.0;

    double seaLevel = 0;
    // Run the mesh simplification and the sea level detail on the boundary
    LineSimplifyByError simplify(heightTolerance / 2);
    KeepSeaLevelInLine forceSeaLevel(seaLevel, seaLevel);
    SelectBoundary(simplify);
    SelectBoundary(forceSeaLevel);

    SimplifyFlat mesher(std::min(heightTolerance, minimumHeightTolerance));
    mesher.Select(_quantized.data(), _qSize, _qSize);
    SimplifySlope flattener(0.2);
    flattener.Select(_quantized.data(), _qSize, _qSize);
    return true;
}

void SetSceneAttributes(scenegraph::Scene* scene, const std::string &name, const CoordinatesRange &location)
{
    if (!scene)
        return;
    scene->name = name;
    scene->attributes.setAttribute("terrain", true);
    auto sw = location.low();
    auto ne = location.high();
    scene->attributes.setAttribute("origin_lat", sw.latitude().value());
    scene->attributes.setAttribute("origin_lon", sw.longitude().value());
    scene->attributes.setAttribute("sw_lat", sw.latitude().value());
    scene->attributes.setAttribute("sw_lon", sw.longitude().value());
    scene->attributes.setAttribute("ne_lat", ne.latitude().value());
    scene->attributes.setAttribute("ne_lon", ne.longitude().value());
}


/*! \brief Saves the triangulated mesh into the scene

    The PRESERVED points in the DEM elevation grid make a triangle mesh. The boundary
    of the grid is a constraint so that all the edges are included.

    By default, the The (X,Y) coordinates are in tile units, so they range from 0 to 1024. 
    The origin is the SW corner of the tile. Displacement and scaling are provided in the
    _origin and _scale members
    
    The mesh becomes a Face in the scene. 
*/
void DEMesher::Triangulate(scenegraph::Scene* scene)
{
    if (!IsGenerated())
        return;

    // The boundary rectangle is clockwise on a 'reversed' (Y down) axis
    // the best way to build it is: (xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax)
    // The boundary is outside the tile. The Delaunay code does not like when the 
    // boundary falls on the tile edged.
    ctl::PointList boundary {
        ctl::Point(-1,       -1, 0),   
        ctl::Point(_qSize,   -1, 0),  
        ctl::Point(_qSize,  _qSize, 0),
        ctl::Point(-1,      _qSize, 0),  
        ctl::Point(-1,      -1, 0)  };

    // All the points that have not been removed are constraints

    // First, insert the constrainted outside contour
    std::vector<ctl::Point> contour;

    auto AddToContour = [&](int x, int y)
    {
        if (!_elev.IsRemoved(y,x))
            contour.push_back(ctl::Point(x,y,_elev.Elevation(y ,x)));
    };
    int x=0, y =_qSize-1;
    for (; x < _qSize - 1; x++)
        AddToContour(x, y);
    for (; y>0; y--)
        AddToContour(x, y);
    for (; x>0; x--)
        AddToContour(x, y);
    // Note that the first point is repeated at the end
    for (; y <= _qSize-1; y++)
        AddToContour(x, y);

    ctl::DelaunayTriangulation dt(boundary, 200, 1e-6, 3e-5, 10000, 0);
       // ctl::DelaunayTriangulation::INTERPOLATE_EDGES);//| ctl::DelaunayTriangulation::INTERPOLATE_FACES);

    int failedPoints = 0;
    if (!dt.InsertConstrainedLineString(contour)) 
        failedPoints = (int)contour.size();
    int insertedPoints = (int)contour.size();
    for (y=1; y<_qSize-1; y++)
        for (x = 1; x < _qSize-1; x++)
            if (!_elev.IsRemoved(y,x))
            {
                insertedPoints++;
                auto p = ctl::Point(x, y, _elev.Elevation(y, x));
                dt.InsertWorkingPoint(p);
            }
    int subsample = _nominalTileSize / (_qSize-1);
    for (auto& polygon : _polygons)
    {
        ctl::PointList clipped = *polygon;

        for (auto& point : clipped)
        {
            point.x /= subsample;
            point.y /= subsample;
            point.z = _elev.Elevation(int(point.y + 0.5), int(point.x + 0.5))+20;
        }
        dt.InsertConstrainedPolygon(clipped);
    }

    std::vector<ctl::Edge*> edges;
    edges = dt.GatherTriangles(ctl::PointList());
    ctl::TIN container(&dt, edges);
    int useTexture = true;
    if (scene)
    {
        SetSceneAttributes(scene, _tile.getFilename(""), _tile.getCoordinates());
        scene->hasVertexNormals = true;
        // The sfa:Matrix is very well thought out
        scene->matrix.PushTranslate(sfa::Point(_origin.x, _origin.y, _origin.z));
        scene->matrix.PushScale(sfa::Point(_scale.x, _scale.y, _scale.z));
        scene->faces.emplace_back();
        scenegraph::Face& face = scene->faces.back();
        face.primaryColor = scenegraph::Color((double)_faceRGB[0]/255.0, (double)_faceRGB[1]/255.0,(double)_faceRGB[2]/255.0);
        face.alternateColor = face.primaryColor;

        Tile imageTile = _tile;
        imageTile.setDataset(Dataset::Imagery);
        std::string imageName = ccl::joinPaths(_root, imageTile.getFilename(".jp2"));

        scenegraph::MappedTexture* texture = nullptr;
        
        if (useTexture && ccl::fileExists(imageName))
        {
            std::cout << "Image tile " << _nominalTileSize << ": " << imageName << std::endl;
            face.textures.emplace_back();
            texture = &(face.textures.back());
            texture->SetTextureName(imageName);
        }

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
            if ((v.x<=0 || v.x>=_qSize) && (v.y<=0 || v.y>=_qSize))
            {
                if (numCorners<4)
                    corners[numCorners] = corners[3] = numPoints;
                numCorners++;
            }
            // The tile coordinate system has the origin on the bottom-right, but the images have the
            // origin at the top-left. Flip the Y
            v.x *= subsample;
            v.y = (_qSize-1-v.y)*subsample;
            face.addVert(sfa::Point(v.x, v.y, v.z));
            face.vertexNormals.push_back(sfa::Point(-normals[i].x, -normals[i].y, -normals[i].x));
            if (texture)
            {
                // UV is tricky: uv = 0.5/1024 samples the first row or column of the image.
                // But the image has only 1024 pixels, so 1024.5/1024 would be out of range
                // So we use as the limits 0.5/1025 and 1024.5/1025
                texture->uvs.push_back(sfa::Point((v.x+0.5) / (_nominalTileSize+1), (v.y+0.5) / (_nominalTileSize+1)));
            }
            numPoints++;
        }

        int numTriangles = 0;
        auto triangles = container.triangles;
        auto IsCorner = [&](int index)->bool { return index<=corners[3] && (index==corners[3] || index == corners[0] || index == corners[1] || index == corners[2]); };
        for (size_t i = 0; i+2 < triangles.size(); i += 3)
        {
            if (!IsCorner(triangles[i]) && !IsCorner(triangles[i+1]) && !IsCorner(triangles[i+2]))
            {
                // Note that we swap two vertices to get the triangle orientation correct, because 
                // we inverted the Y coordinate above.
                success = success && face.AddFacet(triangles[i+0], triangles[i+2], triangles[i+1]);
                ++numTriangles;
            }
        }
        
        // Add a new face with blue color that marks the sea level
        if (numCorners >= 4 && _addSeaLevel)
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
                seaSurface.addVert(sfa::Point(v.x*subsample, v.y*subsample, v.z));
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

std::vector<Tile> TestTileEnum(int lod, double west, double east, double south, double north)
{
    CoordinatesRange selector(west, east, south, north);

    auto tiles = generate_tiles(selector, Dataset::Elevation, lod);
    std::sort(tiles.begin(), tiles.end());
    for (auto t : tiles)
        std::cout << t.getFilename() << std::endl;
    assert(tiles.size() > 0);
    double southInTile = tiles[0].getCoordinates().low().latitude().value();
    double northInTile = tiles[0].getCoordinates().high().latitude().value();
    double westInTile = tiles[0].getCoordinates().low().longitude().value();
    LOD level(lod);

    for (int i = 0; i < tiles.size(); i++)
    {
        Tile current = tiles[i];
        auto sw = current.getCoordinates().low();
        auto ne = current.getCoordinates().high();

        double cellHeight = 1.0 / level.rows();
        double cellWidth = (double)get_tile_width(TileLatitude(sw.latitude().value())) / level.cols();


        assert((ne.latitude().value() - sw.latitude().value()) == cellHeight);
        assert((ne.longitude().value() - sw.longitude().value()) == cellWidth);

        assert(southInTile < north || (south == north));
        if (current.getCoordinates().low().latitude().value() != southInTile)
        {
            assert(westInTile >= west);
            assert(westInTile >= east);
            
            westInTile = sw.longitude().value();
            southInTile = sw.latitude().value();
            assert(southInTile == northInTile);
            northInTile = ne.latitude().value();
        }
        assert(sw.longitude().value() == westInTile);
        assert(westInTile < east || (east == west));
        assert(ne.longitude().value() >= west);
        westInTile = ne.longitude().value();
    }
    return tiles;
}

void AddTestPolygon(DEMesher* object, CoordinatesRange targetArea)
{
    // Start with a polygon horizontal, and 0.0002 degrees high, and then
    // apply rotations
    Coordinates sw = targetArea.low();
    Coordinates ne = targetArea.high();
    double latitude = (sw.latitude().value() + ne.latitude().value()) / 2;
    double longitude = (sw.longitude().value() + ne.longitude().value()) / 2;
    double width = ne.longitude().value() - sw.longitude().value();
    double height = 0.001;
    std::vector <sfa::Point> centerPoints(5);
    centerPoints[0] = sfa::Point(-width/2, -height/2);
    centerPoints[1] = sfa::Point(width/2, -height/2);
    centerPoints[2] = sfa::Point(width/2, height/2);
    centerPoints[3] = sfa::Point(-width/2, height/2);
    centerPoints[4] = centerPoints[0];


    ctl::PointList points(5);
    for (int r = 0; r < 4; r++)
    {
        double rotationAngle = 30.0 / 180.0 *  M_PI * r;
        sfa::Matrix rotation;
        rotation.PushRotate(sfa::Point(0, 0, 1), rotationAngle);
        rotation.PushTranslate(sfa::Point(longitude, latitude));
        for (int i = 0; i < points.size(); i++)
        {
            sfa::Point rotated = rotation * centerPoints[i];
            points[i] = ctl::Point(rotated.X(), rotated.Y(), rotated.Z());
        }
        object->AddPolygon(points);
    }
}

void TestPolygons()
{
    ctl::PointList polygon = { ctl::Point(0,0,1), ctl::Point(1,0,0), ctl::Point(0,1,0) };
    ctl::PointList extended = polygon;
    ctl::PointList line(polygon.begin(), polygon.begin()+2);
    ctl::PointList point(1, polygon[0]);
    ctl::PointList empty;
    ctl::Point outside(4, 0, 0);
    ctl::Point center = polygon[0] + polygon[1] + polygon[2];
    center.x /= 3.0;
    center.y /= 3.0;
    center.z /= 3.0;
    extended.push_back(extended.front());
    
    ctl::PointList clip = polygon;
    for (auto& point : clip)
        point = point + ctl::Point(0.2, 0.2, 0.2);
    ctl::PointList clipExtended = clip;
    clipExtended.push_back(clipExtended.front());
    

    double  area = PArea2D(polygon);
    assert(area == -0.5); 
    assert(area == PArea2D(extended));

    assert(PArea2D(line)==0);
    assert(PArea2D(point)==0);
    assert(PArea2D(empty)==0);

//    double  area3 = PArea3D(polygon);
//    assert(area3 == PArea3D(extended));


    assert(PointInPolygon(center, polygon));
    assert(PointInPolygon(center, extended));

    assert(!PointInPolygon(outside, polygon));
    assert(!PointInPolygon(outside, extended));

    double epsilon = 1e-5;
    auto section = ClipToPolygon(polygon, clip, epsilon);
    auto section2 = ClipToPolygon(polygon, clipExtended, epsilon);
    assert(section == section2);

    std::reverse(clip.begin(), clip.end());
    std::reverse(clipExtended.begin(), clipExtended.end());
    section = ClipToPolygon(polygon, clip, epsilon);
    section2 = ClipToPolygon(polygon, clipExtended, epsilon);
    assert(section == section2);
}

int main(int argc, char **argv)
{
    if (0)
    {
        TestPolygons();
       auto t1 = TestTileEnum(0, 1, 1, 0, 2);
       assert(t1.size() == 2);
       auto t2 = TestTileEnum(1,
           t1[0].getCoordinates().low().longitude().value(),
           t1[0].getCoordinates().high().longitude().value(),
           t1[0].getCoordinates().low().latitude().value(),
           t1[0].getCoordinates().high().latitude().value());
       assert(t2.size() == 4);

        TestTileEnum(0, 1, 1, 48, 52);
        TestTileEnum(1, 1, 2, 48, 52);
        TestTileEnum(1, 1, 1, 0, 2);
        TestTileEnum(1, 1, 1.5, 0, 2);
        //TestTileEnum(0, 179, -179, 0, 2);
    }

    cognitics::gdal::init(argv[0]);
    char const* cdb = "C:/ocb/CDB_Yemen_4.0.0";
    LOD lod = 6;
    double increment = 2.0 / lod.rows();
    //Coordinates target(12.75, 45);
    //Coordinates target(12.75, 44.97);
    Coordinates target(12.78, 45.0);
    
    //Coordinates target(12.98, 45.0);
    CoordinatesRange corner(target.longitude().value()-increment, target.longitude().value()+increment,
        target.latitude().value()-increment, target.latitude().value()+increment);
    auto tiles = generate_tiles(corner, Dataset::Elevation, lod);

    CoordinatesRange swCorner(target.longitude().value(), target.longitude().value()+0.0000001,
        target.latitude().value(), target.latitude().value()+0.0000001);

    auto centerTiles = generate_tiles(swCorner, Dataset::Elevation, lod);
    std::unique_ptr<scenegraph::Scene> scene(new scenegraph::Scene());
    scene.get()->name = centerTiles[0].Filename();
    
    auto fullArea = centerTiles[0].getCoordinates();
    for (auto& t : tiles )
    {
        fullArea.Expand(t.getCoordinates().low());
        fullArea.Expand(t.getCoordinates().high());
    }

    auto sw = fullArea.low();
    auto ne = fullArea.high();
    cts::FlatEarthProjection projection(sw.latitude().value(), sw.longitude().value());
    double xMeters = projection.convertGeoToLocalX(ne.longitude().value());
    double yMeters = projection.convertGeoToLocalY(ne.latitude().value());

    auto distance = ConvertGeoToTileUnits(lod, sw, ne);
    auto sceneScale = ctl::Point(xMeters / distance.x, yMeters / distance.y, 1);
    Coordinates sceneOrigin(projection.convertLocalToGeoLat(yMeters/2),
        projection.convertLocalToGeoLon(xMeters/2));
    auto ts_start = std::chrono::steady_clock::now();

    for (int i=0; i<tiles.size(); i++)
    {
//        if (i != 10 && i != 6)  continue;
        Tile& current = tiles[i];
        // Calculate distance to the origin of the new tile
        ctl::Point offset = ConvertGeoToTileUnits(lod, sceneOrigin, current.getCoordinates().low());
        DEMesher test(cdb, current, offset, sceneScale);
        if (!test.Generate())
            continue;
        AddTestPolygon(&test, fullArea);

        std::cout << "Origin: " << offset.x << "  " << offset.y << std::endl;
        scenegraph::Scene* tileScene = new scenegraph::Scene(scene.get());
        test.Triangulate(tileScene);
    }
    auto ts_stop = std::chrono::steady_clock::now();
    double time = std::chrono::duration<double>(ts_stop - ts_start).count();
    std::cout<< ("runtime: " + std::to_string(tiles.size()/time)+" tiles/second")<<std::endl;

    std::string location = 
        std::to_string(fullArea.low().latitude().value()) + "_" +
        std::to_string(fullArea.low().longitude().value()) + "_" +
        std::to_string(fullArea.high().latitude().value()) + "_" +
        std::to_string(fullArea.high().longitude().value());
    SetSceneAttributes(scene.get(), "SC_"+location, fullArea);

    std::string outputFileName = "c:/build/temp/terrain" + std::to_string(lod.value()) + ".flt";

    std::cout << "File: " << outputFileName << " " << location << std::endl;
    scenegraph::buildOpenFlightFromScene(outputFileName, scene.get());
    scene.reset();
    return 0;
}
