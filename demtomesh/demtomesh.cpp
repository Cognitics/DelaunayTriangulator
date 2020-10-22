// demtomesh.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "cdb_tile/Tile.h"
#include "ccl/FileInfo.h"
#include "ctl/Vector.h"
#include "ccl/gdal.h"
#include "elev/SimpleDEMReader.h"
#include "ctl/DelaunayTriangulation.h"
#include "ctl/TIN.h"
#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>
#undef min



using namespace cognitics::cdb;


//! \brief A polygon with coordinates relative to the tile grid
class TilePolygon
{
    // All the points of the polygon are enclosed within this rectangle.
    int _minx, _maxx, _miny, _maxy;
    ctl::PointList _p; // Points in the polygon 

//! \brief Defines the polygon, and convertes them relative to the tile origin
    TilePolygon(const ctl::PointList& poly, ctl::Point tileOrigin)
    {
        if (!poly.size())
        {
            _minx = _maxx = _miny = _maxy = 0;
            return;
        }
        _minx = _miny = INT_MAX;
        _maxx = _maxy = INT_MIN;
        _p.resize(poly.size());
        for(int i = 1; i < poly.size(); i++)
        {
            _p[i] = poly[i] - tileOrigin;
            if(_minx > _p[i].x)
                _minx = _p[i].x;
            if(_maxx < (int)_p[i].x)
                _maxx = (int)_p[i].x+1;
            if(_miny > _p[i].y)
                _miny = _p[i].y;
            if(_maxy < (int)_p[i].y)
                _maxy = (int)_p[i].y+1;
        }
    }
};


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

class DEMesherLineSelector;

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
        return _quantized.data() + _qSize * y; 
    }

    public: 

 //! \brief The quantized elevations are given in cm, and the low 8 bits 
 //         contain the flags POINT_xxx
    static int ElevToQuantized(double elevation)
    {
        return ((int)(elevation * 100)) <<8;
    }

    static double QuantizedToElevation(int quantized)
    {
        return (double)(quantized >> 8)/100.0;
    }

 //! \brief Quantized elevation flags
 //      To facilitate the traversal of the elevation array this flags indicate 
 //      how each point should be handled
    enum {
        POINT_BOUNDARY = 1,  // The point is the bounday of the target area
        POINT_REMOVED = 2,   // The point was discarded, considered not relevant
        POINT_EDGE = 4,      // The point was identified as an edge between removed areas, so it should be preserved
        POINT_PENDING = 8    // The point has been already pushed in the analysis stack
    };
  

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

    void Triangulate(ctl::TIN &tr, MeshWriter *save=nullptr);

private:
    void SelectBoundary(DEMesherLineSelector& processor, std::vector<int> &constraints);
    void MarkBoundary();
    bool CopyGrid(elev::SimpleDEMReader &reader, int x, int y, int subsample=2);
    void SimplifyFlat(double  toleranceMeters);    
    void SimplifySlope(double  toleranceRatio);    
    void SaveQuantized(const char* fileName);
};

/*! \brief Generic interface to tweak the selected points in a line of quantized elevations
     The line is defined by its length and the increment between values, so it can
     handle lines which are not contiguous in the array

    The clients implement the Process() method that carry out the processing. The intention is
    that the clients should only remove some dots from the triangulation, or preserve them.
*/
class DEMesherLineSelector
{
protected:
    int *_first = nullptr; // Pointer to the first quantized point
    int _n = 0;    // Number of points
    int _step = 1; // Step between points

    DEMesherLineSelector(int lineSize, int step) :
        _n(lineSize), _step(step) {}

 /*! \brief Implements simplification of the line
  
      The function must be thread safe

      The 'constraints' array contains the indexes that MUST be preserved:
         -the array is sorted
         -it is guaranteed to contain valid indexes
         -at least the first and last indexes are include [0,_n-1]
         -indexes can be duplicated

      \param first Points to the first element, the next point
      \param constraints Array of point indexes that must be preserved

       \code
         class Simplifier : public DEMesherSimplifier 
         {
            void Process(int* first, std::vector<int> &constraints);
         };
         // Lines are 65 values long, the distance between them is 513 values
         Simplifier s(65, 513);

         // Force to preserve the first two values and the last two values
         std::vector<int> constraints = {0, 1, 63, 64};
         int *firstPoint = GetFirstPoint();
         s.SimplifyLine(firstPoint, constraints);
    \endcode

 */
    virtual void Process(int* first, std::vector<int> &constraints) = 0;

    //! \brief The points are removed by setting the POINT_REMOVED flag
    static void Remove (int &value) { value |= DEMesher::POINT_REMOVED; };      // remove from triangulation
    static void Preserve (int &value) { value &= (~DEMesher::POINT_REMOVED); }; // Marks the value to keep
    static bool Removed(int value) { return (value & DEMesher::POINT_REMOVED) != 0; }

public:
    int LineSize() const { return _n; }
    int Step() const { return _step; }
    void Normalize(std::vector<int>& constraints)
    {
        std::sort(constraints.begin(), constraints.end());

        while (constraints.size() > 0 && constraints[0] < 0)
            constraints.erase(constraints.begin());
        while (constraints.size() > 0 && constraints.back() > _n - 1)
            constraints.pop_back();
        if (!constraints.size() || constraints[0] != 0)
            constraints.insert(constraints.begin(), 0);
        if (constraints.back() != _n-1)
            constraints.push_back(_n-1);
    }

    void SelectIn(int *line, int n, int step, std::vector<int> &constraints)
    {
        if (line && _n > 2)
        {
            _n = n;
            _step = step;
            Normalize(constraints);
            Process(line, constraints);
        }
    }

    void SimplifyLine(int* first, std::vector<int>& constraints)
    {
        if (first && _n > 2)
        {
            Normalize(constraints);
            Process(first, constraints);
        }
    }
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
            dst[i] = ElevToQuantized(srcRow[i * subsample]);
    }
    return true;
}


/*! \brief Removes triangulation points from flat areas

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

    The algorithm uses a simple index to address points in the quantized buffer. The 
*/
void DEMesher::SimplifyFlat(double toleranceMeters)
{
    int tolerance = ElevToQuantized(toleranceMeters);

    int* data = _quantized.data();

    int limit = _qSize * (_qSize - 1);

    std::vector<int> stack; // The stack contains the positions of the next point pending to check

    // The altitude range of the current region. We need to make sure the range does 
    // not exceed the tolerance
    int low, high;


    auto CheckNeighbor = [&](int index)->bool
    {
        assert(index >= 0 && index < _qSize*_qSize);
        int altitude = data[index] & (~0xff);
        if (altitude == 0 && tolerance != 0)
            return false;
        bool inRange = altitude>=high-tolerance && altitude<=low+tolerance;
        if (inRange)
        {
            if (low > altitude)
                low = altitude;
            else if (high < altitude)
                high = altitude;
            if ((data[index] & (POINT_BOUNDARY | POINT_REMOVED | POINT_EDGE | POINT_PENDING)) == 0)
            {
                data[index] |= POINT_PENDING;
                stack.push_back(index);
            }
        }
        else if (altitude != 0 && tolerance!=0 && data[index] & (POINT_BOUNDARY | POINT_EDGE))
        {
            // If the neighbor is an edge and the altitude is close to the limit, let's
            // discard it anyway
            int highTolerance = (tolerance * 12) / 10;
            inRange = altitude>=high-highTolerance && altitude<=low+highTolerance;
        }
            
        return inRange;
    };

    int removed = 0;
    int edges = 0;

    // Check all the point, skipping the first and last rows because they are boundaries
    for(int testIndex= _qSize; testIndex < limit; testIndex++)
    {
        int v = _quantized[testIndex];
        if (v & (POINT_BOUNDARY | POINT_REMOVED | POINT_EDGE))
            continue;
        low = high = v & (~0xff);
        // For a point at sea level the tolerance is zero so that the area only contains points at sea level
        tolerance = (low == 0) ? 0 : ElevToQuantized(toleranceMeters);
        stack.push_back(testIndex);
        while(stack.size())
        {
            int idx = stack.back();
            stack.pop_back();
            data[idx] &= ~POINT_PENDING;
            if (data[idx] & (POINT_BOUNDARY | POINT_REMOVED | POINT_EDGE))
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
                data[idx] |= POINT_REMOVED;
                removed++;
            }
            else
            {
                data[idx] |= POINT_EDGE;
                edges++;
            }
        }
    }
    std::cout << removed << " removed points, " << edges << " edge points" << std::endl;
}


void DEMesher::SimplifySlope(double tolerance)
{
    auto IsInline = [tolerance](int e1, int e2, int e3)
    {
        if ((e1 | e2 | e3) & POINT_REMOVED)
            return false;
        e2 >>= 8;
        e1 >>= 8;
        e3 >>= 8;
        if (e2 == 0)
            return e1 == 0 && e3 == 0;
        int error = e2-(e1 + e3) / 2;
        int margin = 100 + (int)(tolerance * (e3 - e1));
        if (margin < 0)
            margin = -margin;
        if (error < 0)
            error = -error;
        return error <= margin;
    };
    int removed = 0;
    int accepted = 0;
    for (int y = 0; y < _qSize; y++)
    {
        int* row = ptr(y);
        int* above = (y > 0) ? ptr(y - 1) : row;
        int* below = (y + 1 < _qSize) ? ptr(y + 1) : row;
        for (int x = 1; x < _qSize - 1; x++)
        {
            // There are 8 lines that combine the center pixel with the neighbor
            // Check if the center is aligned with any of them. There must be a better way
            if (IsInline(row[x-1],   row[x], row[x+1]) ||
                IsInline(above[x],   row[x], below[x]) ||
                IsInline(above[x-1], row[x], below[x+1]) ||
                IsInline(above[x+1], row[x], below[x-1]) ||

                IsInline(above[x-1], row[x], row[x+1]) ||
                IsInline(below[x-1], row[x], row[x+1]) ||
                IsInline(row[x-1], row[x], above[x+1]) ||
                IsInline(row[x-1], row[x], below[x+1])
                )
            {
                removed++;
                row[x] |= POINT_REMOVED;
            }
            else if ((row[x] & POINT_REMOVED) == 0)
                accepted++;
        }
    }
    std::cout << removed << " colinear points, " << accepted << " remaining points" << std::endl;
}


/*! \brief Simplifies a line keeping the error below a limit
 
    The function checks how much error would be generated if the point was removed
    and the previous and next segments were joined. If the error is less than the
    objective the point is removed and the two segments joined.

    Point N considers two segments [M,N] and [N,N+1]. The algorithm evaluates the
    error in the range [M, N+1]
*/
class LineSimplifyByError : public DEMesherLineSelector
{
    // Quantized tolerance, removing 8 bits
    int _tolerance;

    // Maximum spacing between preserved points. We want to keep a few points in the
    // line to make triangulation simpler
    static int const MAX_DISTANCE = 32;
    static int const ELEV_SHIFT = 8;

public:
    LineSimplifyByError(double toleranceMeters, int lineSize, int lineStep) :
        DEMesherLineSelector(lineSize, lineStep),
        _tolerance(DEMesher::ElevToQuantized(toleranceMeters)>>ELEV_SHIFT)
    {}

    void Process(int* first, std::vector<int>& constraints)
    {
        auto Item = [&](int x)->int& { return first[x * _step]; }; // Address of a value

        int currentConstraint = 0;
        int xprev = 0, xnext = 0; // index of the previous point preserved
        int x;

        // Calculates the error at point X, given the current segment
        auto ZError = [&](int x)->int
        {
            int zprev = Item(xprev)>>ELEV_SHIFT;
            int znext = Item(xnext)>>ELEV_SHIFT;
            // Note that the multiplication will NOT overflow because the range [xprev,xmax] is less than 32
            int z = (Item(x) >> ELEV_SHIFT);
            int error = z - (zprev + ((x - xprev) * (znext - zprev)) / (xnext - xprev));
            return error >= 0 ? error : -error;
        };

        for (x = 0; x < _n-1; x++)
        {
            if (x >= xprev + MAX_DISTANCE || x == constraints[currentConstraint])
            {
                // Start the segment here
                // Consider duplicate constraints
                while (x == constraints[currentConstraint])
                    ++currentConstraint;
                Preserve(Item(x));
                xprev = x;
                continue;
            }
            xnext = x + 1;
            // Calculate the error for all the points in the range
            int error = 0;
            for (int i = x; i > xprev && error <= _tolerance; i--)
                error = ZError(x);
            if (error <= _tolerance)
                Remove(Item(x));
            else
            {
                Preserve(Item(x));
                xprev = x;
            }
        }
        Preserve(Item(x));
    }
}; // LineSimplifyByError


/*! \brief Preserve sea level boundaries

     This class inserts grid points where there is a transition from sea to land
*/
class KeepSeaLevelInLine : public DEMesherLineSelector
{
    // Maximum and minimum sea level, quantized
    int _maxLevel, _minLevel;
    void Process(int* first, std::vector<int>& constraints)
    {
        auto InRange = [&](int x)->bool
        {
            int z = first[x * _step] >> 8;
            return (z >= _minLevel && z <= _maxLevel);
        };

        auto KeepSample = [&](int x)
        {
            Preserve(first[x * _step]);
        };

        for (int x = 1; x < _n-1; x++)
        {
            if (InRange(x))
            {
                if (!InRange(x+1))
                {
                    KeepSample(x);
                    KeepSample(x+1);

                }
                if (!InRange(x-1))
                {
                    KeepSample(x);
                    KeepSample(x-1);
                }
            }
        }
    }
public: 
    KeepSeaLevelInLine(double maxLevelMeters, double minLevelMeters,  int lineSize, int lineStep) :
        DEMesherLineSelector(lineSize, lineStep),
        _maxLevel(DEMesher::ElevToQuantized(maxLevelMeters) >> 8),
        _minLevel(DEMesher::ElevToQuantized(minLevelMeters) >> 8) {}
};


class TestLineSimple
{
public:
    void DoTest(std::vector<double> values, std::vector<int> constraints, int step)
    {
        std::vector<int> elev(values.size()*step, -10000);
        for (int i = 0; i < values.size(); i++)
            elev[i*step] = DEMesher::ElevToQuantized(values[i]) | DEMesher::POINT_BOUNDARY | DEMesher::POINT_REMOVED;
        LineSimplifyByError obj1(10, (int)elev.size(), 1);
        obj1.SimplifyLine(elev.data(), constraints);
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

TestLineSimple t;

void DEMesher::MarkBoundary()
{
    for (int i = 0; i < _qSize; i++)
    {
        ptr(0)[i] |= POINT_BOUNDARY;
        ptr(_qSize - 1)[i] |= POINT_BOUNDARY;
        ptr(i)[0] |= POINT_BOUNDARY;
        ptr(i)[_qSize - 1] |= POINT_BOUNDARY;
    }
}

void DEMesher::SelectBoundary(DEMesherLineSelector& processor, std::vector<int> &constraints)
{
    // Run the four corners in the standard order (minx, miny) (maxx, miny) (maxx,maxy) (minx,maxy)
    processor.SelectIn(ptr(0),                 _qSize, 1, constraints);
    processor.SelectIn(ptr(0)+_qSize-1,        _qSize, _qSize, constraints);
    processor.SelectIn(ptr(_qSize-1)+_qSize-1, _qSize, -_qSize, constraints);
    processor.SelectIn(ptr(_qSize - 1),        _qSize, -_qSize, constraints);
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
        _quantized.resize(_qSize * (_qSize+1), ElevToQuantized(-10000));
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
    MarkBoundary();
    double tolerance = 20;
    double seaLevel = 0;
    SimplifyFlat(tolerance);
    SimplifySlope(0.2);


    // Run the mesh simplification and the sea level detail on the boundary
    std::vector<int> empty;
    LineSimplifyByError topSimplify(tolerance / 2, _qSize, 1);
    KeepSeaLevelInLine topSeaLevel(seaLevel, seaLevel, _qSize, 1);
    SelectBoundary(topSimplify, empty);
    SelectBoundary(topSeaLevel, empty);


    SaveQuantized("c:\\build\\temp\\quantized.json");
    return;


   
    topSimplify.SimplifyLine(ptr(0), empty);
    topSeaLevel.SimplifyLine(ptr(0), empty);
    topSimplify.SimplifyLine(ptr(_qSize-1), empty);
    topSeaLevel.SimplifyLine(ptr(_qSize-1), empty);
    LineSimplifyByError leftSimplify(tolerance / 2, _qSize, _qSize);
    KeepSeaLevelInLine  leftSeaLevel(seaLevel, seaLevel, _qSize, _qSize);
    leftSimplify.SimplifyLine(ptr(0), empty);
    leftSeaLevel.SimplifyLine(ptr(0), empty);
    leftSimplify.SimplifyLine(ptr(0)+_qSize-1, empty);
    leftSeaLevel.SimplifyLine(ptr(0)+_qSize-1, empty);

}

void DEMesher::Triangulate(ctl::TIN &container, MeshWriter *save)
{
    // The boundary rectangle is counterclockwise, but on a 'normal' (Y up) axis
    // the best way to build it is: (xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax)
    int* top = ptr(0);
    int* bottom = ptr(_qSize - 1);
    ctl::PointList boundary{ 
        ctl::Point(       -1,        -1, QuantizedToElevation(top[0])),
        ctl::Point(_qSize,        -1, QuantizedToElevation(top[_qSize-1])), 
        ctl::Point(_qSize, _qSize, QuantizedToElevation(bottom[_qSize-1])), 
        ctl::Point(-1,        _qSize, QuantizedToElevation(bottom[0])) 
    };

    ctl::DelaunayTriangulation dt(boundary);

    // All the points that have not been removed are constraints

    // First, insert the constrainted outside contour
    std::vector<ctl::Point> contour;

    auto AddToContour = [&](int x, int y)
    {
        int altitude = ptr(y)[x];
        if ((altitude & POINT_REMOVED) == 0)
            contour.push_back(ctl::Point(x, y, QuantizedToElevation(altitude)));
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
            if ((altitude & POINT_REMOVED) == 0)
            {
                insertedPoints++;
                ctl::Point target(x, y, QuantizedToElevation(altitude));
                if (dt.InsertConstrainedPoint(target) == 0)
                    failedPoints++;
            }
        }
    std::vector<ctl::Edge*> edges;
    edges = dt.GatherTriangles(ctl::PointList());
    container.CreateFromDT(&dt, edges);
    if (save)
    {
        bool success = true;
        auto verts = container.verts;

        // Indexes of the corners. The last entry always contains the last enumerated corner.
        int corners[4] = { -1,-1,-1,-1 };
        int numPoints = 0;
        int numCorners = 0;
        for (auto& v : verts)
        {
            if (v.x < 0 || v.x >= _qSize)
            {
                // The corners will generate a square that covers the tile representing seal level
                // Corner altitude just above zero so that the surface hides the mesh below
                v.z = 0.001;
                if (v.x < 0) v.x = 0;
                if (v.x >= _qSize) v.x = _qSize - 1;
                if (v.y < 0) v.y = 0;
                if (v.y >= _qSize) v.y = _qSize - 1;
                if (numCorners<4)
                    corners[numCorners] = corners[3] = numPoints;
                numCorners++;
            }
            success = success && save->AddVertex(v);
            numPoints++;
        }
        // Add the two triangles that mark the sea level
        success = success && save->AddFacet(corners[0], corners[1], corners[2]);
        success = success && save->AddFacet(corners[0], corners[2], corners[3]);

        auto triangles = container.triangles;
        auto IsCorner = [&](int index)->bool { return index<=corners[3] && (index==corners[3] || index == corners[0] || index == corners[1] || index == corners[2]); };
        for (int i = 0; i+2 < triangles.size(); i += 3)
        {
            if (!IsCorner(triangles[i]) && !IsCorner(triangles[i+1]) && !IsCorner(triangles[i+2]))
                success = success && save->AddFacet(triangles[i+0], triangles[i+1], triangles[i+2]);
        }

        success = success && save->Close(); // Make sure the data is saved
        if (!success)
            std::cout << "Cannot save the mesh to file" << std::endl;
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
    ObjMeshWriter writer("c:/build/temp/terrainmesh.obj");
    ctl::TIN container;
    test.Triangulate(container, &writer);


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
