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
    {
    }
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

int badFUnction() {
    return true;
}


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


 //! \brief The quantized elevations are given in cm, and the low 8 bits 
 //         contain the flags POINT_xxx
    static int ElevToInteger(double elevation)
    {
        return ((int)(elevation * 100)) <<8;
    }

    static double IntegerToElevation(int elevation)
    {
        return (double)(elevation >> 8)/100.0;
    }

 //! \brief Quantized elevation flags
 //      To facilitate the traversal of the elevation array we need to include 
 //      information about what is in each data point
    enum {
        POINT_BOUNDARY = 1,  // The point is the bounday of the target area
        POINT_REMOVED = 2,   // The point was discarded, considered not relevant
        POINT_EDGE = 4,      // The point was identified as an edge between removed areas, so it should be preserved
        POINT_PENDING = 8    // The point has been already pushed in the analysis stack
    };

   
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

    void Triangulate(ctl::TIN &tr, MeshWriter *save=nullptr);

private:
    void MarkBoundary();
    bool CopyGrid(elev::SimpleDEMReader &reader, int x, int y, int subsample=2);
    void SimplifyFlat(double  toleranceMeters);    
    void SimplifySlope(double  toleranceRatio);    
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
            dst[i] = ElevToInteger(srcRow[i * subsample]);
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
    int tolerance = ElevToInteger(toleranceMeters);

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
        tolerance = (low == 0) ? 0 : ElevToInteger(toleranceMeters);
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


/*! \brief Given an array of quantized elevations, remove redundant points
 
    The function computes the error which would be generated by the triangulation
    by removing each point, and if 

*/
void SimplifyEdge(double toleranceMeters, int* firstInput)
{

}


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
        _quantized.resize(_qSize * (_qSize+1), ElevToInteger(-10000));
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
    SimplifyFlat(20);
    SimplifySlope(0.2);
    SaveQuantized("c:\\build\\temp\\quantized.json");
}

void DEMesher::Triangulate(ctl::TIN &container, MeshWriter *save)
{
    // The boundary rectangle is counterclockwise, but on a 'normal' (Y up) axis
    // the best way to build it is: (xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax)
    int* top = ptr(0);
    int* bottom = ptr(_qSize - 1);
    ctl::PointList boundary{ 
        ctl::Point(       -1,        -1, IntegerToElevation(top[0])),
        ctl::Point(_qSize,        -1, IntegerToElevation(top[_qSize-1])), 
        ctl::Point(_qSize, _qSize, IntegerToElevation(bottom[_qSize-1])), 
        ctl::Point(-1,        _qSize, IntegerToElevation(bottom[0])) 
    };

    ctl::DelaunayTriangulation dt(boundary);

    // All the points that have not been removed are constraints
    int failedPoints = 0;
    int insertedPoints = 0;
    for (int y=0; y<_qSize; y++)
        for (int x = 0; x < _qSize; x++)
        {
            int altitude = ptr(y)[x];
            if ((altitude & POINT_REMOVED) == 0)
            {
                insertedPoints++;
                ctl::Point target(x, y, IntegerToElevation(altitude));
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

        // Indexes of the corners. The last entry contains the last enumerated corner
        int corners[4] = { -1,-1,-1,-1 };
        int numPoints = 0;
        int numCorners = 0;
        for (auto& v : verts)
        {
            success = success && save->AddVertex(v);
            if (v.x < 0 || v.x >= _qSize)
            {
                std::cout << "Corner: " << numPoints << " " << v.x << " " << v.y<<std::endl;
                if (numCorners<4)
                    corners[numCorners] = corners[3] = numPoints;
                numCorners++;
            }
            numPoints++;
        }
        auto triangles = container.triangles;
        auto IsCorner = [&](int index)->bool { return index<=corners[3] && (index==corners[3] || index == corners[0] || index == corners[1] || index == corners[2]); };
        for (int i = 0; i + 2 < triangles.size(); i += 3)
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
