/****************************************************************************
Copyright 2020, Cognitics Inc.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
****************************************************************************/
#define DELAUNAY_TRIANGULATION_EXPORTS
#include "ctl.h"
#include "DelaunayTriangulationWrapper.h"
#include <cstddef>
#include <fstream>
#include <iostream>
#include <locale>
#include <codecvt>

#ifdef WIN32
#define WRAPPER_API __declspec(dllexport)
#else
#include <locale>
#include <codecvt>
#define WRAPPER_API __attribute__((visibility("default")))
#endif

//! \brief Builds a PointList from the raw vector
ctl::PointList BuildPointList(double const* points, int begin, int end)
{
    ctl::PointList pList;
    if (points && end > begin && begin >= 0)
    {
        int const VECTORSIZE = 3;
        pList.resize(end - begin);
        points += (size_t)begin * VECTORSIZE;
        for (int i = 0; i < (end - begin); i++, points += VECTORSIZE)
            pList[i] = ctl::Vector(points[0], points[1], points[2]);
    }
    return pList;
}

class DelaunayWrapper
{
    int settings_;
public:
    enum Error
    {
        INVALID_PARAMETER       = 0x04,			// A parameter in the last call had invalid value
    };

    enum Option
    {
        REMOVE_HOLES = 0x1000,          // Remove the triangles that are inside the holes
        REMOVE_EXTERIOR = 0x2000        // Remove the triangles that are outside the exterior
    };

    // The wrapper remembers the "holes" added to the triangulation so that it can remove the triangles
    // that are within the holes. See the CallGatherTriangle function
    // The first point list is the outside contour.
    std::vector<ctl::PointList> polygons_;
    ctl::PointList boundary_;

    ctl::DelaunayTriangulation tr_;
    ctl::TIN store_;

    DelaunayWrapper(const ctl::PointList& boundary, int resizeIncrement, double epsilon, double areaEpsilon, int maxEdgeFlips, int settings) :
        boundary_(boundary),
        settings_(settings),
        tr_(boundary, resizeIncrement, epsilon, areaEpsilon, maxEdgeFlips, settings&0xFFF) {}

    bool RemovingHoles() const { return (settings_ & REMOVE_HOLES) != 0; }
    bool RemovingExterior() const { return (settings_ & REMOVE_EXTERIOR) != 0; }

    //! \brief Builds a PointList from a point vector and adds it to the list of known polygons
    ctl::PointList& AddPolygon(const double* points, int begin, int end)
    {
        polygons_.emplace_back(BuildPointList(points, begin, end));
        return polygons_.back();
    }

    //! \brief Checks if a point is inside any of the known holes
    bool PointInHole(const ctl::Point& centroid)
    {
        for (int p=1; p<polygons_.size(); p++)
            if (ctl::PointInPolygon(centroid, polygons_[p]))
                return true;
        return false;
    }

    //! \brief generates a JSON file with the results of the triangulation
    void DumpToFile(std::ofstream &fd)
    {
        auto Label = [&](std::string attribute)->std::string 
        {
            return "\""+attribute+"\" : ";
        };

        fd << "{\n";
        fd.precision(4);
        std::string options;

        fd << Label("settings") << "[ "<< settings_ << " ],\n";
        fd << Label("boundary");
        DumpVector(fd, boundary_, 0);

        fd << Label("polygons") << "[\n";
        for (int p=0; p<polygons_.size(); p++)
            DumpVector(fd, polygons_[p], 1, p==polygons_.size()-1 );
        fd << " ],\n";

        fd << Label("verts");
        DumpVector(fd, store_.verts, 1);

        fd << Label("triangles");
        auto& triangles = store_.triangles;
        ctl::PointList trPoints(triangles.size()/3);
        for (int i = 0; i + 2 < triangles.size(); i += 3)
            trPoints[i/3] = ctl::Vector(triangles[i], triangles[i + 1], triangles[i + 2]);
        DumpVector(fd, trPoints, 1, true);
        fd << "\n}";
    }

private:
    void DumpVector(std::ofstream& fd, const ctl::PointList& points, int level, bool lastInArray=false)
    {
        std::string starter;
        for (int i = 0; i < level; i++)
            starter += "    ";
        fd << starter << "[ ";
        for (int i = 0; i < points.size(); i++)
        {
            if (i > 0)
                fd << ",  ";
            // Break long vectors into lines
            if ((i & 0xF) == 0xF)
                fd << "\n" << starter;
            fd << points[i].x << "," << points[i].y << "," << points[i].z;
        }
        if (lastInArray)
            fd << " ]";
        else
            fd << " ],\n";
    }
};


static std::ofstream OpenDumpFile(wchar_t const* fileName)
{
    std::ofstream fd;
    // The file is binary, we don't want to deal with EOL conversions.
    auto flags = std::fstream::out | std::fstream::trunc | std::fstream::binary;
#ifdef WIN32_
    // Windows supports native UTF-16 strings, Linux supports only UTF-8
    fd.open(fileName, flags);
#else
    // UTF-16 to UTF-8, as per C++ 11
    std::string u8FileName = std::wstring_convert<
        std::codecvt_utf8_utf16<wchar_t>, wchar_t>{}.to_bytes(fileName);

    std::cout << "Converted file name: '" << u8FileName << "'" << std::endl;
    fd.open(u8FileName, flags);
#endif
    return fd;
}


//! \brief Gets the object pointer from the handle parameter
static DelaunayWrapper* FromHandle(uintptr_t handle)
{
    return reinterpret_cast<DelaunayWrapper*>(handle);
}


extern "C"
{
WRAPPER_API uintptr_t NewDelaunayTriangulation(double const* boundary, int beginPoint, int endPoint, int resizeIncrement, double epsilon, double areaEpsilon, int maxEdgeFlips, int settings)
{
    if (!boundary)
        return 0;
    ctl::PointList boundaryList = BuildPointList(boundary, beginPoint, endPoint);
    if (!boundaryList.size())
        return 0;  // The triangulator does not like empty boundaries
    if (resizeIncrement <= 0)
        resizeIncrement = 100000;
    if (epsilon <= 0)
        epsilon = 1e-5;
    if (areaEpsilon <= 0)
        areaEpsilon = 3e-5;
    if (maxEdgeFlips <= 0)
        maxEdgeFlips = 10000;
    if (settings < 0)
        settings = ctl::DelaunayTriangulation::CLIPPING;

    auto object = new DelaunayWrapper(boundaryList, resizeIncrement, epsilon, areaEpsilon, maxEdgeFlips, settings);
    return reinterpret_cast<uintptr_t>(object);
}

WRAPPER_API void FreeDelaunayTriangulation(uintptr_t handle)
{
    auto obj = FromHandle(handle);
    if (obj)
        delete obj;
}

WRAPPER_API int DumpDelaunayTriangulation(uintptr_t handle, const wchar_t *fileName)
{
    auto obj = FromHandle(handle);
    if (!obj || !fileName)
        return DelaunayWrapper::INVALID_PARAMETER;
    try
    {
        std::ofstream fd = OpenDumpFile(fileName);
        obj->DumpToFile(fd);
        fd.close();
        return 0;
    }
    catch (std::exception e)
    {
        // Could not write to the file
        return DelaunayWrapper::INVALID_PARAMETER;
    }
}

WRAPPER_API int ErrorDelaunayTriangulation(uintptr_t handle)
{
    auto obj = FromHandle(handle);
    return obj ? obj->tr_.error() : DelaunayWrapper::INVALID_PARAMETER;
}

WRAPPER_API int CallGatherTriangles(uintptr_t handle, const double* polygonLimits, int beginPoint, int endPoint)
{
    auto obj = FromHandle(handle);
    if (!obj)
        return DelaunayWrapper::INVALID_PARAMETER;

    std::vector<ctl::Edge*> edges;
    if (obj->RemovingExterior() && obj->polygons_.size()>0)
        edges = obj->tr_.GatherTriangles(obj->polygons_[0]);
    else
    {
        ctl::PointList polygon = BuildPointList(polygonLimits, beginPoint, endPoint);
        edges = obj->tr_.GatherTriangles(polygon);
    }
    obj->store_.CreateFromDT(&obj->tr_, edges);

    if (!obj->RemovingHoles())
        return obj->tr_.error();

    // Remove the triangles in the holes
    auto& triangles = obj->store_.triangles;
    auto& verts = obj->store_.verts;
    int numValid = 0;
    for (size_t i = 0; i < triangles.size(); i += 3)
    {
        ctl::Point centroid = verts[triangles[i + 0]] + verts[triangles[i + 1]] + verts[triangles[i + 2]];

        centroid.x /= 3;
        centroid.y /= 3;
        centroid.z /= 3;

        if (!obj->PointInHole(centroid))
        {
            // Copy the new triangle to the location of the last discarded triangle
            if (numValid < i)
            {
                triangles[numValid + 0] = triangles[i + 0];
                triangles[numValid + 1] = triangles[i + 1];
                triangles[numValid + 2] = triangles[i + 2];
            }
            numValid += 3;
        }
    }
    triangles.resize(numValid);
    return obj->tr_.error();
}

WRAPPER_API int CallGetTriangles(uintptr_t handle, int* outList, int capacity)
{
    auto obj = FromHandle(handle);
    if (!obj)
        return 0;

    auto& triangles = obj->store_.triangles;
    if (!outList || !capacity)
        return (int)(triangles.size() / 3);

    // The number of values to return is 3x the number of triangles
    capacity = std::min(capacity * 3, (int)triangles.size());
    std::copy(triangles.begin(), triangles.begin() + capacity, outList);
    return capacity / 3;
}

WRAPPER_API int CallGetVertices(uintptr_t handle, double* outList, int capacity)
{
    auto obj = FromHandle(handle);
    if (!obj)
        return 0;

    auto& verts = obj->store_.verts;
    if (!outList || !capacity)
        return (int)verts.size();

    capacity = std::min(capacity, (int)verts.size());
    for (int i = 0; i < capacity; i++, outList += 3)
    {
        outList[0] = verts[i].x;
        outList[1] = verts[i].y;
        outList[2] = verts[i].z;
    }
    return capacity;
}

WRAPPER_API int CallInsertConstrainedLineString(uintptr_t handle, double const* constraint, int beginPoint, int endPoint)
{
    auto obj = FromHandle(handle);
    if (!obj)
        return DelaunayWrapper::INVALID_PARAMETER;
    obj->tr_.InsertConstrainedLineString(obj->AddPolygon(constraint, beginPoint, endPoint));
    return obj->tr_.error();
}

WRAPPER_API int CallInsertConstrainedPolygon(uintptr_t handle, double const* constraint, int beginPoint, int endPoint)
{
    auto obj = FromHandle(handle);
    if (!obj)
        return DelaunayWrapper::INVALID_PARAMETER;
    obj->tr_.InsertConstrainedPolygon(obj->AddPolygon(constraint, beginPoint, endPoint));
    return obj->tr_.error();
}

WRAPPER_API int CallInsertPoints(uintptr_t handle, int constrained, double const* points, int beginPoint, int endPoint)
{
    auto obj = FromHandle(handle);
    if (!obj)
        return DelaunayWrapper::INVALID_PARAMETER;
    ctl::PointList pointList = BuildPointList(points, beginPoint, endPoint);
    for (ctl::Point& p : pointList)
    {
        if (constrained)
            obj->tr_.InsertConstrainedPoint(p);
        else
            obj->tr_.InsertWorkingPoint(p);

    }
    return obj->tr_.error();
}

};// extern "C"