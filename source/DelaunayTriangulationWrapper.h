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
#include <cstddef>

#ifdef DELAUNAY_TRIANGULATION_EXPORTS
#ifdef WIN32
#define WRAPPER_API __declspec(dllexport)
#else
#include <locale>
#include <codecvt>
#define WRAPPER_API __attribute__((visibility("default")))
#endif
#else
#define WRAPPER_API
#endif

extern "C"
{

/*! \brief Allocates a new triangulation

    To make integration of this API simple, it uses simple arrays to transfer parameters. 

    An array of ctl::Point is received as a vector of double values, each point uses three
    values. 
    
    The points in use are selected by the [beginPoint,endPoint) indexes (endPoint is NOT included). 
    So it is possible to put all the points in a single array and select a subset of them for each call.

    See the DelaunayTriangulation constructor for other parameter details. Parameters
    can be set to -1 to select default values.

    \code
    double square[] = { -1,-1,0,    1,-1,0,  1,1,0,   -1,1,0,  -1,-1, 0 };
    double bigSquare[] = { -2,-2,0,    2,-2,0,  2,2,0,   -2,2,0,  -2,-2, 0 };
    double boundary[] = { -3,-3,0,    3,-3,0,  3,3,0,   -3,3,0,  -3,-3, 0 };
    uintptr_t tr = NewDelaunayTriangulation(boundary, 0, 4, -1, -1, -1, -1, ctl::DelaunayTriangulation::CLIPPING);

    CallInsertConstrainedLineString(tr, bigSquare, 0, 5);  // Insert the outline
    CallInsertConstrainedLineString(tr, square, 0, 5);     // Insert the holes
    CallGatherTriangles(tr, NULL, 0, 0);
    int Nt = CallGetTriangles(tr, NULL, 0);
    int* triangles = new int[Nt * 3];
    CallGetTriangles(tr, triangles, Nt);
    int Nv = CallGetVertices(tr, NULL, 0);
    double* vertices = new double[Nv * 3];
    CallGetVertices(tr, vertices, Nv);
    DumpDelaunayTriangulation(tr, "triangles.json");
    FreeDelaunayTriangulation(tr);

    \endcode 

    \param boundary Contains groups of three values (X,Y,Z) each making up a point.
    \param beginPoint,endPoint  Range of point indexes to use for the bounday.
    \return A pointer to the triangulation, NULL if the boundary is not provided.
            The pointer must be provided to all the other calls. 
            It MUST be destroyed with FreeDelaunayTriangulation
*/
WRAPPER_API uintptr_t NewDelaunayTriangulation(double const* boundary, int beginPoint, int endPoint, int resizeIncrement, double epsilon, double areaEpsilon, int maxEdgeFlips, int settings);

WRAPPER_API void FreeDelaunayTriangulation(uintptr_t handle);

//! \brief After GatherTriangles, generates a JSON file with the input and output of the triangulation
//  \param fileName is UTF-16 encoded
//! \return Returns 0 on success, or INVALID_PARAMETER if it cannot write to the file
WRAPPER_API int DumpDelaunayTriangulation(uintptr_t handle, const wchar_t* fileName);

//! \brief Returns the errors reported by the DelaunayTriangulaion, IVALID_PARAMETER if not set
WRAPPER_API int ErrorDelaunayTriangulation(uintptr_t handle);
 
/*! \brief Calculates the triangulation with the contraints added so far

    The triangles and vertices are retrieved using CallGetTriangles() and CallGetVertices()

    The function can remove some triangles from this triangulation:
        -if the REMOVE_HOLES flag is set, it will remove all triangles inside the hole constraints.
        -if the REMOVE_EXTERIOR flag is set, it will remove all triangles outside the 
         outer constraint
        -if the polygonLimits polygon is provided and the REMOVE_EXTERIOR constraint is not
         set, it will remove all triangles outside the polygonLimits.

    \param polygonLimits A vector of points, each point comprieses three double values
    \param beginPoint,endPoint Range of points within polygonLimits which should be used (endPoint is NOT inclusive)
    \param excludeHoles If not zero, the triangles that fall inside the holes are not returned
    \return Current error status
*/
WRAPPER_API int CallGatherTriangles(uintptr_t handle, const double* polygonLimits, int beginPoint, int endPoint);

/*! \brief Returns the triangles, and/or the number of triangles
    \param outList For each triangle returns three indexes into the vertices array
    \param capacity Number of triangles which can be stored in outList. Note that
             each triangle needs three numbers, so the size of the array must be
             three times the number of triangles.
             Set the capacity to zero to retrieve the total number of triangles             
    \return  The number of triangles copied to outList, or the total number of triangles
*/
WRAPPER_API int CallGetTriangles(uintptr_t handle, int* outList, int capacity);

/*! \brief Returns the vertices from the triangulation, and/or the number of vertices
    \param outList For each vertex returns the three coordinates
    \param capacity Number of vertices which can be stored in outList. Note that
             each vertex needs three numbers, so the size of the array must be
             three times the number of vertices.
             Set the capacity to zero to retrieve the total number of vertices
    \return  The number of vertices copied to outList, or the total number of vertices
*/
WRAPPER_API int CallGetVertices(uintptr_t handle, double* outList, int capacity);

/*! \brief Wrapper for DelaunayTriangulation::InsertConstrainedLineString

    \param constraint A vector of points, each point comprieses three double values
    \param beginPoint,endPoint Range of points which should be used (endPoint is NOT inclusive)
    \return Current error status
*/
WRAPPER_API int CallInsertConstrainedLineString(uintptr_t handle, double const* constraint, int beginPoint, int endPoint);

/*! \brief Wrapper for DelaunayTriangulation::InsertConstrainedPolygon

    \param constraint A vector of points, each point comprieses three double values
    \param beginPoint,endPoint Range of points which should be used (endPoint is NOT inclusive)
    \return Current error status
*/
WRAPPER_API int CallInsertConstrainedPolygon(uintptr_t handle, double const* constraint, int beginPoint, int endPoint);


/*! \brief Wrapper for DelaunayTriangulation::InsertWorkingPoint and InsertConstrainedPoint

    \param constrained If true insert as constrained points, if false insert as working points
    \param points A vector of points, each point comprises three double values
    \param beginPoint,endPoint Range of points which should be used (endPoint is NOT inclusive)
    \return Current error status
*/
WRAPPER_API int CallInsertPoints(uintptr_t handle, int constrained, double const* points, int beginPoint, int endPoint);


};// extern "C"