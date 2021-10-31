//
// Project: Delaunay
// File: Delaunay.hpp
//
// Copyright (c) 2021 Miika 'Lehdari' Lehtimäki
// You may use, distribute and modify this code under the terms
// of the licence specified in file LICENSE which is distributed
// with this source code package.
//

#ifndef DELAUNAY_DELAUNAY_HPP
#define DELAUNAY_DELAUNAY_HPP


#include <Eigen/Dense>
#include <vector>
#include <cstdint>


#ifdef DELAUNAY_INLINE
    #error "DELAUNAY_INLINE already defined"
#endif
#define DELAUNAY_INLINE inline __attribute__((always_inline))


namespace {

    struct Triangle {
        int neighbours[3]; // indices of neighbouring triangles
        int vertices[3]; // indices of corner vertices

        Triangle(int t1, int t2, int t3, int v1, int v2, int v3) :
            neighbours  {t1, t2, t3},
            vertices    {v1, v2, v3}
        {}
    };

    template <typename T_Scalar>
    DELAUNAY_INLINE bool ccw(
        const Eigen::Matrix<T_Scalar, 2, 1>& a,
        const Eigen::Matrix<T_Scalar, 2, 1>& b,
        const Eigen::Matrix<T_Scalar, 2, 1>& c)
    {
        Eigen::Matrix<T_Scalar, 3, 3>   m;
        m <<
            a.transpose(),  1.0,
            b.transpose(),  1.0,
            c.transpose(),  1.0;
        return m.determinant() > 0.0;
    }

    // is d inside of circumcircle of abc?
    template <typename T_Scalar>
    DELAUNAY_INLINE bool inCircle(
        const Eigen::Matrix<T_Scalar, 2, 1>& a,
        const Eigen::Matrix<T_Scalar, 2, 1>& b,
        const Eigen::Matrix<T_Scalar, 2, 1>& c,
        const Eigen::Matrix<T_Scalar, 2, 1>& d)
    {
        Eigen::Matrix<T_Scalar, 4, 4>   m;
        m <<
            a.transpose(),  a.dot(a),   1.0,
            b.transpose(),  b.dot(b),   1.0,
            c.transpose(),  c.dot(c),   1.0,
            d.transpose(),  d.dot(d),   1.0;
        return m.determinant() > 0.0;
    }

    // initial primitive construction functions for edges and triangles
    template <template <typename, typename> class T_Vector, typename T_Scalar, typename T_Allocator>
    DELAUNAY_INLINE void createEdge(const T_Vector<Eigen::Matrix<T_Scalar, 2, 1>, T_Allocator>& points,
        std::vector<Triangle>& triangles, int v1, int v2, int& firstTriangle, int& lastTriangle)
    {
        int s = triangles.size();
        triangles.emplace_back(s+1, s+1, s+1, v1, v2, -1); // edge is presented by 2 ghost triangles (3rd vertex -1)
        triangles.emplace_back(s, s, s, v2, v1, -1);

        firstTriangle = s;
        lastTriangle = s;
    }

    template <template <typename, typename> class T_Vector, typename T_Scalar, typename T_Allocator>
    DELAUNAY_INLINE void createTriangle(
        T_Vector<Eigen::Matrix<T_Scalar, 2, 1>, T_Allocator>& points,
        std::vector<Triangle>& triangles, int v1, int v2, int v3, int& firstTriangle, int& lastTriangle)
    {
        int s = triangles.size();
        if (!ccw(points[v1], points[v2], points[v3])) {
            std::swap(v2, v3);
            lastTriangle = s+2;
        }
        else
            lastTriangle = s+3;
        triangles.emplace_back(s+1, s+2, s+3, v1, v2, v3);
        triangles.emplace_back(s, s+3, s+2, v2, v1, -1);
        triangles.emplace_back(s, s+1, s+3, v3, v2, -1);
        triangles.emplace_back(s, s+2, s+1, v1, v3, -1);

        firstTriangle = s+3;
    }

    // update neighbour of a triangle
    // t: triangle neighbour of which is to be updated
    // current: current neighbour triangle id to be changed
    // updated: triangle id of the new neighbour
    DELAUNAY_INLINE void updateNeighbour(Triangle& t, int current, int updated)
    {
        if (t.neighbours[0] == current) {
            t.neighbours[0] = updated;
            return;
        }
        if (t.neighbours[1] == current) {
            t.neighbours[1] = updated;
            return;
        }
        if (t.neighbours[2] == current) {
            t.neighbours[2] = updated;
            return;
        }
    }

    // rotate triangle indices so that vertex id of -1 is at index 2, no-op for non-ghost or correct ghost triangles
    DELAUNAY_INLINE void correctGhost(Triangle& t)
    {
        if (t.vertices[0] == -1) {
            int tempVertex = t.vertices[2];
            int tempNeighbour = t.neighbours[2];
            t.vertices[2] = t.vertices[0];
            t.neighbours[2] = t.neighbours[0];
            t.vertices[0] = t.vertices[1];
            t.neighbours[0] = t.neighbours[1];
            t.vertices[1] = tempVertex;
            t.neighbours[1] = tempNeighbour;
        }
        else if (t.vertices[1] == -1) {
            int tempVertex = t.vertices[2];
            int tempNeighbour = t.neighbours[2];
            t.vertices[2] = t.vertices[1];
            t.neighbours[2] = t.neighbours[1];
            t.vertices[1] = t.vertices[0];
            t.neighbours[1] = t.neighbours[0];
            t.vertices[0] = tempVertex;
            t.neighbours[0] = tempNeighbour;
        }
    }

    // flip an edge between two triangles
    inline void flip(std::vector<Triangle>& triangles, int t1, int t2)
    {
        Triangle& triangle1 = triangles[t1];
        Triangle& triangle2 = triangles[t2];

        int t1EdgeId = 0; // id of the t1 edge connecting to t2
        if (triangle1.neighbours[1] == t2)
            t1EdgeId = 1;
        else if (triangle1.neighbours[2] == t2)
            t1EdgeId = 2;

        int t2EdgeId = 0; // id of the t2 edge connecting to t1
        if (triangle2.neighbours[1] == t1)
            t2EdgeId = 1;
        else if (triangle2.neighbours[2] == t1)
            t2EdgeId = 2;

        Triangle triangle1temp = triangle1;
        Triangle triangle2temp = triangle2;

        // update neighbour indexing
        updateNeighbour(triangles[triangle1temp.neighbours[(t1EdgeId+1)%3]], t1, t2);
        updateNeighbour(triangles[triangle2temp.neighbours[(t2EdgeId+1)%3]], t2, t1);

        // update vertex and neighbour indexing of t1 and t2
        triangle1.vertices[(t1EdgeId+1)%3] = triangle2temp.vertices[(t2EdgeId+2)%3];
        triangle1.neighbours[t1EdgeId] = triangle2temp.neighbours[(t2EdgeId+1)%3];
        triangle1.neighbours[(t1EdgeId+1)%3] = t2;

        triangle2.vertices[(t2EdgeId+1)%3] = triangle1temp.vertices[(t1EdgeId+2)%3];
        triangle2.neighbours[t2EdgeId] = triangle1temp.neighbours[(t1EdgeId+1)%3];
        triangle2.neighbours[(t2EdgeId+1)%3] = t1;

        // correct ghost indexing in case either of the triangles is a ghost
        correctGhost(triangle1);
        correctGhost(triangle2);
    }

    // get vertex of triangle t opposing neighbour neighbour
    DELAUNAY_INLINE int getOpposingVertex(Triangle& t, int neighbour)
    {
        if (t.neighbours[0] == neighbour)
            return t.vertices[2];
        if (t.neighbours[1] == neighbour)
            return t.vertices[0];
        if (t.neighbours[2] == neighbour)
            return t.vertices[1];

        // should never be reached
        return -1; // neighbour not neighbour of t
    }

    // get the triangle ID of right neighbour w.r.t to a vertex id
    DELAUNAY_INLINE int getRightNeighbour(Triangle& t, int v)
    {
        if (t.vertices[0] == v)
            return t.neighbours[2];
        if (t.vertices[1] == v)
            return t.neighbours[0];
        if (t.vertices[2] == v)
            return t.neighbours[1];

        // should never be reached
        return -1; // v not a vertex of t
    }

    // get the triangle ID of left neighbour w.r.t to a vertex id
    DELAUNAY_INLINE int getLeftNeighbour(Triangle& t, int v)
    {
        if (t.vertices[0] == v)
            return t.neighbours[0];
        if (t.vertices[1] == v)
            return t.neighbours[1];
        if (t.vertices[2] == v)
            return t.neighbours[2];

        // should never be reached
        return -1; // v not a vertex of t
    }

    // ghost triangle adding functions required in ends of a "seam"
    DELAUNAY_INLINE int addGhostTriangle(std::vector<Triangle>& triangles, int v1, int v2, int n2, int n3)
    {
        int s = triangles.size();
        triangles.emplace_back(-1, n2, n3, v1, v2, -1);
        triangles[n2].neighbours[2] = s;
        triangles[n3].neighbours[1] = s;
        return s;
    }

    DELAUNAY_INLINE int addGhostTriangle(std::vector<Triangle>& triangles, int v1, int v2, int n1, int n2, int n3)
    {
        int s = triangles.size();
        triangles.emplace_back(n1, n2, n3, v1, v2, -1);
        triangles[n2].neighbours[2] = s;
        triangles[n3].neighbours[1] = s;
        return s;
    }

    template <template <typename, typename> class T_Vector, typename T_Scalar, typename T_Allocator>
    void merge(T_Vector<Eigen::Matrix<T_Scalar, 2, 1>, T_Allocator>& points,
        std::vector<Triangle>& triangles, int& firstLeft, int lastLeft, int firstRight, int& lastRight)
    {
        // keep track of end vertices in case the end triangles change
        int firstVertex = triangles[firstLeft].vertices[0];
        int lastVertex = triangles[lastRight].vertices[1];

        // find the lower common tangent
        while (true) {
            int nOps = 0;
            int leftNext = triangles[lastLeft].neighbours[1]; // next ghost triangle on the left mesh boundary
            int rightNext = triangles[firstRight].neighbours[2]; // next ghost triangle on the right mesh boundary
            if (ccw(
                points[triangles[leftNext].vertices[0]],
                points[triangles[leftNext].vertices[1]],
                points[triangles[firstRight].vertices[0]])) {
                lastLeft = leftNext;
                ++nOps;
            }
            if (ccw(
                points[triangles[rightNext].vertices[0]],
                points[triangles[rightNext].vertices[1]],
                points[triangles[lastLeft].vertices[1]])) {
                firstRight = rightNext;
                ++nOps;
            }
            if (nOps == 0) {
                // opposing triangles are right of the edge of both ghost triangles,
                // lower common tangent has been found
                break;
            }
        }

        bool leftValid = true;
        bool rightValid = true;
        int baseLeftVertex = triangles[lastLeft].vertices[1];
        int baseRightVertex = triangles[firstRight].vertices[0];

        // add new ghost triangle to the beginning of the seam
        int firstGhost = addGhostTriangle(triangles, baseRightVertex, baseLeftVertex,
            triangles[lastLeft].neighbours[1], triangles[firstRight].neighbours[2]);
        // check if the ghost overrides either of the end triangles - in that case update their indices
        if (triangles[firstGhost].vertices[0] == firstVertex)
            firstLeft = firstGhost;
        if (triangles[firstGhost].vertices[1] == lastVertex)
            lastRight = firstGhost;

        int lastConnectedSide = -1; // required for adding the end ghost, 0: left, 1: right
        int lastConnectedTriangle = -1;

        // merge loop
        while (leftValid || rightValid) {
            auto& baseLeft = points[baseLeftVertex];
            auto& baseRight = points[baseRightVertex];

            // find connective candidates for both sides, delete non-delaunay edges by flipping
            int leftCand = triangles[lastLeft].vertices[0];
            leftValid = ccw(baseLeft, baseRight, points[leftCand]);
            if (leftValid) {
                int leftNeighbour = getRightNeighbour(triangles[lastLeft], baseLeftVertex);
                int leftNeighbourCand = getOpposingVertex(triangles[leftNeighbour], lastLeft);
                while (leftNeighbourCand != -1 &&
                    inCircle(baseLeft, baseRight, points[leftCand], points[leftNeighbourCand])) {
                    flip(triangles, lastLeft, leftNeighbour);
                    lastLeft = leftNeighbour;
                    leftCand = leftNeighbourCand;
                    leftNeighbour = getRightNeighbour(triangles[lastLeft], baseLeftVertex);
                    leftNeighbourCand = getOpposingVertex(triangles[leftNeighbour], lastLeft);
                }
            }

            int rightCand = triangles[firstRight].vertices[1];
            rightValid = ccw(baseLeft, baseRight, points[rightCand]);
            if (rightValid) {
                int rightNeighbour = getLeftNeighbour(triangles[firstRight], baseRightVertex);
                int rightNeighbourCand = getOpposingVertex(triangles[rightNeighbour], firstRight);
                while (rightNeighbourCand != -1 &&
                    inCircle(baseLeft, baseRight, points[rightCand], points[rightNeighbourCand])) {
                    flip(triangles, firstRight, rightNeighbour);
                    rightCand = rightNeighbourCand;
                    rightNeighbour = getLeftNeighbour(triangles[firstRight], baseRightVertex);
                    rightNeighbourCand = getOpposingVertex(triangles[rightNeighbour], firstRight);
                }
            }

            if (!leftValid && !rightValid)
                break;

            if (!leftValid || (rightValid && inCircle(points[leftCand], baseLeft, baseRight, points[rightCand]))) {
                int rightNext = triangles[firstRight].neighbours[1];

                if (firstGhost != -1) { // update first ghost
                    triangles[firstGhost].neighbours[0] = firstRight;
                    triangles[firstRight].neighbours[2] = firstGhost;
                    firstGhost = -1;
                }
                else { // connect left side triangle in case it was previously added
                    if (lastConnectedSide == 0) {
                        triangles[lastConnectedTriangle].neighbours[2] = firstRight;
                        triangles[firstRight].neighbours[2] = lastConnectedTriangle;
                    }
                }

                // connect edge
                triangles[firstRight].vertices[2] = baseLeftVertex;
                lastConnectedTriangle = firstRight;
                firstRight = rightNext;
                baseRightVertex = rightCand;

                lastConnectedSide = 1;
            }
            else {
                int leftNext = triangles[lastLeft].neighbours[2];

                if (firstGhost != -1) { // update first ghost
                    triangles[firstGhost].neighbours[0] = lastLeft;
                    triangles[lastLeft].neighbours[1] = firstGhost;
                    firstGhost = -1;
                }
                else { // connect right side triangle in case it was previously added
                    if (lastConnectedSide == 1) {
                        triangles[lastConnectedTriangle].neighbours[1] = lastLeft;
                        triangles[lastLeft].neighbours[1] = lastConnectedTriangle;
                    }
                }

                // connect edge
                triangles[lastLeft].vertices[2] = baseRightVertex;
                lastConnectedTriangle = lastLeft;
                lastLeft = leftNext;
                baseLeftVertex = leftCand;

                lastConnectedSide = 0;
            }
        }

        // add new ghost triangle to the end of the seam
        int lastGhost;
        if (lastConnectedSide == 0) {
            int lastConnect = triangles[lastLeft].neighbours[1];
            lastGhost = addGhostTriangle(triangles, baseLeftVertex, baseRightVertex, lastConnect, firstRight, lastLeft);
            triangles[lastConnect].neighbours[2] = lastGhost;
        }
        else {
            int lastConnect = triangles[firstRight].neighbours[2];
            lastGhost = addGhostTriangle(triangles, baseLeftVertex, baseRightVertex, lastConnect, firstRight, lastLeft);
            triangles[lastConnect].neighbours[1] = lastGhost;
        }

        // check if the ghost overrides either of the end triangles - in that case update their indices
        if (triangles[lastGhost].vertices[0] == firstVertex)
            firstLeft = lastGhost;
        if (triangles[lastGhost].vertices[1] == lastVertex)
            lastRight = lastGhost;
    }

    template <template <typename, typename> class T_Vector, typename T_Scalar, typename T_Allocator>
    void construct(T_Vector<Eigen::Matrix<T_Scalar, 2, 1>, T_Allocator>& points,
        int begin, int end, std::vector<Triangle>& triangles, int& firstTriangle, int& lastTriangle)
    {
        int d = end-begin;
        if (d == 2) { // initial primitives, edge
            createEdge(points, triangles, begin, begin+1, firstTriangle, lastTriangle);
        }
        else if (d == 3) { // initial primitives, triangle
            createTriangle(points, triangles, begin, begin+1, begin+2, firstTriangle, lastTriangle);
        }
        else {
            // recursive construction and merge
            auto half = begin+d/2;
            int firstLeft; int lastLeft;
            construct(points, begin, half, triangles, firstLeft, lastLeft);
            int firstRight; int lastRight;
            construct(points, half, end, triangles, firstRight, lastRight);

            merge(points, triangles, firstLeft, lastLeft, firstRight, lastRight);

            firstTriangle = firstLeft;
            lastTriangle = lastRight;
        }
    }

}

template <template <typename, typename> class T_Vector, typename T_Scalar, typename T_Allocator>
std::vector<int32_t> delaunayTriangulate(T_Vector<Eigen::Matrix<T_Scalar, 2, 1>, T_Allocator>& points)
{
    // sort points primarily along x-axis (and along y-axis in case of equal x)
    std::sort(points.begin(), points.end(), [](
        const Eigen::Matrix<T_Scalar, 2, 1>& a,
        const Eigen::Matrix<T_Scalar, 2, 1>& b){
        return a(0) == b(0) ? a(1) < b(1) : a(0) < b(0);
    });

    std::vector<Triangle> triangles;
    triangles.reserve(2*points.size());

    int firstTriangle, lastTriangle;
    construct(points, 0, points.size(), triangles, firstTriangle, lastTriangle);

    // list indices into output vector
    std::vector<int32_t> indices;
    indices.reserve(triangles.size()*3);
    for (auto& t : triangles) {
        indices.emplace_back(t.vertices[0]);
        indices.emplace_back(t.vertices[1]);
        indices.emplace_back(t.vertices[2]);
    }
    return indices;
}


#endif //DELAUNAY_DELAUNAY_HPP
