//
// Project: Delaunay
// File: main.cpp
//
// Copyright (c) 2021 Miika 'Lehdari' Lehtimäki
// You may use, distribute and modify this code under the terms
// of the licence specified in file LICENSE which is distributed
// with this source code package.
//

#include <random>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <Eigen/Dense>


// 3x3 determinant
template <typename T_Scalar>
T_Scalar inline __attribute__((always_inline)) determinant(
    const T_Scalar& a, const T_Scalar& b, const T_Scalar& c,
    const T_Scalar& d, const T_Scalar& e, const T_Scalar& f,
    const T_Scalar& g, const T_Scalar& h, const T_Scalar& i)
{
    return a*e*i + b*f*g + c*d*h - c*e*g - b*d*i - a*f*h;
}

#ifndef DELAUNAY_BACKEND
// example implementation of a custom backend

// minimal 2D vector struct
template <typename T_Scalar>
struct Vec2 {
    T_Scalar  data[2];

    // access with functor operator required only for demo dataset construction
    // so it is identical with eigen
    T_Scalar& operator()(int d) { return data[d]; }
    const T_Scalar& operator()(int d) const { return data[d]; }

    Vec2(T_Scalar x, T_Scalar y) :
        data    {x,y}
    {}
};

// macros required for backend specification
#define DELAUNAY_VEC Vec2<T_Scalar>
#define DELAUNAY_VEC_ACCESS(V,D) V.data[D]
#define DELAUNAY_DETERMINANT determinant

#endif // ifndef DELAUNAY_BACKEND


#include "Delaunay.hpp"


#define RND ((rnd()%10000001)*0.0000001)


#if DELAUNAY_BACKEND == DELAUNAY_BACKEND_EIGEN

using Vec2d = Eigen::Matrix<double, 2, 1>;
using Vec2f = Eigen::Matrix<float, 2, 1>;
template <typename T>
using Vector = std::vector<T, Eigen::aligned_allocator<T>>;

#else // custom backend

using Vec2f = Vec2<float>;
using Vec2d = Vec2<double>;
template <typename T>
using Vector = std::vector<T>;

#endif


// generate some non-uniform data (not very fancy but does the job)
void createPoints(Vector<Vec2d>& v, const Vec2d& min, const Vec2d& max, int depth = 0)
{
    static std::default_random_engine rnd(7155178);

    double terminationProbability = 0.02*depth*depth;
    if (RND < terminationProbability) {
        v.emplace_back(min(0) + RND*(max(0)-min(0)), min(1) + RND*(max(1)-min(1)));
        v.emplace_back(min(0) + RND*(max(0)-min(0)), min(1) + RND*(max(1)-min(1)));
    }
    else {
        Vec2d half(min(0)+(0.3+0.4*RND)*(max(0)-min(0)), min(1)+(0.3+0.4*RND)*(max(1)-min(1)));
        createPoints(v, min, half, depth+1);
        createPoints(v, Vec2d(half(0), min(1)), Vec2d(max(0), half(1)), depth+1);
        createPoints(v, Vec2d(min(0), half(1)), Vec2d(half(0), max(1)), depth+1);
        createPoints(v, half, max, depth+1);
    }
}


template <typename T>
void createPointsUniform(Vector<T>& v, const T& min, const T& max, int n)
{
    std::default_random_engine rnd(1507715517);
    v.reserve(n);

    for (int i=0; i<n; ++i) {
        v.emplace_back(min(0) + RND*(max(0)-min(0)), min(1) + RND*(max(1)-min(1)));
    }
}


template <typename T>
void visualize(const Vector<T>& points, std::vector<int32_t>& triangulation,
    const std::string& windowName, const cv::Scalar& lineColor, int wait=0)
{
    cv::Mat img(1024, 1024, CV_8UC3, cv::Scalar(0, 0, 0));

    for (int i=0; i<triangulation.size(); i+=3) {
        cv::line(img,
            cv::Point(points[triangulation[i]](0), 1024-(int)points[triangulation[i]](1)),
            cv::Point(points[triangulation[i+1]](0), 1024-(int)points[triangulation[i+1]](1)),
            lineColor);
        if (triangulation[i+2] != -1) {
            cv::line(img,
                cv::Point(points[triangulation[i+1]](0), 1024-(int)points[triangulation[i+1]](1)),
                cv::Point(points[triangulation[i+2]](0), 1024-(int)points[triangulation[i+2]](1)),
                lineColor);
            cv::line(img,
                cv::Point(points[triangulation[i+2]](0), 1024-(int)points[triangulation[i+2]](1)),
                cv::Point(points[triangulation[i]](0), 1024-(int)points[triangulation[i]](1)),
                lineColor);
        }
    }

    for (auto& p : points) {
        img.at<cv::Vec3b>(1023-(int)p(1), (int)p(0)) = cv::Vec3b(255, 255, 255);
    }

    cv::imshow(windowName, img);
    cv::waitKey(wait);
}


template <typename T>
inline __attribute__((always_inline)) bool inCircle(
    const T& a, const T& b, const T& c, const T& d)
{
    T dd(d(0)*d(0), d(1)*d(1));
    return determinant(
        a(0)-d(0), a(1)-d(1), (a(0)*a(0)-dd(0))+(a(1)*a(1)-dd(1)),
        b(0)-d(0), b(1)-d(1), (b(0)*b(0)-dd(0))+(b(1)*b(1)-dd(1)),
        c(0)-d(0), c(1)-d(1), (c(0)*c(0)-dd(0))+(c(1)*c(1)-dd(1))) > 0.0;
}


template <typename T>
bool checkTriangulation(const Vector<T>& points, std::vector<int32_t>& triangulation)
{
    for (size_t i=0; i<triangulation.size(); i+=3) {
        int p1 = triangulation[i];
        int p2 = triangulation[i+1];
        int p3 = triangulation[i+2];
        if (p3 == -1)
            continue;

        for (size_t j=0; j<points.size(); ++j) {
            if (j == p1 || j == p2 || j == p3)
                continue;
            if (inCircle(points[p1], points[p2], points[p3], points[j]))
                return false;
        }
    }

    return true;
}


void benchmark(void) {
    int nPoints = 10;
    bool breakLoop = false;
    while (true) {
        Vector<Vec2f> pointsFloat;
        Vector<Vec2d> pointsDouble;
        pointsFloat.clear();
        pointsDouble.clear();
        createPointsUniform(pointsFloat, Vec2f(0.0, 0.0), Vec2f(1024.0, 1024.0), nPoints);
        createPointsUniform(pointsDouble, Vec2d(0.0, 0.0), Vec2d(1024.0, 1024.0), nPoints);

        auto t1 = std::chrono::high_resolution_clock::now();
        auto triangulationFloat = delaunayTriangulate(pointsFloat);
        auto t2 = std::chrono::high_resolution_clock::now();
        auto triangulationDouble = delaunayTriangulate(pointsDouble);
        auto t3 = std::chrono::high_resolution_clock::now();

        printf("nPoints: %d, tFloat: %0.5f, tDouble: %0.5f\n", nPoints,
            std::chrono::duration<double, std::milli>(t2-t1).count(),
            std::chrono::duration<double, std::milli>(t3-t2).count());

        if (!checkTriangulation(pointsFloat, triangulationFloat)) {
            printf("Triangulation is non-delaunay!\n");
            breakLoop = true;
        }

        visualize(pointsDouble, triangulationDouble, "triangulationDouble", cv::Scalar(120, 120, 0), 20);
        visualize(pointsFloat, triangulationFloat, "triangulationFloat", cv::Scalar(0, 80, 160), 20);

        if (breakLoop) {
            cv::waitKey(0);
            break;
        }

        nPoints *= 1.1;
    }
}


int main()
{
#if 0
    Vector<Vec2d> points;
    createPoints(points, Vec2d(0.0f, 0.0f), Vec2d(1024.0f, 1024.0f));

    auto triangulation = delaunayTriangulate(points);

    visualize(points, triangulation, "delaunay demo", cv::Scalar(120, 120, 0));
#else
    benchmark();
#endif

    return 0;
}
