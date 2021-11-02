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

// 3x3 determinant
template <typename T_Scalar>
T_Scalar determinant(
    T_Scalar a, T_Scalar b, T_Scalar c,
    T_Scalar d, T_Scalar e, T_Scalar f,
    T_Scalar g, T_Scalar h, T_Scalar i)
{
    return a*e*i + b*f*g + c*d*h - c*e*g - b*d*i - a*f*h;
}

// macros required for backend specification
#define DELAUNAY_VEC Vec2<T_Scalar>
#define DELAUNAY_VEC_ACCESS(V,D) V.data[D]
#define DELAUNAY_DETERMINANT determinant

#endif // ifndef DELAUNAY_BACKEND


#include "Delaunay.hpp"


#define RND ((rnd()%10000001)*0.0000001)


#if DELAUNAY_BACKEND == DELAUNAY_BACKEND_EIGEN

using Vec2d = Eigen::Matrix<double, 2, 1>;
template <typename T>
using Vector = std::vector<T, Eigen::aligned_allocator<T>>;

#else // custom backend

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
void visualize(const Vector<T>& points, std::vector<int32_t>& triangulation,
    const std::string& windowName, const cv::Scalar& lineColor, int wait=0)
{
    cv::Mat img(1024, 1024, CV_8UC3);

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


int main()
{
    Vector<Vec2d> points;
    createPoints(points, Vec2d(0.0f, 0.0f), Vec2d(1024.0f, 1024.0f));

    auto triangulation = delaunayTriangulate(points);

    visualize(points, triangulation, "delaunay demo", cv::Scalar(120, 120, 0));

    return 0;
}
