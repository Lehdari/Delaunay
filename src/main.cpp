//
// Project: Delaunay
// File: main.cpp
//
// Copyright (c) 2021 Miika 'Lehdari' Lehtimäki
// You may use, distribute and modify this code under the terms
// of the licence specified in file LICENSE which is distributed
// with this source code package.
//

#include "Delaunay.hpp"
#include <random>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>


#define RND ((rnd()%10000001)*0.0000001)


using Vec2d = Eigen::Matrix<double, 2, 1>;
template <typename T>
using Vector = std::vector<T, Eigen::aligned_allocator<T>>;


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


int main()
{
    Vector<Vec2d> points;
    createPoints(points, Vec2d(0.0f, 0.0f), Vec2d(1024.0f, 1024.0f));

    auto triangulation = delaunayTriangulate(points);

    cv::Mat img(1024, 1024, CV_8UC3);

    for (int i=0; i<triangulation.size(); i+=3) {
        if (triangulation[i+2] != -1) {
            cv::line(img,
                cv::Point(points[triangulation[i]](0), 1024-(int)points[triangulation[i]](1)),
                cv::Point(points[triangulation[i+1]](0), 1024-(int)points[triangulation[i+1]](1)),
                cv::Scalar(120, 120, 0));
            cv::line(img,
                cv::Point(points[triangulation[i+1]](0), 1024-(int)points[triangulation[i+1]](1)),
                cv::Point(points[triangulation[i+2]](0), 1024-(int)points[triangulation[i+2]](1)),
                cv::Scalar(120, 120, 0));
            cv::line(img,
                cv::Point(points[triangulation[i+2]](0), 1024-(int)points[triangulation[i+2]](1)),
                cv::Point(points[triangulation[i]](0), 1024-(int)points[triangulation[i]](1)),
                cv::Scalar(120, 120, 0));
        }
        else {
            cv::line(img,
                cv::Point(points[triangulation[i]](0), 1024-(int)points[triangulation[i]](1)),
                cv::Point(points[triangulation[i+1]](0), 1024-(int)points[triangulation[i+1]](1)),
                cv::Scalar(120, 120, 0));
        }
    }

    for (auto& p : points) {
        img.at<cv::Vec3b>(1023-(int)p(1), (int)p(0)) = cv::Vec3b(255, 255, 255);
    }

    cv::imshow("delaunay demo", img);
    cv::waitKey();

    return 0;
}
