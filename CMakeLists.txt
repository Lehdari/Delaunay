cmake_minimum_required(VERSION 3.10)
project(Delaunay)


set(CMAKE_CXX_STANDARD 20)


option(BUILD_DELAUNAY_DEMO "Build the demo, requires OpenCV" ON)


find_package(Eigen3 REQUIRED)

if (BUILD_DELAUNAY_DEMO)
    find_package(OpenCV REQUIRED)
endif()


add_library(delaunay INTERFACE)

target_include_directories(delaunay
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(delaunay
    INTERFACE
        Eigen3::Eigen
)


if (BUILD_DELAUNAY_DEMO)

add_executable(delaunay_demo src/main.cpp)

target_link_libraries(delaunay_demo
    PUBLIC
        delaunay
        ${OpenCV_LIBS}
)

endif()
