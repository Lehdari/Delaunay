cmake_minimum_required(VERSION 3.10)
project(Delaunay)


set(CMAKE_CXX_STANDARD 20)


option(BUILD_DELAUNAY_DEMO "Build the demo, requires OpenCV" ON)
option(DELAUNAY_USE_EIGEN "Use Eigen3 as backend" ON)


if (DELAUNAY_USE_EIGEN)
    find_package(Eigen3 REQUIRED)
endif()

if (BUILD_DELAUNAY_DEMO)
    find_package(OpenCV REQUIRED)
endif()


add_library(delaunay INTERFACE)

target_include_directories(delaunay
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

if (DELAUNAY_USE_EIGEN)
    target_link_libraries(delaunay
        INTERFACE
            Eigen3::Eigen
    )

    target_compile_definitions(delaunay
        INTERFACE
            DELAUNAY_BACKEND=DELAUNAY_BACKEND_EIGEN
    )
endif()


if (BUILD_DELAUNAY_DEMO)

add_executable(delaunay_demo src/main.cpp)

target_link_libraries(delaunay_demo
    PUBLIC
        delaunay
        ${OpenCV_LIBS}
)

endif()
