cmake_minimum_required(VERSION 3.16)

project(EllipticCurve LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)

find_package(PkgConfig REQUIRED)
pkg_check_modules(GMP REQUIRED IMPORTED_TARGET gmpxx)

add_executable(EllipticCurve main.cpp curve.cpp curve.h)

target_link_libraries(EllipticCurve PkgConfig::GMP)
