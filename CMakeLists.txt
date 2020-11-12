project(j3AstroIO)

cmake_minimum_required(VERSION 2.8.12)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set (OpenCV_DIR "/usr/local/Cellar/opencv/*/lib/cmake/opencv4/")

FIND_PACKAGE( OpenCV REQUIRED )



set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra -Wshadow -Wnon-virtual-dtor -Wunused -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wno-unused")

add_compile_options(-std=c++11)

if(DEFINED $ENV{CI})
  set(BUILD_SHARED_LIBS=OFF)
endif()

add_library( j3AstroIO  j3AstroIO.cpp )
set_property(TARGET j3AstroIO PROPERTY CXX_STANDARD 11)
target_include_directories(j3AstroIO PRIVATE ${OpenCV_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR} /usr/local/include )

#set(MY_OpenCV_LIBS "opencv_core;opencv_highgui;")

target_link_libraries( j3AstroIO  ${OpenCV_LIBS} raw cfitsio)

install(TARGETS j3AstroIO DESTINATION lib)
#install(FILE j3AstroIO.hpp DESTINATION include)


set(CPACK_GENERATOR "TGZ")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Joachim Janz") # required
set(CPACK_PACKAGE_VENDOR "joxda")
set(CPACK_PACKAGE_VERSION_MAJOR 1)
set(CPACK_PACKAGE_VERSION_MINOR 0)
set(CPACK_PACKAGE_VERSION_PATCH 0)
set(CPACK_PACKAGE_DESCRIPTION_FILE ${CMAKE_SOURCE_DIR}/README.md)
set(CPACK_RESOURCE_FILE_README ${CMAKE_SOURCE_DIR}/README.md)
set(CPACK_RESOURCE_FILE_LICENSE ${CMAKE_SOURCE_DIR}/LICENSE)

include(CPack)