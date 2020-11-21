project(j3AstroIOraw)

cmake_minimum_required(VERSION 2.8.12)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

link_directories(/usr/local/lib/)
link_directories(/usr/lib/)

FIND_PACKAGE( OpenCV REQUIRED )
find_package (Threads REQUIRED) #threading
FIND_PACKAGE( Eigen3 3.3 NO_MODULE )


set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra -Wshadow -Wnon-virtual-dtor -Wunused -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wno-unused")

add_compile_options(-std=c++11)

add_library( j3AstroIOraw  j3AstroIOraw.cpp )
set_property(TARGET j3AstroIOraw PROPERTY CXX_STANDARD 11)
set_property(TARGET j3AstroIOraw PROPERTY PUBLIC_HEADER j3AstroIOraw.hpp)
target_include_directories(j3AstroIOraw PRIVATE ${OpenCV_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR} /usr/local/include )

IF( TARGET Eigen3::Eigen )
 target_link_libraries( j3AstroIOraw  ${OpenCV_LIBS} raw magic Eigen3::Eigen )
ELSE()
 target_link_libraries( j3AstroIOraw  ${OpenCV_LIBS} raw magic)
ENDIF()

install(TARGETS j3AstroIOraw
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  PUBLIC_HEADER DESTINATION include)

add_executable( j3AstroIOconvert  j3AstroIOconvert.cpp )
target_include_directories(j3AstroIOconvert PRIVATE ${OpenCV_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR} /usr/local/include )

target_link_libraries( j3AstroIOconvert PRIVATE ${OpenCV_LIBS} j3AstroIOraw )

install(TARGETS j3AstroIOconvert DESTINATION bin)

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
