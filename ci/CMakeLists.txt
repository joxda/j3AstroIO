project(j3AstroIO)

cmake_minimum_required(VERSION 2.8.12)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

link_directories(/usr/local/lib/)

message(STATUS "THIS IS A CI RUN")
set(CMAKE_FIND_STATIC FIRST) 
set(BUILD_SHARED_LIBS=OFF)
SET(CMAKE_FIND_LIBRARY_SUFFIXES .a $(CMAKE_FIND_LIRBARY_SUFFIXES))

FIND_PACKAGE( OpenCV REQUIRED )
find_package (Threads REQUIRED) #threading
find_package(exiv2 REQUIRED CONFIG NAMES exiv2) 
FIND_PACKAGE( Eigen3 3.3 NO_MODULE )

set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra -Wshadow -Wnon-virtual-dtor -Wunused -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wno-unused")

add_compile_options(-std=c++11)
 

add_library( j3AstroIO  ../j3AstroIO.cpp )
set_property(TARGET j3AstroIO PROPERTY CXX_STANDARD 11)
set_property(TARGET j3AstroIO PROPERTY PUBLIC_HEADER ../j3AstroIO.hpp)
target_include_directories(j3AstroIO PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${OpenCV_INCLUDE_DIRS} /usr/local/include )

IF( (APPLE) )
  target_link_libraries( j3AstroIO  /usr/local/lib/libexiv2.a iconv /usr/local/lib/libpng.a z  ${OpenCV_LIBS} cfitsio  /usr/local/lib/libmagic.a)
ELSE() 
  target_link_libraries( j3AstroIO  ${OpenCV_LIBS} cfitsio  magic /usr/local/lib/libexiv2.a pthread png z)
ENDIF()

install(TARGETS j3AstroIO DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(TARGETS j3AstroIO PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

add_executable( j3AstroIOconvert  ../j3AstroIOconvert.cpp )

IF( (APPLE) )
  target_link_libraries( j3AstroIOconvert PRIVATE j3AstroIO z)
ELSE()
  target_link_libraries( j3AstroIOconvert PRIVATE j3AstroIO )
ENDIF()

install(TARGETS j3AstroIOconvert DESTINATION bin)

set(CPACK_GENERATOR "TGZ")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Joachim Janz") # required
set(CPACK_PACKAGE_VENDOR "joxda")
set(CPACK_PACKAGE_VERSION_MAJOR 1)
set(CPACK_PACKAGE_VERSION_MINOR 0)
set(CPACK_PACKAGE_VERSION_PATCH 0)
set(CPACK_PACKAGE_DESCRIPTION_FILE ${CMAKE_SOURCE_DIR}/../README.md)
set(CPACK_RESOURCE_FILE_README ${CMAKE_SOURCE_DIR}/../README.md)
set(CPACK_RESOURCE_FILE_LICENSE ${CMAKE_SOURCE_DIR}/../LICENSE)

include(CPack)