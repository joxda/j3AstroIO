project(j3AstroIOraw)

cmake_minimum_required(VERSION 2.8.12)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

link_directories(/usr/local/lib/)

if(DEFINED $ENV{CI})
  message(STATUS "THIS IS A CI RUN")
  set(CMAKE_FIND_STATIC FIRST) 
  set(BUILD_SHARED_LIBS=OFF)
  SET(CMAKE_FIND_LIBRARY_SUFFIXES .a $(CMAKE_FIND_LIRBARY_SUFFIXES))
endif()

FIND_PACKAGE( OpenCV REQUIRED )
find_package (Threads REQUIRED) #threading
FIND_PACKAGE( Eigen3 3.3 NO_MODULE )

#IF( (DEFINED $ENV{CI}) AND (APPLE) )
# find_package( FindIconv ) 
#ENDIF()

set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra -Wshadow -Wnon-virtual-dtor -Wunused -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wno-unused")

add_compile_options(-std=c++11)


add_library( j3AstroIOraw  j3AstroIOraw.cpp )
set_property(TARGET j3AstroIOraw PROPERTY CXX_STANDARD 11)
set_property(TARGET j3AstroIOraw PROPERTY PUBLIC_HEADER j3AstroIO.hpp)
target_include_directories(j3AstroIOraw PRIVATE ${OpenCV_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR} /usr/local/include )

IF( (DEFINED $ENV{CI}) AND (APPLE) )
  target_link_libraries( j3AstroIOraw  ${OpenCV_LIBS} /usr/local/lib/libraw.a /usr/local/lib/libmagic.a )
ELSEIF( DEFINED $ENV{CI} )
  target_link_libraries( j3AstroIOraw  ${OpenCV_LIBS} /usr/local/lib/libraw.a magic )
ELSEIF( TARGET Eigen3::Eigen )
  target_link_libraries( j3AstroIOraw  ${OpenCV_LIBS} raw magic Eigen3::Eigen )
ELSE()
  target_link_libraries( j3AstroIOraw  ${OpenCV_LIBS} raw magic)
ENDIF()

install(TARGETS j3AstroIOraw
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  PUBLIC_HEADER DESTINATION include)

add_executable( j3AstroIOconvert  j3AstroIOconvert.cpp )
if( (DEFINED $ENV{CI}) AND (APPLE) )
  target_link_libraries( j3AstroIOconvert PRIVATE j3AstroIOraw /usr/local/lib/libiconv.dylib)
ELSE()
  target_link_libraries( j3AstroIOconvert PRIVATE j3AstroIOraw )
ENDIF()

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
