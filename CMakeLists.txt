project(j3AstroIO)

cmake_minimum_required(VERSION 2.8.12)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

link_directories(/usr/local/lib/)

FIND_PACKAGE( OpenCV REQUIRED )
find_package (Threads REQUIRED) #threading
find_package(exiv2 REQUIRED CONFIG NAMES exiv2) 

set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra -Wshadow -Wnon-virtual-dtor -Wunused -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wno-unused")

add_compile_options(-std=c++11)
 
add_library( j3AstroIO  j3AstroIO.cpp )
set_property(TARGET j3AstroIO PROPERTY CXX_STANDARD 11)
set_property(TARGET j3AstroIO PROPERTY PUBLIC_HEADER j3AstroIO.hpp)

target_include_directories(j3AstroIO  PRIVATE ${OpenCV_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR} /usr/local/include )

target_link_libraries( j3AstroIO   ${OpenCV_LIBS} cfitsio  magic exiv2)

install(TARGETS j3AstroIO  DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(TARGETS j3AstroIO  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

add_executable( j3AstroIOconvert  j3AstroIOconvert.cpp )
target_link_libraries( j3AstroIOconvert PRIVATE j3AstroIO  )

install(TARGETS j3AstroIOconvert DESTINATION bin)
