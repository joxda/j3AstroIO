project(j3AstroIOfits)

cmake_minimum_required(VERSION 2.8.12)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

link_directories(/usr/local/lib/)

FIND_PACKAGE( OpenCV REQUIRED )
find_package (Threads REQUIRED) #threading

set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra -Wshadow -Wnon-virtual-dtor -Wunused -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wno-unused")

add_compile_options(-std=c++11)

add_library( j3AstroIOfits  j3AstroIOfits.cpp )
set_property(TARGET j3AstroIOfits  PROPERTY CXX_STANDARD 11)
set_property(TARGET j3AstroIOfits  PROPERTY PUBLIC_HEADER j3AstroIOfits.hpp)
target_include_directories(j3AstroIOfits  PRIVATE ${OpenCV_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR} /usr/local/include )

target_link_libraries( j3AstroIOfits   ${OpenCV_LIBS} cfitsio  magic)

install(TARGETS j3AstroIOfits  DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(TARGETS j3AstroIOfits  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

add_executable( j3AstroIOconvert  j3AstroIOconvert.cpp )
target_link_libraries( j3AstroIOconvert PRIVATE j3AstroIOfits  )

install(TARGETS j3AstroIOconvert DESTINATION bin)

