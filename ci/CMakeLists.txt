project(j3AstroIOexiv)

cmake_minimum_required(VERSION 2.8.12)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

link_directories(/usr/local/lib/)

  message(STATUS "THIS IS A CI RUN")
  set(CMAKE_FIND_STATIC FIRST) 
  set(BUILD_SHARED_LIBS=OFF)
  SET(CMAKE_FIND_LIBRARY_SUFFIXES .a $(CMAKE_FIND_LIRBARY_SUFFIXES))

find_package (Threads REQUIRED) #threading
find_package(exiv2 REQUIRED CONFIG NAMES exiv2) 


set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra -Wshadow -Wnon-virtual-dtor -Wunused -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wno-unused")

add_compile_options(-std=c++11)
 

add_library( j3AstroIOexiv  ../j3AstroIOexiv.cpp )
set_property(TARGET j3AstroIOexiv PROPERTY CXX_STANDARD 11)
set_property(TARGET j3AstroIOexiv PROPERTY PUBLIC_HEADER ../j3AstroIOexiv.hpp)
target_include_directories(j3AstroIOexiv PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} /usr/local/include )

IF( (APPLE) )
  target_link_libraries( j3AstroIOexiv  expat /usr/local/lib/libexiv2.a /usr/local/lib/libexiv2-xmp.a iconv /usr/local/lib/libpng.a z)
ELSE() 
  target_link_libraries( j3AstroIOexiv   /usr/local/lib/libexiv2.a /usr/local/lib/libexiv2-xmp.a pthread png z /usr/lib/x86_64-linux-gnu/libexpat.a )
ENDIF()

install(TARGETS j3AstroIOexiv DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(TARGETS j3AstroIOexiv PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

add_executable( j3AstroIOconvert  ../j3AstroIOconvert.cpp )
target_link_libraries( j3AstroIOconvert PRIVATE j3AstroIOexiv )
 
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