project(j3AstroIOraw)

cmake_minimum_required(VERSION 2.8.12)

# TBD
message("!!!NEED TO INCLUDE THE MAGIC FILE!!!")


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

link_directories(/usr/local/lib/)
link_directories(/usr/lib/)

  message(STATUS "THIS IS A CI RUN")
  set(CMAKE_FIND_STATIC FIRST) 
  set(BUILD_SHARED_LIBS=OFF)
  SET(CMAKE_FIND_LIBRARY_SUFFIXES .a $(CMAKE_FIND_LIRBARY_SUFFIXES))
  SET(LINK_DEPENDS_NO_SHARED=ON)

FIND_PACKAGE( OpenCV REQUIRED )
find_package (Threads REQUIRED) #threading
FIND_PACKAGE( Eigen3 3.3 NO_MODULE )


set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra -Wshadow -Wnon-virtual-dtor -Wunused -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wno-unused")

add_compile_options(-std=c++11)

add_library( j3AstroIOraw  ../j3AstroIOraw.cpp )
set_property(TARGET j3AstroIOraw PROPERTY CXX_STANDARD 11)
set_property(TARGET j3AstroIOraw PROPERTY PUBLIC_HEADER ../j3AstroIOraw.hpp)
target_include_directories(j3AstroIOraw PRIVATE ${OpenCV_INCLUDE_DIRS} ../${CMAKE_CURRENT_SOURCE_DIR} /usr/local/include )

IF( APPLE )
  message("APPLE CI")   
 # FIND_PACKAGE(ZLIB) ${ZLIB_LIBRARIES} 
  FIND_PACKAGE(JPEG)
  target_link_libraries( j3AstroIOraw PRIVATE ${OpenCV_LIBS} /usr/local/lib/libraw.a /usr/local/lib/libmagic.a bz2 /usr/local/lib/libjpeg.a)
  install(FILES /usr/local/opt/libmagic/share/misc/magic.mgc DESTINATION .)
ELSE() 
  message("LINUX CI")
  target_link_libraries( j3AstroIOraw PRIVATE ${OpenCV_LIBS} /usr/local/lib/libraw.a magic )
ENDIF()

install(TARGETS j3AstroIOraw
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  PUBLIC_HEADER DESTINATION include)

add_executable( j3AstroIOconvert  ../j3AstroIOconvert.cpp )
target_include_directories(j3AstroIOconvert PRIVATE ${OpenCV_INCLUDE_DIRS} ../${CMAKE_CURRENT_SOURCE_DIR} /usr/local/include )

if( APPLE )
  target_link_libraries( j3AstroIOconvert PRIVATE  ${OpenCV_LIBS} j3AstroIOraw ) #${ZLIB_LIBRARIES}
ELSE()
  target_link_libraries( j3AstroIOconvert PRIVATE ${OpenCV_LIBS} j3AstroIOraw )
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
