project(j3AstroIO)

cmake_minimum_required(VERSION 3.11.99)


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

link_directories(/usr/local/lib/)

message(STATUS "THIS IS A CI RUN")
set(CMAKE_FIND_STATIC FIRST) 
set(BUILD_SHARED_LIBS=OFF)
SET(CMAKE_FIND_LIBRARY_SUFFIXES .a $(CMAKE_FIND_LIRBARY_SUFFIXES))
SET(LINK_DEPENDS_NO_SHARED=ON)

find_package (Threads REQUIRED) #threading
find_package(exiv2 REQUIRED CONFIG NAMES exiv2) 
FIND_PACKAGE( OpenCV REQUIRED )
FIND_PACKAGE( Eigen3 3.3 NO_MODULE )

link_directories(/usr/lib/)

add_compile_definitions(CIRUN)

include_directories(embed)
add_subdirectory(embed)

embed_resources( someResources magic.mgc )

set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra -Wshadow -Wnon-virtual-dtor -Wunused -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wno-unused")

add_library( j3AstroIO  ../j3AstroIO.cpp ${someResources} )
set_property(TARGET j3AstroIO PROPERTY CXX_STANDARD 11)
set_property(TARGET j3AstroIO PROPERTY PUBLIC_HEADER ../j3AstroIO.hpp)
target_include_directories(j3AstroIO PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../ /usr/local/include ${OpenCV_INCLUDE_DIRS} )
set_property(TARGET j3AstroIO PROPERTY CXX_STANDARD 11) 

IF( (APPLE) )
  message("APPLE CI")   
 # FIND_PACKAGE(ZLIB) ${ZLIB_LIBRARIES} 
  FIND_PACKAGE(JPEG)
  target_link_libraries( j3AstroIO  PRIVATE /usr/local/lib/libexiv2.a iconv /usr/local/lib/libpng.a bz2 curl cfitsio  /usr/local/lib/libmagic.a ${OpenCV_LIBS} /usr/local/lib/libraw.a /usr/local/lib/libjpeg.a )
#  install(FILES /usr/local/opt/libmagic/share/misc/magic.mgc DESTINATION .)
ELSE() 
  message("LINUX CI")
  target_link_libraries( j3AstroIO PRIVATE ${OpenCV_LIBS} /usr/local/lib/libraw.a cfitsio  magic /usr/local/lib/libexiv2.a pthread )
ENDIF()

install(TARGETS j3AstroIO DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(TARGETS j3AstroIO PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

add_executable( j3AstroIOconvert  ../j3AstroIOconvert.cpp )
set_property(TARGET j3AstroIOconvert PROPERTY CXX_STANDARD 11) 

target_include_directories(j3AstroIOconvert PRIVATE ${OpenCV_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/../ /usr/local/include )
target_link_libraries( j3AstroIOconvert PRIVATE j3AstroIO )


#install(TARGETS j3AstroIOraw
#  ARCHIVE DESTINATION lib
#  LIBRARY DESTINATION lib
#  PUBLIC_HEADER DESTINATION include)

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
