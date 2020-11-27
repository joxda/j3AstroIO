project(j3AstroIO)


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

set(CMAKE_CXX_STANDARD 11)

file(WRITE   ${CMAKE_BINARY_DIR}/CMakeFiles/test.cpp "#include <exiv2/exiv2.hpp>\nint main(){Exiv2::Image::UniquePtr EXimage;\nreturn 0;}")
try_compile(COMPILE_SUCCEEDED ${CMAKE_BINARY_DIR}/CMakeFiles/compile_tests  ${CMAKE_BINARY_DIR}/CMakeFiles/test.cpp CXX_STANDARD 11)

if(COMPILE_SUCCEEDED)
  add_compile_definitions(UNIQ)
endif()


add_library( j3AstroIO  j3AstroIO.cpp )
set_property(TARGET j3AstroIO PROPERTY CXX_STANDARD 11)
set_property(TARGET j3AstroIO PROPERTY PUBLIC_HEADER j3AstroIO.hpp)

target_include_directories(j3AstroIO  PRIVATE ${OpenCV_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR} /usr/local/include )

IF( TARGET Eigen3::Eigen )
 target_link_libraries( j3AstroIO  ${OpenCV_LIBS} raw magic Eigen3::Eigen cfitsio  magic exiv2)
ELSE()
 target_link_libraries( j3AstroIO  ${OpenCV_LIBS} raw magic cfitsio  magic exiv2)
ENDIF()
install(TARGETS j3AstroIO  ARCHIVE  DESTINATION lib  PUBLIC_HEADER DESTINATION include)
#install(TARGETS j3AstroIO  DESTINATION ${CMAKE_INSTALL_LIBDIR})
#install(TARGETS j3AstroIO  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

add_executable( j3AstroIOconvert  j3AstroIOconvert.cpp )
target_include_directories(j3AstroIOconvert PRIVATE ${OpenCV_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/../ /usr/local/include )
target_link_libraries( j3AstroIOconvert PRIVATE j3AstroIO  )

install(TARGETS j3AstroIOconvert DESTINATION bin)


#install(TARGETS j3AstroIO
#  ARCHIVE DESTINATION lib
#  LIBRARY DESTINATION lib
#  PUBLIC_HEADER DESTINATION include)

