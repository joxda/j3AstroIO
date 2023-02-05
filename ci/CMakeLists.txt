cmake_minimum_required(VERSION 3.11.99)

project(j3AstroIO LANGUAGES CXX VERSION 0.1.0)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/../cmake/")

include(GNUInstallDirs)

message(STATUS "THIS IS A CI RUN")
set(CMAKE_FIND_STATIC FIRST) 
set(BUILD_SHARED_LIBS=OFF)
SET(CMAKE_FIND_LIBRARY_SUFFIXES .a $(CMAKE_FIND_LIRBARY_SUFFIXES))
SET(LINK_DEPENDS_NO_SHARED=ON)

find_package( OpenCV REQUIRED )
find_package( LibRaw REQUIRED )
find_package( LibMagic REQUIRED )
find_package( Exiv2 REQUIRED )
find_package( CFITSIO REQUIRED )

set(CMAKE_CXX_STANDARD 11)

include_directories(embed)
add_subdirectory(embed)

embed_resources( someResources magic.mgc )

set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra -Wshadow -Wnon-virtual-dtor -Wunused -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wno-unused")

set(CPACK_PACKAGE_NAME ${PROJECT_NAME})
set(CPACK_GENERATOR "TGZ")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Joachim Janz") # required
set(CPACK_PACKAGE_VENDOR "j3")
set(CPACK_VERBATIM_VARIABLES ture)
set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
set (CPACK_RESOURCE_FILE_LICENSE  
    "${CMAKE_CURRENT_SOURCE_DIR}/../LICENSE")
    set(CPACK_PACKAGE_DESCRIPTION_FILE ${CMAKE_SOURCE_DIR}/../README.md)
set(CPACK_RESOURCE_FILE_README ${CMAKE_SOURCE_DIR}/../README.md)
include(CPack)

file(WRITE   ${CMAKE_BINARY_DIR}/CMakeFiles/test.cpp "#include <exiv2/exiv2.hpp>\nint main(){Exiv2::Image::UniquePtr EXimage;\nreturn 0;}")
try_compile(COMPILE_SUCCEEDED ${CMAKE_BINARY_DIR}/CMakeFiles/compile_tests  ${CMAKE_BINARY_DIR}/CMakeFiles/test.cpp CMAKE_FLAGS "-DINCLUDE_DIRECTORIES:STRING=/usr/local/include" CXX_STANDARD 11 OUTPUT_VARIABLE OUTPUT)

add_library( j3AstroIO  ../j3AstroIO.cpp ${someResources} )

target_compile_definitions( j3AstroIO PRIVATE CIRUN=1 )
if(COMPILE_SUCCEEDED)
  add_compile_definitions(UNIQ)
endif()

set_property(TARGET j3AstroIO PROPERTY CXX_STANDARD 11)
set_property(TARGET j3AstroIO PROPERTY PUBLIC_HEADER ../j3AstroIO.hpp)
target_include_directories( j3AstroIO SYSTEM PRIVATE ${OpenCV_INCLUDE_DIRS} ${CFITSIO_INCLUDES})
target_include_directories( j3AstroIO PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}> )

target_link_libraries( j3AstroIO  ${OpenCV_LIBS} ${LibRaw_LIBRARIES} ${LIBMAGIC_LIBRARY} ${CFITSIO_LIBRARIES} ${Exiv2_LIBRARIES})

#IF( (APPLE) )
#  message("APPLE CI")   
#  FIND_PACKAGE(ZLIB) ${ZLIB_LIBRARIES} 
#  FIND_PACKAGE(JPEG)
##  target_link_libraries( j3AstroIO  PRIVATE expat /usr/local/lib/libexiv2.a /usr/local/lib/libexiv2-xmp.a iconv /usr/local/lib/libpng.a bz2 curl cfitsio  /usr/local/lib/libmagic.a ${OpenCV_LIBS} /usr/local/lib/libraw.a /usr/local/lib/libjpeg.a )
#  target_link_libraries( j3AstroIO  PRIVATE exiv2 cfitsio  magic ${OpenCV_LIBS} raw )
#  install(FILES /usr/local/opt/libmagic/share/misc/magic.mgc DESTINATION .)
#ELSE() 
#  message("LINUX CI")
#  target_link_libraries( j3AstroIO PRIVATE ${OpenCV_LIBS} /usr/local/lib/libraw.a cfitsio  magic /usr/local/lib/libexiv2.a /usr/local/lib/libexiv2-xmp.a /usr/lib/x86_64-linux-gnu/libexpat.a pthread )
#ENDIF()

set_target_properties(j3AstroIO PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

install(TARGETS j3AstroIO
    EXPORT j3AstroIO-targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(EXPORT j3AstroIO-targets
    FILE j3AstroIO-targets.cmake
    NAMESPACE j3AstroIO::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/j3AstroIO
)

add_library(j3AstroIO::j3AstroIO STATIC IMPORTED)

set_target_properties(j3AstroIO::j3AstroIO PROPERTIES
    CMAKE_CXX_STANDARD 11
    IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/lib
)

include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${CMAKE_SOURCE_DIR}/../cmake/j3AstroIO-config.cmake.in
    ${CMAKE_BINARY_DIR}/cmake/j3AstroIO-config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/j3AstroIO
)

write_basic_package_version_file(
    ${CMAKE_BINARY_DIR}/cmake/j3AstroIO-config-version.cmake
    VERSION ${j3AstroIO_VERSION}
    COMPATIBILITY AnyNewerVersion
)

install(
    FILES
        ${CMAKE_BINARY_DIR}/cmake/j3AstroIO-config.cmake
        ${CMAKE_BINARY_DIR}/cmake/j3AstroIO-config-version.cmake
        ${CMAKE_MODULE_PATH}/FindExiv2.cmake
        ${CMAKE_MODULE_PATH}/FindCFITSIO.cmake
        ${CMAKE_MODULE_PATH}/FindLibMagic.cmake
        ${CMAKE_MODULE_PATH}/FindLibRaw.cmake
        ${CMAKE_MODULE_PATH}/LibFindMacros.cmake
        ${CMAKE_MODULE_PATH}/FindExiv2.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/j3AstroIO
)

install(
    DIRECTORY
        ${CMAKE_CURRENT_SOURCE_DIR}/../licenses
    DESTINATION ${CMAKE_INSTALL_DATADIR}/j3AstroIO
)



# To install, for example, MSVC runtime libraries:
#include(InstallRequiredSystemLibraries)

# To install other/non-system 3rd party required libraries:
#configure_file(
#  ${CMAKE_MODULE_PATH}/FixBundle.cmake.in
#  ${CMAKE_CURRENT_BINARY_DIR}/FixBundle.cmake
#  @ONLY
#  )

#install(SCRIPT ${CMAKE_CURRENT_BINARY_DIR}/FixBundle.cmake)


export(EXPORT j3AstroIO-targets
    FILE ${CMAKE_BINARY_DIR}/cmake/j3AstroIO-targets.cmake
    NAMESPACE j3AstroIO::
)

add_executable(j3AstroIOconvert ${CMAKE_CURRENT_SOURCE_DIR}/../j3AstroIOconvert.cpp)

if(COMPILE_SUCCEEDED)
  target_compile_definitions(j3AstroIOconvert PUBLIC UNIQ)
endif()

target_include_directories(j3AstroIOconvert SYSTEM PRIVATE ${OpenCV_INCLUDE_DIRS})
target_link_libraries(j3AstroIOconvert PRIVATE j3AstroIO ${OpenCV_LIBS})

include(GetPrerequisites)

set(BINARY_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/j3AstroIOconvert")
get_prerequisites(${BINARY_LOCATION} DEPENDENCIES 0 0 "" "")
message(${DEPENDENCIES})
list(FILTER DEPENDENCIES INCLUDE REGEX ".libopencv*")
message(${OpenCV_INSTALL_PATH})
foreach(DEPENDENCY_FILE ${DEPENDENCIES})
  get_filename_component(DEPENDENCY_NAME "${DEPENDENCY_FILE}" NAME)
  get_filename_component(DEPENDENCY "${DEPENDENCY_FILE}" REALPATH)
  get_filename_component(DEPENDENCY_PATH "${DEPENDENCY}" DIRECTORY)
  get_filename_component(DEP_FILE ${OpenCV_INSTALL_PATH}/lib//${DEPENDENCY_NAME} REALPATH)
  message(${DEP_FILE})
  install(FILES ${DEP_FILE} DESTINATION "${CMAKE_INSTALL_LIBDIR}" )
endforeach()

install(
    DIRECTORY
      ${OpenCV_INCLUDE_DIRS}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
