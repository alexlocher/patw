cmake_minimum_required (VERSION 2.8.9) 
project (PATW)

list(APPEND CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/thirdparty")
#set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "${CMAKE_SOURCE_DIR}/thirdparty/lib" )
#set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "${CMAKE_SOURCE_DIR}/thirdparty/include" )

#message("MODULE PATH: ${CMAKE_MODULE_PATH}")

ADD_SUBDIRECTORY(lib/hpmvs)
ADD_SUBDIRECTORY(lib/yaply)

set(BUILD_TESTS OFF CACHE BOOL "enable tests of  Sophus")
ADD_SUBDIRECTORY(lib/Sophus)

# get the required variables
get_directory_property(HPMVS_LIBRARIES DIRECTORY lib/hpmvs DEFINITION HPMVS_LIBRARIES)
list(APPEND HPMVS_LIBRARIES HPMVS_static)
get_directory_property(HPMVS_INCLUDE_DIRS DIRECTORY lib/hpmvs DEFINITION HPMVS_INCLUDE_DIRS)
set(SOPHUS_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/lib/Sophus )
set(YAPLY_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/lib/yaply/include )


find_package(PCL 1.8 REQUIRED COMPONENTS common octree)
find_package(Theia REQUIRED )

message("HPMVS libraries: ${HPMVS_INCLUDE_DIRS}, libs: ${HPMVS_LIBRARIES}")
message("Sophus header: ${SOPHUS_INCLUDE_DIRS}, libs: ${Sophus_LIBRARIES}")


include_directories( 
	${PROJECT_SOURCE_DIR}/include
	${THEIA_INCLUDE_DIRS}
	${HPMVS_INCLUDE_DIRS}
	${YAPLY_INCLUDE_DIRS}
	${SOPHUS_INCLUDE_DIRS} 
	${PCL_INCLUDE_DIRS}
)

ADD_DEFINITIONS(
    -std=c++11
    -O3
)


add_executable (patw src/patw.cpp)
target_link_libraries(patw ${HPMVS_LIBRARIES} ${THEIA_LIBRARIES} ${PCL_COMMON_LIBRARIES} ${PCL_OCTREE_LIBRARIES})

add_executable (nvm2ply src/nvm2ply.cpp)
target_link_libraries(nvm2ply ${HPMVS_LIBRARIES} ${THEIA_LIBRARIES} )

install(TARGETS patw nvm2ply RUNTIME DESTINATION bin)
