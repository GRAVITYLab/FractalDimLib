cmake_minimum_required (VERSION 2.8)
project( FractalDimLib )

# Use MPI compiler and linker
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})

# Use mpicxx to replace 
SET( CMAKE_C_COMPILER mpicc )
SET( CMAKE_CXX_COMPILER mpicxx )

# Set ITL path
set(ITL_DIR "../.." CACHE PATH "The path to ITL")

# Set OSUFlow path
set( OSUFLOW_INSTALL_DIR "/home/abon/install/osuflow/osuflow" CACHE PATH "Installed path of osuflow" )

# Set OpenCV path
set( OPENCV_DIR "/home/abon/install/opencv-2.2.0" CACHE PATH "Path to opencv" )

# Set OpenCV path
set( MY_INSTALL_DIR "/home/abon/code/chaudhua/install" CACHE PATH "Path to my installed libraries" )


# Include header file directory to the search path
include_directories ( "${OPENCV_DIR}/include/opencv" )
include_directories ( "${OPENCV_DIR}/include/opencv2" )
include_directories ( "${OPENCV_DIR}/include" )
#####
include_directories ( "${OSUFLOW_INSTALL_DIR}/include" )
#####
include_directories( "${ITL_DIR}/include" )
#####
include_directories( "${MY_INSTALL_DIR}/FlowVisLib/include" )
#####
include_directories ( "${PROJECT_SOURCE_DIR}/include" )

add_library( FractalDimLib STATIC 
                  ${PROJECT_SOURCE_DIR}/src/FDL_streamline.cpp
				  ${PROJECT_SOURCE_DIR}/src/FDL_osuflow.cpp	
				  ${PROJECT_SOURCE_DIR}/src/FDL_boxcountdim.cpp
				  ${PROJECT_SOURCE_DIR}/src/FDL_curvature.cpp
				  ${PROJECT_SOURCE_DIR}/src/FDL_featuremap.cpp
				  ${PROJECT_SOURCE_DIR}/src/FDL_featurehistogram.cpp
				  ${PROJECT_SOURCE_DIR}/src/FDL_algorithm_bcd.cpp
				  ${PROJECT_SOURCE_DIR}/src/FDL_algorithm_curvature.cpp )
     
# Installation Instructions
# lib is interpreted as relative to CMAKE_INSTALL_PREFIX
FILE(GLOB headers "${PROJECT_SOURCE_DIR}/include/*.h")
install( FILES ${headers} DESTINATION include )
install( TARGETS FractalDimLib DESTINATION lib )
