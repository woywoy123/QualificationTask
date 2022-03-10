# Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
#
# This file describes to CTest how to build an entire release, and upload
# the build results to a CDash server.
#

# Set some basic properties for the build:
set( CTEST_SOURCE_DIRECTORY /home/tnom6927/QualificationTask/PostAnalysis )
set( CTEST_BINARY_DIRECTORY /home/tnom6927/QualificationTask/build )

set( CTEST_CUSTOM_MAXIMUM_NUMBER_OF_ERRORS   1000000 )
set( CTEST_CUSTOM_MAXIMUM_NUMBER_OF_WARNINGS 1000000 )

set( CTEST_COMMAND "ctest -D Experimental" )
set( CTEST_CMAKE_COMMAND "cmake" )

set( CTEST_USE_LAUNCHERS 1 )

set( CTEST_CMAKE_GENERATOR "Unix Makefiles" )
set( CTEST_BUILD_CONFIGURATION RelWithDebInfo )
if( ${CTEST_CMAKE_GENERATOR} STREQUAL "Unix Makefiles" )
   set( CTEST_BUILD_FLAGS -j8 )
endif()

# Identify the build:
set( CTEST_BUILD_NAME x86_64-centos7-gcc62-opt )
set( CTEST_SITE localhost.localdomain )

# Start an experimental build:
ctest_start( Experimental )

# Upload the project description file to CDash:
ctest_submit( FILES ${CTEST_BINARY_DIRECTORY}/Project.xml )

# Execute the project configuration:
ctest_configure( BUILD ${CTEST_BINARY_DIRECTORY}
   SOURCE ${CTEST_SOURCE_DIRECTORY}
   OPTIONS "-DCTEST_USE_LAUNCHERS=${CTEST_USE_LAUNCHERS}" )

# Execute the project build:
ctest_build( BUILD ${CTEST_BINARY_DIRECTORY} )

# Execute the test(s) of project:
ctest_test( BUILD ${CTEST_BINARY_DIRECTORY} )

# Submit all the results all at once:
ctest_submit()
