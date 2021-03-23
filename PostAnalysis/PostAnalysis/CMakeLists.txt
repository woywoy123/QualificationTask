########################################################################################
############################## Package: PostAnalysis ###################################
########################################################################################

# Setting the name of the package 
atlas_subdir(PostAnalysis)

# Minimum CMake version 
cmake_minimum_required( VERSION 3.2 FATAL_ERROR )

# Include ROOT for the analysis 
find_package( ROOT COMPONENTS REQUIRED Graf Gpad RooFitCore Hist RIO MathCore RooFit Threads Tree)

atlas_add_root_dictionary(PackageLib PostAnalysisDictSource
  ROOT_HEADERS PostAnalysis/*.h Root/LinkDef.h
  EXTERNAL_PACKAGES ROOT)


# Creating the library and linking the files
atlas_add_library( PackageLib
  PostAnalysis/*.h Root/*.cxx Root/*.h ${PostAnalysisDictSource}
  INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
  LINK_LIBRARIES ${ROOT_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

atlas_add_executable( PostAnalysis util/PostAnalysis.cxx
  INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
  LINK_LIBRARIES ${ROOT_LIBRARIES} PackageLib)
