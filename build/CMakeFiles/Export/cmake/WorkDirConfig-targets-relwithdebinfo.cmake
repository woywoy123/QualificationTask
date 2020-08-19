#----------------------------------------------------------------
# Generated CMake target import file for configuration "RelWithDebInfo".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "WorkDir::PackageLib" for configuration "RelWithDebInfo"
set_property(TARGET WorkDir::PackageLib APPEND PROPERTY IMPORTED_CONFIGURATIONS RELWITHDEBINFO)
set_target_properties(WorkDir::PackageLib PROPERTIES
  IMPORTED_LOCATION_RELWITHDEBINFO "${_IMPORT_PREFIX}/lib/libPackageLib.so"
  IMPORTED_SONAME_RELWITHDEBINFO "libPackageLib.so"
  )

list(APPEND _IMPORT_CHECK_TARGETS WorkDir::PackageLib )
list(APPEND _IMPORT_CHECK_FILES_FOR_WorkDir::PackageLib "${_IMPORT_PREFIX}/lib/libPackageLib.so" )

# Import target "WorkDir::PostAnalysis" for configuration "RelWithDebInfo"
set_property(TARGET WorkDir::PostAnalysis APPEND PROPERTY IMPORTED_CONFIGURATIONS RELWITHDEBINFO)
set_target_properties(WorkDir::PostAnalysis PROPERTIES
  IMPORTED_LOCATION_RELWITHDEBINFO "${_IMPORT_PREFIX}/bin/PostAnalysis"
  )

list(APPEND _IMPORT_CHECK_TARGETS WorkDir::PostAnalysis )
list(APPEND _IMPORT_CHECK_FILES_FOR_WorkDir::PostAnalysis "${_IMPORT_PREFIX}/bin/PostAnalysis" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
