#----------------------------------------------------------------
# Generated CMake target import file for configuration "None".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "gmsh::shared" for configuration "None"
set_property(TARGET gmsh::shared APPEND PROPERTY IMPORTED_CONFIGURATIONS NONE)
set_target_properties(gmsh::shared PROPERTIES
  IMPORTED_LOCATION_NONE "${_IMPORT_PREFIX}/lib/loongarch64-linux-gnu/libgmsh.so.4.15.2"
  IMPORTED_SONAME_NONE "libgmsh.so.4.15"
  )

list(APPEND _cmake_import_check_targets gmsh::shared )
list(APPEND _cmake_import_check_files_for_gmsh::shared "${_IMPORT_PREFIX}/lib/loongarch64-linux-gnu/libgmsh.so.4.15.2" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
