#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "uwlcm::uwlcm" for configuration "Release"
set_property(TARGET uwlcm::uwlcm APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(uwlcm::uwlcm PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/uwlcm"
  )

list(APPEND _IMPORT_CHECK_TARGETS uwlcm::uwlcm )
list(APPEND _IMPORT_CHECK_FILES_FOR_uwlcm::uwlcm "${_IMPORT_PREFIX}/bin/uwlcm" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
