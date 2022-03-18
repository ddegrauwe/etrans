############################################################################################
# etrans project
# Collection of stripped down libraries from Aladin
# - etrans

cmake_minimum_required( VERSION 3.12 FATAL_ERROR )
find_package( ecbuild 3.4 REQUIRED HINTS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../ecbuild )

project( etrans LANGUAGES C Fortran )

#ecbuild_enable_fortran( REQUIRED NO_MODULE_DIRECTORY )

ecbuild_enable_fortran( REQUIRED MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/module )

# dependencies:
ecbuild_find_package( fiat REQUIRED )
ecbuild_find_package( ectrans REQUIRED )

############################################################################################
# build etrans library

add_subdirectory( src )

ecbuild_install_project( NAME etrans )

ecbuild_print_summary()
