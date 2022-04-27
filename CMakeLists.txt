############################################################################################
# etrans project
# Collection of stripped down libraries from Aladin
# - etrans

cmake_minimum_required( VERSION 3.12 FATAL_ERROR )
find_package( ecbuild 3.4 REQUIRED HINTS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../ecbuild )

project( etrans LANGUAGES C Fortran CXX )

#ecbuild_enable_fortran( REQUIRED NO_MODULE_DIRECTORY )

ecbuild_enable_fortran( REQUIRED MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/module )

# dependencies:
ecbuild_find_package( fiat REQUIRED )
ecbuild_find_package( ectrans REQUIRED )
#find_package( OpenACC REQUIRED )

#enable_language( CUDA )
#set( OpenACC_Fortran_FLAGS "-acc" )
#set( OpenACC_Fortran_FLAGS "${OpenACC_Fortran_FLAGS} -ta=tesla:lineinfo,deepcopy,fastmath" )
#set( OpenACC_Fortran_FLAGS "${OpenACC_Fortran_FLAGS} -Mcuda" )

############################################################################################
# build etrans library

add_subdirectory( src )

ecbuild_add_option( FEATURE DRIVER
                    DEFAULT ON
                    DESCRIPTION "Compile driver program for LAM transforms"
				)

#if( HAVE_DRIVER )
#  add_subdirectory( driver )
#endif()

ecbuild_install_project( NAME etrans )

ecbuild_print_summary()
