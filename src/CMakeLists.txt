
### Create library etrans_${prec}.so

message( "ectrans_HAVE_SINGLE_PRECISION = ${ectrans_HAVE_SINGLE_PRECISION}")
message( "ectrans_HAVE_DOUBLE_PRECISION = ${ectrans_HAVE_DOUBLE_PRECISION}")
message( "ectrans_HAVE_FFTW = ${ectrans_HAVE_FFTW}")


if( ectrans_HAVE_SINGLE_PRECISION )
	list( APPEND precs sp )
endif()

if( ectrans_HAVE_DOUBLE_PRECISION )
	list( APPEND precs dp )
endif()

message( "precisions : ${precs}" )

ecbuild_list_add_pattern( LIST etrans_src
                          GLOB
                                biper/*/*.F90
                                etrans/*/*.F90
                                ifsaux/*/*.F90
                          QUIET
                        )

ecbuild_list_add_pattern( LIST etrans_hdr
                          GLOB
                                etrans/interface/*.h
                          QUIET
                        )
						
foreach ( prec IN LISTS precs )

	ecbuild_add_library( TARGET etrans_${prec}

	  PUBLIC_INCLUDES
				$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/etrans/interface>
				$<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
				$<INSTALL_INTERFACE:include>
				$<INSTALL_INTERFACE:include/etrans>

	  PUBLIC_LIBS  
				trans_gpu_static_CA_${prec}
				fiat
				parkind_${prec}

	  SOURCES
				${etrans_src}

      INSTALL_HEADERS_LIST
				${etrans_hdr}
				
      HEADER_DESTINATION include/etrans
		TYPE             STATIC
    )


	if( ectrans_HAVE_FFTW )
	  target_link_libraries( etrans_${prec} PRIVATE ${FFTW_LIBRARIES} )
	  target_include_directories( etrans_${prec} PRIVATE ${FFTW_INCLUDE_DIRS} )
	  target_compile_definitions( etrans_${prec} PRIVATE WITH_FFTW )
	endif()
	
	target_compile_definitions( etrans_${prec} PRIVATE USE_CUDA_AWARE_MPI_EFTDIR )

    set_target_properties( etrans_${prec} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/module/etrans_${prec} )
    target_include_directories( etrans_${prec} PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/module/etrans_${prec}> )

    target_include_directories(  etrans_${prec} PUBLIC $<INSTALL_INTERFACE:module/etrans_${prec}> )
    install( DIRECTORY ${CMAKE_BINARY_DIR}/module/etrans_${prec}/
	  	     DESTINATION module/etrans_${prec}
		     COMPONENT modules )
	  
endforeach()
