# build driver program for LAM spectral transforms

set( HAVE_dp ${ectrans_HAVE_DOUBLE_PRECISION} )
set( HAVE_sp ${ectrans_HAVE_SINGLE_PRECISION} )

foreach( prec sp dp )
	if( HAVE_${prec} )
	  ecbuild_add_executable(TARGET  driver_lam_${prec}
							 SOURCES driver_lam.F90
							 INCLUDES
								  ${MPI_Fortran_INCLUDE_PATH}
							 LIBS 
								  etrans_${prec}
								  fiat parkind_${prec}
								  ${MPI_Fortran_LIBRARIES}
								  trans_gpu_static_${prec}
								  gpu
								  OpenACC::OpenACC_Fortran 
								  ${LAPACK_LIBRARIES}
								  nvToolsExt
								  
						  )

	  #target_link_libraries( driver-spectrans PRIVATE OpenACC::OpenACC_Fortran )
	  set_property( TARGET driver_lam_${prec} PROPERTY CUDA_ARCHITECTURES 70 )
	  target_compile_options( driver_lam_${prec} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-g -acc -Minfo=acc -gpu=cc70,lineinfo,deepcopy,fastmath,nordc -cudalib=cufft,cublas -fpic>  )
	  set_target_properties(driver_lam_${prec} PROPERTIES LINK_FLAGS "-acc -cudalib=cufft,cublas -fpic")
	  if ( NVTX_ROOT )
		target_link_libraries( driver_lam_${prec} "${NVTX_ROOT}/lib/libnvtx.a" )
	  endif()
	  message("Building ${prec} GPU driver")
	endif()
endforeach()
