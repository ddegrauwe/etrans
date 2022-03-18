
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

foreach ( prec IN LISTS precs )

	ecbuild_add_library( TARGET etrans_${prec}

	  PUBLIC_INCLUDES
				$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/etrans/interface>
				$<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
				$<INSTALL_INTERFACE:include>
				$<INSTALL_INTERFACE:include/etrans>

	  PUBLIC_LIBS  
				trans_${prec}
				fiat
				parkind_${prec}

	  SOURCES
				biper/interface/etibihie.h
				biper/interface/fpbipere.h
				biper/interface/horiz_field.h
				biper/external/etibihie.F90
				biper/external/fpbipere.F90
				biper/external/horiz_field.F90
				biper/module/espline_mod.F90
				biper/module/esmoothe_mod.F90
				biper/module/extper_mod.F90
				biper/module/ewindowe_mod.F90
				etrans/external/egath_grid.F90
				etrans/external/edist_spec.F90
				etrans/external/egpnorm_trans.F90
				etrans/external/esetup_trans.F90
				etrans/external/egath_spec.F90
				etrans/external/etrans_inq.F90
				etrans/external/etrans_end.F90
				etrans/external/etrans_release.F90
				etrans/external/edist_grid.F90
				etrans/external/edir_trans.F90
				etrans/external/einv_trans.F90
				etrans/external/especnorm.F90
				etrans/interface/einv_trans.h
				etrans/interface/edir_trans.h
				etrans/interface/edist_spec.h
				etrans/interface/einv_transad.h
				etrans/interface/especnorm.h
				etrans/interface/esetup_trans.h
				etrans/interface/etrans_inq.h
				etrans/interface/edir_transad.h
				etrans/interface/egath_grid.h
				etrans/interface/edist_grid.h
				etrans/interface/egath_spec.h
				etrans/interface/etrans_release.h
				etrans/interface/etrans_end.h
				etrans/interface/egpnorm_trans.h
				etrans/module/eupdspb_mod.F90
				etrans/module/easre1b_mod.F90
				etrans/module/efsc_mod.F90
				etrans/module/eftinv_ctl_mod.F90
				etrans/module/suemp_trans_preleg_mod.F90
				etrans/module/egath_spec_control_mod.F90
				etrans/module/evdtuv_mod.F90
				etrans/module/eprfi1b_mod.F90
				etrans/module/suefft_mod.F90
				etrans/module/edist_spec_control_mod.F90
				etrans/module/einv_trans_ctl_mod.F90
				etrans/module/tpmald_geo.F90
				etrans/module/cpl_int_mod.F90
				etrans/module/eftdir_ctl_mod.F90
				etrans/module/esetup_geom_mod.F90
				etrans/module/tpmald_tcdis.F90
				etrans/module/eset_resol_mod.F90
				etrans/module/eprfi2_mod.F90
				etrans/module/tpmald_fields.F90
				etrans/module/euvtvd_comm_mod.F90
				etrans/module/suemp_trans_mod.F90
				etrans/module/espnormc_mod.F90
				etrans/module/euvtvd_mod.F90
				etrans/module/eltdir_ctl_mod.F90
				etrans/module/eledir_mod.F90
				etrans/module/eltdir_mod.F90
				etrans/module/suemplatb_mod.F90
				etrans/module/esetup_dims_mod.F90
				etrans/module/tpmald_distr.F90
				etrans/module/eltinv_ctl_mod.F90
				etrans/module/tpmald_fft.F90
				etrans/module/eupdsp_mod.F90
				etrans/module/tpmald_dim.F90
				etrans/module/espnormd_mod.F90
				etrans/module/eltinv_mod.F90
				etrans/module/edealloc_resol_mod.F90
				etrans/module/eleinv_mod.F90
				etrans/module/eprfi1_mod.F90
				etrans/module/edir_trans_ctl_mod.F90
				etrans/module/espnorm_ctl_mod.F90
				etrans/module/suemplat_mod.F90
				etrans/module/eprfi2b_mod.F90
				etrans/module/suestaonl_mod.F90
				etrans/module/espnsde_mod.F90
				ifsaux/fa/ellips.F90

      INSTALL_HEADERS_LIST
				etrans/interface/etrans_end.h
				etrans/interface/etrans_release.h
				etrans/interface/esetup_trans.h
				etrans/interface/edist_grid.h
				etrans/interface/egath_spec.h
				etrans/interface/edir_trans.h
				etrans/interface/edist_spec.h
				etrans/interface/etrans_inq.h
				etrans/interface/egath_grid.h
				etrans/interface/einv_trans.h
				etrans/interface/egpnorm_trans.h
				etrans/interface/especnorm.h
				
      HEADER_DESTINATION include/etrans
    )


	if( ectrans_HAVE_FFTW )
	  target_link_libraries( etrans_${prec} PRIVATE ${FFTW_LIBRARIES} )
	  target_include_directories( etrans_${prec} PRIVATE ${FFTW_INCLUDE_DIRS} )
	  target_compile_definitions( etrans_${prec} PRIVATE WITH_FFTW )
	endif()

    set_target_properties( etrans_${prec} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/module/etrans_${prec} )
    target_include_directories( etrans_${prec} PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/module/etrans_${prec}> )

    target_include_directories(  etrans_${prec} PUBLIC $<INSTALL_INTERFACE:module/etrans_${prec}> )
    install( DIRECTORY ${CMAKE_BINARY_DIR}/module/etrans_${prec}/
	  	     DESTINATION module/etrans_${prec}
		     COMPONENT modules )
	  
endforeach()