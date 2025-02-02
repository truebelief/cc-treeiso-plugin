cmake_minimum_required( VERSION 3.10 )

option( PLUGIN_STANDARD_QTREEISO "Install qTreeIso plugin" OFF )

if( PLUGIN_STANDARD_QTREEISO )
	
	project( QTREEISO_PLUGIN )
	
	AddPlugin( NAME ${PROJECT_NAME} )	

	find_package(Eigen3 QUIET)
	if(NOT Eigen3_FOUND)
		set( EIGEN_ROOT_DIR "" CACHE PATH "Eigen root (contains the Eigen directory)" )
		if ( NOT EIGEN_ROOT_DIR )
			message( SEND_ERROR "No Eigen root directory specified (EIGEN_ROOT_DIR)" )
		endif()
		include_directories(${EIGEN_ROOT_DIR})		
	else()
		target_link_libraries(${PROJECT_NAME} Eigen3::Eigen)
	endif()

	add_subdirectory( include )
	add_subdirectory( src )
	add_subdirectory( ui )

endif()
