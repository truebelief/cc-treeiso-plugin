cmake_minimum_required( VERSION 3.10 )

option( PLUGIN_STANDARD_QTREEISO "Install qTreeIso plugin" OFF )

if( PLUGIN_STANDARD_QTREEISO )
	project( QTREEISO_PLUGIN )
	
	find_package(Boost 1.67 REQUIRED COMPONENTS graph)
	find_package(Eigen3 REQUIRED)	

	add_subdirectory( include )
	add_subdirectory( src )
	add_subdirectory( ui )

	target_link_libraries(${PROJECT_NAME} Boost::graph Eigen3::Eigen)
endif()
