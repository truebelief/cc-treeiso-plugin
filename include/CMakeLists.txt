
target_sources( ${PROJECT_NAME}
	PRIVATE
		${CMAKE_CURRENT_LIST_DIR}/ccTreeIsoDlg.h
		${CMAKE_CURRENT_LIST_DIR}/block.hpp
		${CMAKE_CURRENT_LIST_DIR}/maxflow.hpp
		${CMAKE_CURRENT_LIST_DIR}/cp_d0_dist.hpp
		${CMAKE_CURRENT_LIST_DIR}/cut_pursuit_d0.hpp
		${CMAKE_CURRENT_LIST_DIR}/cut_pursuit.hpp
		${CMAKE_CURRENT_LIST_DIR}/knncpp.h
		${CMAKE_CURRENT_LIST_DIR}/qTreeIso.h
		${CMAKE_CURRENT_LIST_DIR}/qTreeIsoCommands.h		
		${CMAKE_CURRENT_LIST_DIR}/TreeIso.h
		${CMAKE_CURRENT_LIST_DIR}/TreeIsoHelper.hpp
)

target_include_directories( ${PROJECT_NAME}
	PRIVATE
		${CMAKE_CURRENT_SOURCE_DIR}
)
