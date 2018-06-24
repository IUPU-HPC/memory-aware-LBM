add_definitions(-DPANEL)
add_definitions(-DADDPAPI)

add_executable(unsteady_panel_papi ${CMAKE_SOURCE_DIR}/unsteady.c
	${CMAKE_SOURCE_DIR}/lb.c
	${CMAKE_SOURCE_DIR}/boundaries.c)

target_link_libraries(unsteady_panel_papi m ${PAPILIB})
