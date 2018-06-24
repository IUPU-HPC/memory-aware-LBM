add_definitions(-DCOMBINE)
add_definitions(-DADDPAPI)

add_executable(unsteady_combine_papi ${CMAKE_SOURCE_DIR}/unsteady.c
	${CMAKE_SOURCE_DIR}/lb.c
	${CMAKE_SOURCE_DIR}/boundaries.c)

target_link_libraries(unsteady_combine_papi m ${PAPILIB})
