add_definitions(-DTIGHT)
add_definitions(-DADDPAPI)

add_executable(unsteady_tight_papi ${CMAKE_SOURCE_DIR}/unsteady.c
	${CMAKE_SOURCE_DIR}/lb.c
	${CMAKE_SOURCE_DIR}/boundaries.c)

target_link_libraries(unsteady_tight_papi m ${PAPILIB})
