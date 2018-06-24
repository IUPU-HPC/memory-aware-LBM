add_definitions(-DTIGHT_BLOCK)

add_executable(unsteady_tight_block ${CMAKE_SOURCE_DIR}/unsteady.c
	${CMAKE_SOURCE_DIR}/lb.c
	${CMAKE_SOURCE_DIR}/boundaries.c)

target_link_libraries(unsteady_tight_block m)
