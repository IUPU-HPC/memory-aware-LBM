add_executable(unsteady_origin ${CMAKE_SOURCE_DIR}/unsteady.c
	${CMAKE_SOURCE_DIR}/lb.c
	${CMAKE_SOURCE_DIR}/boundaries.c)

target_link_libraries(unsteady_origin m)
