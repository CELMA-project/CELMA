add_executable(divSource divSource.cxx)
target_link_libraries(divSource bout++::bout++)
set_target_properties(divSource PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
