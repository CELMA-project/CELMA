add_executable(celmaWBA celmaWBA.cxx)
target_link_libraries(celmaWBA bout++::bout++)
set_target_properties(celmaWBA PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
