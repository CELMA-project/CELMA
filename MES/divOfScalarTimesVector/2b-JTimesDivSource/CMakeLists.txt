add_executable(JTimesDivSource JTimesDivSource.cxx)
target_link_libraries(JTimesDivSource bout++::bout++)
set_target_properties(JTimesDivSource PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
