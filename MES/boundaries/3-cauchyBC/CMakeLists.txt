add_executable(cauchyBC cauchyBC.cxx)
target_link_libraries(cauchyBC bout++::bout++)
set_target_properties(cauchyBC PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
