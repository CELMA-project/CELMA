add_executable(divPerp divPerp.cxx)
target_link_libraries(divPerp bout++::bout++)
set_target_properties(divPerp PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
