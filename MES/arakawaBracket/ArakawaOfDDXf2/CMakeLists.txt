add_executable(ArakawaOfDDX ArakawaOfDDX.cxx)
target_link_libraries(ArakawaOfDDX bout++::bout++)
set_target_properties(ArakawaOfDDX PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
