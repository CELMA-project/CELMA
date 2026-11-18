add_executable(VolIntCyl VolIntCyl.cxx)
target_link_libraries(VolIntCyl bout++::bout++)
set_target_properties(VolIntCyl PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
