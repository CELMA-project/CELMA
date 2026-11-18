add_executable(DDXCylinder DDXCylinder.cxx)
target_link_libraries(DDXCylinder bout++::bout++)
set_target_properties(DDXCylinder PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
