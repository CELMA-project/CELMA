add_executable(DDZCylinder DDZCylinder.cxx)
target_link_libraries(DDZCylinder bout++::bout++)
set_target_properties(DDZCylinder PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
