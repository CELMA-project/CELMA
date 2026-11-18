add_executable(DDXOverJCylinder DDXOverJCylinder.cxx)
target_link_libraries(DDXOverJCylinder bout++::bout++)
set_target_properties(DDXOverJCylinder PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
