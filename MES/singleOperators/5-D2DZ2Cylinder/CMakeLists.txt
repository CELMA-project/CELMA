add_executable(D2DZ2Cylinder D2DZ2Cylinder.cxx)
target_link_libraries(D2DZ2Cylinder bout++::bout++)
set_target_properties(D2DZ2Cylinder PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
