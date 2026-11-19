add_executable(D2DX2Cylinder D2DX2Cylinder.cxx)
target_link_libraries(D2DX2Cylinder bout++::bout++)
set_target_properties(D2DX2Cylinder PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
