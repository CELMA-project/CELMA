add_executable(DDXCylinderHighResCenter DDXCylinderHighResCenter.cxx)
target_link_libraries(DDXCylinderHighResCenter bout++::bout++)
set_target_properties(DDXCylinderHighResCenter PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
