add_executable(uEParSheath uEParSheath.cxx)
target_link_libraries(uEParSheath bout++::bout++)
set_target_properties(uEParSheath PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
