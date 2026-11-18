add_executable(jParSheath jParSheath.cxx)
target_link_libraries(jParSheath bout++::bout++)
set_target_properties(jParSheath PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
