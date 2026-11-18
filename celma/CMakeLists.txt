add_executable(celma celma.cxx)
target_link_libraries(celma bout++::bout++)
set_target_properties(celma PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
