add_executable(JTimesDivPerp JTimesDivPerp.cxx)
target_link_libraries(JTimesDivPerp bout++::bout++)
set_target_properties(JTimesDivPerp PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
