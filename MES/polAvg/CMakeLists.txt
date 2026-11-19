add_executable(polAvg polAvg.cxx)
target_link_libraries(polAvg bout++::bout++)
set_target_properties(polAvg PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
