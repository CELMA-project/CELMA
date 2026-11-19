add_executable(yExtrapolation yExtrapolation.cxx)
target_link_libraries(yExtrapolation bout++::bout++)
set_target_properties(yExtrapolation PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
