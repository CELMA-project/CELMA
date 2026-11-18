add_executable(radialLowPass radialLowPass.cxx)
target_link_libraries(radialLowPass bout++::bout++)
set_target_properties(radialLowPass PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
