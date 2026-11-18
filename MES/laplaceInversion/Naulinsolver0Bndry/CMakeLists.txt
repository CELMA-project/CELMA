add_executable(NaulinSolver NaulinSolver.cxx)
target_link_libraries(NaulinSolver bout++::bout++)
set_target_properties(NaulinSolver PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
