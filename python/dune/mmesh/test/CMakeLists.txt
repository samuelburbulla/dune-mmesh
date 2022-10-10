add_subdirectory(grids)

add_python_targets(test __main__)

set(TESTS
  adaptation
  interfaceindicator
  normals
  distance
  laplace
  monolithicsolve
)

foreach(test ${TESTS})
  dune_python_add_test(NAME ${test}
                       SCRIPT ${test}.py
                       WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                       LABELS python)

 add_python_targets(test ${test})
endforeach()
