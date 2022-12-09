add_subdirectory(grids)

add_python_targets(test __main__)

# Add one parallel test without dune-fem
dune_python_add_test(NAME grid
                     SCRIPT grid.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     MPI_RANKS 1 2 4
                     TIMEOUT 100000
                     LABELS python)
add_python_targets(test grid)

set(TESTS
  adaptation
  coupledsolve
  boundary
  interfaceindicator
  normals
  distance
  solve
)

foreach(test ${TESTS})
  dune_python_add_test(NAME ${test}
                       SCRIPT ${test}.py
                       WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                       LABELS python)

 add_python_targets(test ${test})
endforeach()


if(dune-fem_FOUND)

  # Python tests with dune-fem
  set(TESTS
    adaptation
    coupledsolve
    boundary
    interfaceindicator
    normals
    distance
    monolithicsolve
  )

  foreach(test ${TESTS})
    dune_python_add_test(NAME ${test}
                         SCRIPT ${test}.py
                         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                         LABELS python)

   add_python_targets(test ${test})
  endforeach()

  # Parallel tests
  set(PARALLEL_TESTS
    laplace
  )

  foreach(test ${PARALLEL_TESTS})
    dune_python_add_test(NAME ${test}
                         SCRIPT ${test}.py
                         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                         MPI_RANKS 1 2 4
                         TIMEOUT 100000
                         LABELS python)

   add_python_targets(test ${test})
  endforeach()

endif() # dune-fem_FOUND
