set(GRIDFILES
  line
  plane
)

foreach(test ${GRIDFILES})
 add_python_targets(test ${GRIDFILES})
endforeach()
