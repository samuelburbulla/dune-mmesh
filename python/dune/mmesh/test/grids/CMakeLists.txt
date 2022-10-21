set(GRIDFILES
  line
  plane
  tjunction
)

foreach(test ${GRIDFILES})
 add_python_targets(test ${GRIDFILES})
endforeach()
