set(GRIDFILES
  line
  plane
  tjunction
  vertical
)

foreach(test ${GRIDFILES})
 add_python_targets(test ${GRIDFILES})
endforeach()
