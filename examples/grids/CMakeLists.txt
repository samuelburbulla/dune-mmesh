set(FILES
  cube2d.dgf
  cube2dinterval.dgf
  cube2d.msh
  cube3d.msh
  circle.msh
)

foreach (file ${FILES})
  configure_file(${file} ${file} COPYONLY)
endforeach()
