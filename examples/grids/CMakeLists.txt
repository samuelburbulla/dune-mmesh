set(FILES
  cube2d.dgf
  cube2d.geo
  cube2d.msh
  cube2dinterval.dgf
  cube3d.geo
  cube3d.msh
  horizontal2d.geo
  horizontal2d.msh
  horizontal3d.geo
  horizontal3d.msh
)

foreach (file ${FILES})
  configure_file(${file} ${file} COPYONLY)
endforeach()
