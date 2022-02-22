add_python_targets(mmesh
  __init__
  _grids
  _move
  _skeletontrace
  _solve
  utility
)

set(PRECOMPILE
  _grid2d
  _grid3d
  _interfacegrid2d
  _interfacegrid3d
  _interfaceindicator2d
  _interfaceindicator3d
  _normals2d
  _normals3d
  _distance2d
  _distance3d
)

foreach(obj ${PRECOMPILE})
  dune_add_pybind11_module(NAME ${obj})
  set_property(TARGET ${obj} PROPERTY LINK_LIBRARIES dunecommon dunegeometry dunegrid APPEND)
  if(SKBUILD)
    install(TARGETS ${obj} LIBRARY DESTINATION python/dune/mmesh)
  endif()
endforeach()
