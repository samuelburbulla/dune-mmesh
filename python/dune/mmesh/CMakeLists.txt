add_python_targets(mmesh
  __init__
  _grids
  _move
  _skeletontrace
  _solve
  _utility
)

set(PRECOMPILE
  _grid2d
  _grid3d
  _interfacegrid2d
  _interfacegrid3d
)

foreach(obj ${PRECOMPILE})
  dune_add_pybind11_module(NAME ${obj})
  set_property(TARGET ${obj} PROPERTY LINK_LIBRARIES dunecommon dunegeometry dunegrid dunemmesh APPEND)
  if(SKBUILD)
    install(TARGETS ${obj} LIBRARY DESTINATION python/dune/mmesh)
  endif()
endforeach()
