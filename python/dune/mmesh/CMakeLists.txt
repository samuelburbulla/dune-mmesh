add_python_targets(mmesh
  __init__
  __main__
  _grids
  _move
  _skeletontrace
  _solve
  _utility
)

dune_symlink_to_source_files(FILES
  tutorial
)

set(PRECOMPILE
  _mmesh2
  _mmesh3
  _mimesh2
  _mimesh3
)

foreach(obj ${PRECOMPILE})
  dune_add_pybind11_module(NAME ${obj})
  set_property(TARGET ${obj} PROPERTY LINK_LIBRARIES dunecommon dunegeometry dunegrid dunemmesh APPEND)
  install(TARGETS ${obj} LIBRARY DESTINATION python/dune/mmesh)
endforeach()
