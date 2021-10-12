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

dune_add_pybind11_module(NAME _mmesh)
set_property(TARGET _mmesh PROPERTY LINK_LIBRARIES dunecommon dunegeometry dunegrid dunemmesh APPEND)
install(TARGETS _mmesh LIBRARY DESTINATION python/dune/mmesh)
