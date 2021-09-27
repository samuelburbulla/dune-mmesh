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
