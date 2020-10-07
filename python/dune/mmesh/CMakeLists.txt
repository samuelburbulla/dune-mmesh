add_python_targets(mmesh
  __init__
  _grids
  _skeletontrace
  _solve
  _utility
)

dune_add_pybind11_module(NAME _mmesh)
