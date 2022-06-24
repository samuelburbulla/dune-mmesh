#include <config.h>

#include <dune/mmesh/mmesh.hh>
#include <dune/python/grid/hierarchical.hh>
#include <dune/python/mmesh/grid.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

PYBIND11_MODULE( _grid3d, module )
{
  auto cls = Dune::Python::insertClass< Dune::MovingMesh<3> >( module,
    "HierarchicalGrid",
    pybind11::dynamic_attr(),
    Dune::Python::GenerateTypeName("Dune::MovingMesh<3>"),
    Dune::Python::IncludeFiles{"dune/mmesh/mmesh.hh", "dune/python/grid/hierarchical.hh"}
  ).first;
  Dune::Python::MMGrid::registerHierarchicalGrid( module, cls );
}
