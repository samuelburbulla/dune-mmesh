#include <config.h>

#include <dune/mmesh/mmesh.hh>
#include <dune/python/grid/hierarchical.hh>
#include <dune/python/mmesh/grid.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

PYBIND11_MODULE( _interfacegrid3d, module )
{
  auto cls = Dune::Python::insertClass< typename Dune::MovingMesh<3>::InterfaceGrid >( module,
    "HierarchicalGrid",
    pybind11::dynamic_attr(),
    Dune::Python::GenerateTypeName("typename Dune::MovingMesh<3>::InterfaceGrid"),
    Dune::Python::IncludeFiles{"dune/mmesh/mmesh.hh", "dune/python/grid/hierarchical.hh"}
  ).first;
  Dune::Python::MMIFGrid::registerHierarchicalGrid( module, cls );
}
