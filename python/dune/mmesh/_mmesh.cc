#include <config.h>

#include <dune/mmesh/mmesh.hh>
#include <dune/python/grid/hierarchical.hh>
#include <dune/python/mmesh/grid.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

PYBIND11_MODULE( _mmesh, module )
{
  auto mmesh2 = Dune::Python::insertClass< Dune::MovingMesh<2> >( module,
    "mmesh2",
    pybind11::dynamic_attr(),
    Dune::Python::GenerateTypeName("Dune::MovingMesh<2>"),
    Dune::Python::IncludeFiles{"dune/mmesh/mmesh.hh", "dune/python/grid/hierarchical.hh"}
  ).first;
  Dune::Python::MMGrid::registerHierarchicalGrid( module, mmesh2 );


  auto mmesh3 = Dune::Python::insertClass< Dune::MovingMesh<3> >( module,
    "mmesh3",
    pybind11::dynamic_attr(),
    Dune::Python::GenerateTypeName("Dune::MovingMesh<3>"),
    Dune::Python::IncludeFiles{"dune/mmesh/mmesh.hh", "dune/python/grid/hierarchical.hh"}
  ).first;
  Dune::Python::MMGrid::registerHierarchicalGrid( module, mmesh3 );


  auto mimesh2 = Dune::Python::insertClass< typename Dune::MovingMesh<2>::InterfaceGrid >( module,
    "mimesh2",
    pybind11::dynamic_attr(),
    Dune::Python::GenerateTypeName("typename Dune::MovingMesh<2>::InterfaceGrid"),
    Dune::Python::IncludeFiles{"dune/mmesh/mmesh.hh", "dune/python/grid/hierarchical.hh"}
  ).first;
  Dune::Python::MMIFGrid::registerHierarchicalGrid( module, mimesh2 );


  auto mimesh3 = Dune::Python::insertClass< typename Dune::MovingMesh<3>::InterfaceGrid >( module,
    "mimesh3",
    pybind11::dynamic_attr(),
    Dune::Python::GenerateTypeName("typename Dune::MovingMesh<3>::InterfaceGrid"),
    Dune::Python::IncludeFiles{"dune/mmesh/mmesh.hh", "dune/python/grid/hierarchical.hh"}
  ).first;
  Dune::Python::MMIFGrid::registerHierarchicalGrid( module, mimesh3 );

}
