#include <config.h>

#include <dune/mmesh/mmesh.hh>
#include <dune/python/grid/hierarchical.hh>
#include <dune/python/mmesh/utility.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

PYBIND11_MODULE( _normals3d, module )
{
  auto cls = Dune::Python::insertClass< Dune::Fem::Normals< typename Dune::MovingMesh<3>::InterfaceGrid::LeafGridView > >( module,
    "Normals",
    Dune::Python::GenerateTypeName("Dune::Fem::Normals<typename Dune::MovingMesh<3>::InterfaceGrid::LeafGridView>"),
    Dune::Python::IncludeFiles{"dune/mmesh/mmesh.hh", "dune/python/grid/hierarchical.hh", "dune/python/mmesh/utility.hh"}
  ).first;
  Dune::Fem::registerNormals( module, cls );
}