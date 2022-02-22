#include <config.h>

#include <dune/mmesh/mmesh.hh>
#include <dune/python/grid/hierarchical.hh>
#include <dune/python/mmesh/distance.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

PYBIND11_MODULE( _distance3d, module )
{
  auto cls = Dune::Python::insertClass< Dune::Fem::Distance< typename Dune::MovingMesh<3>::LeafGridView > >( module,
    "Distance",
    Dune::Python::GenerateTypeName("Dune::Fem::Distance<typename Dune::MovingMesh<3>::LeafGridView>"),
    Dune::Python::IncludeFiles{"dune/mmesh/mmesh.hh", "dune/python/grid/hierarchical.hh", "dune/python/mmesh/distance.hh"}
  ).first;
  Dune::Fem::registerDistance( module, cls );
}
