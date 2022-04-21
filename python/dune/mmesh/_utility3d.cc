#include <config.h>

#if HAVE_DUNE_FEM

#include <dune/mmesh/mmesh.hh>
#include <dune/python/grid/hierarchical.hh>
#include <dune/python/mmesh/utility.hh>
#include <dune/python/mmesh/distance.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

PYBIND11_MODULE( _utility3d, module )
{
  // InterfaceIndicator
  auto clsIndicator = Dune::Python::insertClass< Dune::Fem::InterfaceIndicator< typename Dune::MovingMesh<3>::LeafGridView > >(
    module,
    "InterfaceIndicator",
    Dune::Python::GenerateTypeName("Dune::Fem::InterfaceIndicator<typename Dune::MovingMesh<3>::LeafGridView>"),
    Dune::Python::IncludeFiles{"dune/mmesh/mmesh.hh", "dune/python/grid/hierarchical.hh", "dune/python/mmesh/utility.hh"}
  ).first;
  Dune::Fem::registerInterfaceIndicator( module, clsIndicator );

  // Normals
  auto clsNormals = Dune::Python::insertClass< Dune::Fem::Normals< typename Dune::MovingMesh<3>::InterfaceGrid::LeafGridView > >(
    module,
    "Normals",
    Dune::Python::GenerateTypeName("Dune::Fem::Normals<typename Dune::MovingMesh<3>::InterfaceGrid::LeafGridView>"),
    Dune::Python::IncludeFiles{"dune/mmesh/mmesh.hh", "dune/python/grid/hierarchical.hh", "dune/python/mmesh/utility.hh"}
  ).first;
  Dune::Fem::registerNormals( module, clsNormals );

  // Distance
  auto clsDistance = Dune::Python::insertClass< Dune::Fem::Distance< typename Dune::MovingMesh<3>::LeafGridView > >(
    module,
    "Distance",
    Dune::Python::GenerateTypeName("Dune::Fem::Distance<typename Dune::MovingMesh<3>::LeafGridView>"),
    Dune::Python::IncludeFiles{"dune/mmesh/mmesh.hh", "dune/python/grid/hierarchical.hh", "dune/python/mmesh/distance.hh"}
  ).first;
  Dune::Fem::registerDistance( module, clsDistance );
}

#endif //HAVE_DUNE_FEM
