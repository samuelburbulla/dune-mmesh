 // -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_GRID_GMSHREADER_HH
#define DUNE_MMESH_GRID_GMSHREADER_HH

// System includes
#include <memory>
#include <utility>
#include <vector>

// Dune includes
#include <dune/grid/common/intersection.hh>
#include <dune/grid/io/file/gmshreader.hh>

// MMesh includes
#include <dune/mmesh/grid/gmshgridfactory.hh>

namespace Dune
{

  // GmshReader for MMesh
  // --------------------
  template<int dim>
  struct GmshReader< Dune::MovingMesh< dim > >
  {
    using Grid = Dune::MovingMesh< dim >;

    static std::shared_ptr<Grid> read (const std::string& fileName, bool verbose = true, bool insertBoundarySegments = true)
    {
      Dune::GmshGridFactory<Grid> gmshFactory( fileName );
      return gmshFactory.grid();
    }

    static void read (Dune::GridFactory<Grid>& factory, const std::string& fileName,
                      bool verbose = true, bool insertBoundarySegments = true)
    {
      Dune::GmshGridFactory<Grid> gmshFactory( factory, fileName );
    }

    static void read (Dune::GridFactory<Grid>& factory, const std::string& fileName, std::vector<int>& boundaryIDs, std::vector<int>& elementsIDs, bool verbose = false, bool boundarySegments = false)
    {
      Dune::GmshGridFactory<Grid> gmshFactory( factory, fileName );
    }

    std::vector<int> elementsIDs;
  };

} // end namespace Dune

#endif
