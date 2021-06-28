 // -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_INTERFACE_DGFPARSER_HH
#define DUNE_MMESH_INTERFACE_DGFPARSER_HH

// System includes
#include <memory>
#include <utility>
#include <vector>

// Dune includes
#include <dune/grid/common/intersection.hh>
#include <dune/grid/io/file/dgfparser/blocks/projection.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/parser.hh>

// MMesh includes
#include <dune/mmesh/grid/explicitgridfactory.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------
  template< class GridImp, class IntersectionImp >
  class Intersection;

  // dummy specialization of DGFGridFactory for MMeshInterfaceGrid
  // -------------------------------------------------------------
  template< class MMeshImp >
  struct DGFGridFactory< MMeshInterfaceGrid<MMeshImp> >
  {
    typedef MMeshInterfaceGrid<MMeshImp> Grid;
    const static int dimension = Grid::dimension;
    typedef MPIHelper::MPICommunicator MPICommunicatorType;
    typedef typename Grid::template Codim<0>::Entity Element;
    typedef typename Grid::template Codim<dimension>::Entity Vertex;

    explicit DGFGridFactory ( std::istream &input,
                              MPICommunicatorType comm = MPIHelper::getCommunicator() );
    explicit DGFGridFactory ( const std::string &filename,
                              MPICommunicatorType comm = MPIHelper::getCommunicator() );

    //! return grid pointer
    Grid* grid() const
    {
      DUNE_THROW(NotImplemented, "DGFGridFactory for MMeshInterfaceGrid");
      return grid_;
    }

    //! Returns if intersection was inserted
    template< class Intersection >
    bool wasInserted ( const Intersection &intersection ) const
    {
      DUNE_THROW(NotImplemented, "DGFGridFactory for MMeshInterfaceGrid");
      return false;
    }

    //! Return boundary id of intersection
    template< class Intersection >
    int boundaryId ( const Intersection &intersection ) const
    {
      DUNE_THROW(NotImplemented, "DGFGridFactory for MMeshInterfaceGrid");
      return 0;
    }

    //! Returns dgf element parameters for given element
    std::vector< double > parameter ( const Element &element )
    {
      DUNE_THROW(NotImplemented, "DGFGridFactory for MMeshInterfaceGrid");
      return {0.0};
    }

    //! Returns dgf vertex parameters for given vertex
    std::vector< double > parameter ( const Vertex &vertex )
    {
      DUNE_THROW(NotImplemented, "DGFGridFactory for MMeshInterfaceGrid");
      return {0.0};
    }

    //! Returns true if boundary parameters found
    bool haveBoundaryParameters () const
    {
      DUNE_THROW(NotImplemented, "DGFGridFactory for MMeshInterfaceGrid");
      return false;
    }

    template < class GG, class II >
    const DGFBoundaryParameter::type
    boundaryParameter ( const Intersection< GG, II > & intersection ) const
    {
      DUNE_THROW(NotImplemented, "DGFGridFactory for MMeshInterfaceGrid");
      return DGFBoundaryParameter::type();
    }

    template< int codim >
    int numParameters () const
    {
      DUNE_THROW(NotImplemented, "DGFGridFactory for MMeshInterfaceGrid");
      return 0;
    }

  private:
    Grid* grid_;
    DuneGridFormatParser dgf_;
  };


  // Implementation of DGFGridFactory for MMeshInterfaceGrid
  // -------------------------------------------------------

  //! DGFGridFactory for MMeshInterfaceGrid using std::istream
  template< class MMeshImp >
  inline DGFGridFactory< MMeshInterfaceGrid<MMeshImp> >
  ::DGFGridFactory ( std::istream &input, MPICommunicatorType comm )
    : dgf_( 0, 1 )
  {
    DUNE_THROW(NotImplemented, "DGFGridFactory for MMeshInterfaceGrid");
  }

  //! DGFGridFactory for MMeshInterfaceGrid using filename
  template< class MMeshImp >
  inline DGFGridFactory< MMeshInterfaceGrid<MMeshImp> >
  ::DGFGridFactory ( const std::string &filename, MPICommunicatorType comm )
    : dgf_( 0, 1 )
  {
    DUNE_THROW(NotImplemented, "DGFGridFactory for MMeshInterfaceGrid");
  }

}  // namespace Dune

#endif
