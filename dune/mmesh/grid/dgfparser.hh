 // -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_GRID_DGFPARSER_HH
#define DUNE_MMESH_GRID_DGFPARSER_HH

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

  // DGFGridFactory for MMesh
  // ------------------------
  template< class HostGrid, int dim >
  struct DGFGridFactory< MMesh<HostGrid, dim> >
  {
    typedef MMesh<HostGrid, dim> Grid;
    const static int dimension = Grid::dimension;
    typedef MPIHelper::MPICommunicator MPICommunicatorType;
    typedef typename Grid::template Codim<0>::Entity Element;
    typedef typename Grid::template Codim<dimension>::Entity Vertex;
    typedef Dune::MMeshExplicitGridFactory<Grid> GridFactory;

    explicit DGFGridFactory ( std::istream &input,
                              MPICommunicatorType comm = MPIHelper::getCommunicator() );
    explicit DGFGridFactory ( const std::string &filename,
                              MPICommunicatorType comm = MPIHelper::getCommunicator() );

    //! return grid pointer
    Grid* grid() const
    {
      return grid_;
    }

    //! Returns if intersection was inserted
    template< class Intersection >
    bool wasInserted ( const Intersection &intersection ) const
    {
      return factory_.wasInserted( intersection );
    }

    //! Return boundary id of intersection
    template< class Intersection >
    int boundaryId ( const Intersection &intersection ) const
    {
      return intersection.impl().boundaryId();
    }

    //! Returns dgf element parameters for given element
    std::vector< double > &parameter ( const Element &element )
    {
       if( numParameters< 0 >() <= 0 )
      {
        DUNE_THROW( InvalidStateException,
                    "Calling DGFGridFactory::parameter is only allowed if there are parameters." );
      }
      return dgf_.elParams[ factory_.insertionIndex( element ) ];
    }

    //! Returns dgf vertex parameters for given vertex
    std::vector< double > &parameter ( const Vertex &vertex )
    {
      if( numParameters< dimension >() <= 0 )
      {
        DUNE_THROW( InvalidStateException,
                    "Calling DGFGridFactory::parameter is only allowed if there are parameters." );
      }
      return dgf_.vtxParams[ factory_.insertionIndex( vertex ) ];
    }

    //! Returns true if boundary parameters found
    bool haveBoundaryParameters () const
    {
      return dgf_.haveBndParameters;
    }

    template < class GG, class II >
    const DGFBoundaryParameter::type
    boundaryParameter ( const Intersection< GG, II > & intersection ) const { return DGFBoundaryParameter::type(); }

    template< int codim >
    int numParameters () const { return 0; }

  private:
    bool generate( std::istream &input )
    {
      dgf_.element = DuneGridFormatParser::Simplex;
      dgf_.dimgrid = dimension;
      dgf_.dimw = dimension;

      if( !dgf_.readDuneGrid( input, dimension, dimension ) )
        return false;

      // Insert vertices
      for( int n = 0; n < dgf_.nofvtx; ++n )
      {
        typename GridFactory::WorldVector coord;
        for( std::size_t i = 0; i < dimension; ++i )
          coord[ i ] = dgf_.vtx[ n ][ i ];
        factory_.insertVertex( coord );
      }

      // Insert elements
      for( int n = 0; n < dgf_.nofelements; ++n )
        factory_.insertElement( GeometryTypes::simplex(dimension), dgf_.elements[ n ] );

      // Insert boundary segments and ids
      for( const auto& face : dgf_.facemap )
      {
        const auto& entityKey = face.first;
        const std::size_t boundaryId = face.second.first;

        std::vector< unsigned int > vertices;
        for( int i = 0; i < entityKey.size(); ++i )
          vertices.push_back( entityKey[i] );

        std::sort( vertices.begin(), vertices.end() );

        // insert boundary segment
        factory_.insertBoundarySegment( vertices );

        // insert boundary id
        std::size_t index = factory_.boundarySegments().size()-1;
        factory_.addBoundaryId( index, boundaryId );
      }

      grid_ = factory_.createGrid().release();
      return true;
    }

    Grid* grid_;
    GridFactory factory_;
    DuneGridFormatParser dgf_;
  };


  // DGFGridInfo for MMesh
  // ---------------------

  template< class HostGrid, int dim >
  struct DGFGridInfo< MMesh<HostGrid, dim> >
  {
    static int refineStepsForHalf ()
    {
      return 2;
    }

    static double refineWeight ()
    {
      return -1;
    }
  };


  // Implementation of DGFGridFactory for MMesh
  // ------------------------------------------

  //! DGFGridFactory for MMesh using std::istream
  template< class HostGrid, int dim >
  inline DGFGridFactory< MMesh<HostGrid, dim> >
  ::DGFGridFactory ( std::istream &input, MPICommunicatorType comm )
    : dgf_( 0, 1 )
  {
    input.clear();
    input.seekg( 0 );
    if( !input )
      DUNE_THROW(DGFException, "Error resetting input stream." );
    generate( input );
  }

  //! DGFGridFactory for MMesh using filename
  template< class HostGrid, int dim >
  inline DGFGridFactory< MMesh<HostGrid, dim> >
  ::DGFGridFactory ( const std::string &filename, MPICommunicatorType comm )
    : dgf_( 0, 1 )
  {
    std::ifstream input( filename.c_str() );
    if( !input )
      DUNE_THROW( DGFException, "Macrofile " << filename << " not found." );
    if( !generate( input ) )
      DUNE_THROW( DGFException, "Could not generate MMesh from macrofile " << filename );
    input.close();
  }

}  // namespace Dune

#endif
