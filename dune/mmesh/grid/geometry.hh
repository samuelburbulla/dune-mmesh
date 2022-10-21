// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_GRID_GEOMETRY_HH
#define DUNE_MMESH_GRID_GEOMETRY_HH

/** \file
 * \brief The MMeshGeometry class and its specializations
 * Inherits from Dune::AffineGeometry.
 */

// Dune includes
#include <dune/common/fmatrix.hh>
#include <dune/common/typetraits.hh>
#include <dune/grid/common/geometry.hh>
#include <dune/geometry/affinegeometry.hh>

// local includes
#include <dune/mmesh/grid/pointfieldvector.hh>

namespace Dune
{

  /** \brief Geometry
   */

  template<int mydim, int coorddim, class GridImp>
  class MMeshGeometry {};

  //! 2D Geometry

  template<int md, class GridImp>
  class MMeshGeometry<md, 2, GridImp> :
    public AffineGeometry <typename GridImp::ctype, md, 2>
  {
  public:
    static constexpr int mydim = md;
    static constexpr int coorddim = 2;
    typedef AffineGeometry <typename GridImp::ctype, mydim, coorddim> BaseType;
    typedef FieldVector<typename GridImp::ctype, coorddim> FVector;
    typedef typename GridImp::LeafIntersection Intersection;

    enum { dimension = GridImp::dimension };
    enum { dimensionworld = GridImp::dimensionworld };
    enum { coorddimension = coorddim };
    enum { mydimension = mydim };

    explicit MMeshGeometry()
     : BaseType(GeometryTypes::simplex(mydim),
       std::array<FVector, 3>( {
         FVector( { 0.0, 0.0 } ),
         FVector( { 1.0, 0.0 } ),
         FVector( { 0.0, 1.0 } )
       } ) ) {}

    //! Constructor from host geometry with codim 0
    MMeshGeometry(const typename GridImp::template HostGridEntity<0>& hostEntity)
     : BaseType( GeometryTypes::simplex(mydim), getVertices(hostEntity) )
    {}

    MMeshGeometry(const std::array<FVector, 3> points)
     : BaseType( GeometryTypes::simplex(mydim), points )
    {}

    //! Constructor from host geometry with codim 1
    MMeshGeometry(const typename GridImp::template HostGridEntity<1>& hostEntity)
     : BaseType( GeometryTypes::simplex(mydim), getVertices(hostEntity) )
    {}

    //! Constructor of local intersection geometry
    MMeshGeometry( int idx )
     : BaseType( GeometryTypes::simplex(mydim), getLocalVertices_( idx ) )
    {}

    //! Constructor from host geometry with codim 2
    MMeshGeometry(const typename GridImp::template HostGridEntity<2>& hostEntity)
     : BaseType( GeometryTypes::simplex(mydim), std::array<FVector, 1>( { makeFieldVector( hostEntity->point() ) } ) )
    {}

    /** \brief Obtain the circumcenter */
    const FVector circumcenter () const
    {
      return computeCircumcenter(*this);
    }

  private:
    static inline std::array<FVector, 3> getVertices ( const typename GridImp::template HostGridEntity<0> hostEntity )
    {
      const auto& cgalIdx = hostEntity->info().cgalIndex;

      std::array<FVector, 3> vertices ( {
        makeFieldVector( hostEntity->vertex(cgalIdx[0])->point() ),
        makeFieldVector( hostEntity->vertex(cgalIdx[1])->point() ),
        makeFieldVector( hostEntity->vertex(cgalIdx[2])->point() )
      } );

      return vertices;
    }

    static inline std::array<FVector, 2> getVertices ( const typename GridImp::template HostGridEntity<1>& hostEntity )
    {
      const auto& cgalIdx = hostEntity.first->info().cgalIndex;

      auto facetIdx = MMeshImpl::cgalFacetToDuneFacet<2, typename GridImp::template HostGridEntity<1>>( hostEntity );

      std::array<FVector, 2> vertices;
      for ( int k = 0; k < 2; ++k )
        vertices[k] = makeFieldVector( hostEntity.first->vertex( cgalIdx[ MMeshImpl::ref<2>().subEntity(facetIdx, 1, k, 2) ] )->point() );

      return vertices;
    }

    static inline std::array<FVector, 2> getLocalVertices_ ( int k )
    {
      static const std::array<FVector, 3> local = {
        FVector( { 0.0, 0.0 } ),
        FVector( { 1.0, 0.0 } ),
        FVector( { 0.0, 1.0 } )
      };

      return { local[ k==2 ? 1 : 0 ], local[ k==0 ? 1 : 2 ] };
    }
  };

  //! 3D Geometry

  template<int md, class GridImp>
  class MMeshGeometry<md, 3, GridImp> :
    public AffineGeometry <typename GridImp::ctype, md, 3>
  {
  public:
    static constexpr int mydim = md;
    static constexpr int coorddim = 3;
    typedef AffineGeometry <typename GridImp::ctype, mydim, coorddim> BaseType;
    typedef typename GridImp::ctype ctype;
    typedef FieldVector<ctype, coorddim> FVector;
    typedef typename GridImp::LeafIntersection Intersection;

    enum { dimension = GridImp::dimension };
    enum { dimensionworld = GridImp::dimensionworld };
    enum { coorddimension = coorddim };
    enum { mydimension = mydim };

    explicit MMeshGeometry()
     : BaseType(GeometryTypes::simplex(mydim),
       std::array<FVector, 4>( {
         FVector( { 0.0, 0.0, 0.0 } ),
         FVector( { 1.0, 0.0, 0.0 } ),
         FVector( { 0.0, 1.0, 0.0 } ),
         FVector( { 0.0, 0.0, 1.0 } )
       } ) ) {}

    //! Constructor from host geometry with codim 0
    MMeshGeometry(const typename GridImp::template HostGridEntity<0>& hostEntity)
     : BaseType( GeometryTypes::simplex(mydim), getVertices<0>(hostEntity) )
    {}

    MMeshGeometry(const std::array<FVector, 4> points)
     : BaseType( GeometryTypes::simplex(mydim), points )
    {}

    //! Constructor from host geometry with codim 1
    MMeshGeometry(const typename GridImp::template HostGridEntity<1>& hostEntity)
     : BaseType( GeometryTypes::simplex(mydim), getVertices<1>(hostEntity) )
    {}

    //! Constructor of local intersection geometry
    MMeshGeometry( int idx )
     : BaseType( GeometryTypes::simplex(mydim), getLocalVertices_( idx ) )
    {}

    //! Constructor from host geometry with codim 1
    MMeshGeometry(const typename GridImp::template HostGridEntity<2>& hostEntity)
     : BaseType( GeometryTypes::simplex(mydim), getVertices<2>( hostEntity ) )
    {}

    //! Constructor from host geometry with codim 3
    MMeshGeometry(const typename GridImp::template HostGridEntity<3>& hostEntity)
     : BaseType( GeometryTypes::simplex(mydim), std::array<FVector, 1>( { makeFieldVector( hostEntity->point() ) } ) )
    {}

    /** \brief Obtain the circumcenter */
    const FVector circumcenter () const
    {
      return computeCircumcenter(*this);
    }

  private:

      template< int codim, typename Enable = std::enable_if_t< codim == 0 > >
      static inline std::array<FVector, 4> getVertices ( const typename GridImp::template HostGridEntity<0>& hostEntity )
      {
        const auto& cgalIdx = hostEntity->info().cgalIndex;

        std::array<FVector, 4> vertices;
        for ( int i = 0; i < 4; ++i )
          vertices[i] = makeFieldVector( hostEntity->vertex( cgalIdx[i] )->point() );
        return vertices;
      }

      template< int codim, typename Enable = std::enable_if_t< codim == 1 > >
      static inline std::array<FVector, 3> getVertices ( const typename GridImp::template HostGridEntity<1>& hostEntity )
      {
        const auto& cgalIdx = hostEntity.first->info().cgalIndex;
        const auto& cell = hostEntity.first;

        auto facetIdx = MMeshImpl::cgalFacetToDuneFacet<3, typename GridImp::template HostGridEntity<1>>( hostEntity );

        std::array<FVector, 3> vertices;
        for ( int i = 0; i < 3; ++i )
          vertices[i] = makeFieldVector( cell->vertex( cgalIdx[ MMeshImpl::ref<3>().subEntity(facetIdx, 1, i, 3) ] )->point() );
       return vertices;
      }

      template< int codim, typename Enable = std::enable_if_t< codim == 2 > >
      static inline std::array<FVector, 2> getVertices ( const typename GridImp::template HostGridEntity<2>& hostEntity )
      {
        const auto& cgalIdx = hostEntity.first->info().cgalIndex;
        const auto& cell = hostEntity.first;

        auto edgeIdx = MMeshImpl::cgalEdgeToDuneEdge<3, typename GridImp::template HostGridEntity<2>>( hostEntity );

        std::array<FVector, 2> vertices;
        vertices[0] = makeFieldVector( cell->vertex( cgalIdx[ MMeshImpl::ref<3>().subEntity(edgeIdx, 2, 0, 3) ] )->point() );
        vertices[1] = makeFieldVector( cell->vertex( cgalIdx[ MMeshImpl::ref<3>().subEntity(edgeIdx, 2, 1, 3) ] )->point() );

        return vertices;
      }

      static inline std::array<FVector, 3> getLocalVertices_ ( int k )
      {
        static const std::array<FVector, 4> local = {
          FVector( { 0.0, 0.0, 0.0 } ),
          FVector( { 1.0, 0.0, 0.0 } ),
          FVector( { 0.0, 1.0, 0.0 } ),
          FVector( { 0.0, 0.0, 1.0 } )
        };

        return { local[ k<=2 ? 0 : 1 ], local[ k<=1 ? 1 : 2 ], local[ k==0 ? 2 : 3 ] };
      }
  };

}  // namespace Dune

#endif
