// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_INTERFACE_GEOMETRY_HH
#define DUNE_MMESH_INTERFACE_GEOMETRY_HH

/** \file
 * \brief The MMeshInterfaceGridGeometry class and its specializations
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
  class MMeshInterfaceGridGeometry {};

  //! 2D Geometry

  template<int md, class GridImp>
  class MMeshInterfaceGridGeometry<md, 2, GridImp> :
    public AffineGeometry <typename GridImp::ctype, md, 2>
  {
  public:
    static constexpr int mydim = md;
    static constexpr int coorddim = 2;
    typedef AffineGeometry <typename GridImp::ctype, mydim, coorddim> BaseType;
    typedef FieldVector<typename GridImp::ctype, coorddim> FVector;

    enum { dimension = GridImp::dimension };
    enum { dimensionworld = GridImp::dimensionworld };
    enum { coorddimension = coorddim };
    enum { mydimension = mydim };

    //! Constructor from host geometry with codim 0
    MMeshInterfaceGridGeometry(const typename GridImp::template MMeshInterfaceEntity<0>& hostEntity)
     : BaseType( GeometryTypes::simplex(mydim), getVertices(hostEntity) )
    {}

    //! Constructor from vertex array
    MMeshInterfaceGridGeometry(const std::array<FVector, 2>& points)
     : BaseType( GeometryTypes::simplex(mydim), points )
    {}

    //! Constructor from host geometry with codim 1
    MMeshInterfaceGridGeometry(const typename GridImp::template MMeshInterfaceEntity<1>& hostEntity)
     : BaseType( GeometryTypes::simplex(mydim), std::array<FVector, 1>( { makeFieldVector( hostEntity->point() ) } ) )
    {}

    //! Constructor for local geometry of intersection from intersection index for 3D
    MMeshInterfaceGridGeometry(int i)
     : BaseType( GeometryTypes::simplex(1),
                 std::array<FVector, 2>(
                  {
                    FVector ( { i==2 ? 1.0 : 0.0, 0.0 } ),
                    FVector ( { i==0 ? 1.0 : 0.0, i>0 ? 1.0 : 0.0 } )
                  } ) )
    {}

    /** \brief Obtain the circumcenter */
    const FVector circumcenter () const
    {
      return computeCircumcenter(*this);
    }

  private:
    static inline std::array<FVector, 2> getVertices ( const typename GridImp::template MMeshInterfaceEntity<0>& hostEntity )
    {
      const auto& cgalIdx = MMeshInterfaceImpl::computeCGALIndices<typename GridImp::template MMeshInterfaceEntity<0>, 1> ( hostEntity );

      std::array<FVector, 2> vertices ( {
        makeFieldVector( hostEntity.first->vertex( cgalIdx[0] )->point() ),
        makeFieldVector( hostEntity.first->vertex( cgalIdx[1] )->point() )
      } );

      return vertices;
    }
  };

  //! 3D Geometry

  template<int md, class GridImp>
  class MMeshInterfaceGridGeometry<md, 3, GridImp> :
    public AffineGeometry <typename GridImp::ctype, md, 3>
  {
  public:
    static constexpr int mydim = md;
    static constexpr int coorddim = 3;
    typedef AffineGeometry <typename GridImp::ctype, mydim, coorddim> BaseType;
    typedef typename GridImp::ctype ctype;
    typedef FieldVector<ctype, coorddim> FVector;

    enum { dimension = GridImp::dimension };
    enum { dimensionworld = GridImp::dimensionworld };
    enum { coorddimension = coorddim };
    enum { mydimension = mydim };

    //! Constructor from host geometry with codim 0
    MMeshInterfaceGridGeometry(const typename GridImp::template MMeshInterfaceEntity<0>& hostEntity)
     : BaseType( GeometryTypes::simplex(mydim), getVertices<0>(hostEntity) )
    {}

    //! Constructor from vertex list
    MMeshInterfaceGridGeometry(const std::array<FVector, 3>& points)
     : BaseType( GeometryTypes::simplex(mydim), points )
    {}

    //! Constructor from host geometry with codim 1
    MMeshInterfaceGridGeometry(const typename GridImp::template MMeshInterfaceEntity<1>& hostEntity)
     : BaseType( GeometryTypes::simplex(mydim), getVertices<1>( hostEntity ) )
    {}

    //! Constructor from host geometry with codim 2
    MMeshInterfaceGridGeometry(const typename GridImp::template MMeshInterfaceEntity<2>& hostEntity)
     : BaseType( GeometryTypes::simplex(mydim), std::array<FVector, 1>( { makeFieldVector( hostEntity->point() ) } ) )
    {}

    /** \brief Obtain the circumcenter */
    const FVector circumcenter () const
    {
      return computeCircumcenter(*this);
    }

  private:

    template< int codim, typename Enable = std::enable_if_t< codim == 0 > >
    static inline std::array<FVector, 3> getVertices ( const typename GridImp::template MMeshInterfaceEntity<0>& hostEntity )
    {
      const auto& cgalIdx = MMeshInterfaceImpl::computeCGALIndices<typename GridImp::template MMeshInterfaceEntity<0>, 2> ( hostEntity );

      std::array<FVector, 3> vertices;

      const auto& cell = hostEntity.first;

      // use the CGAL index convention to obtain the vertices
      for ( int i = 0; i < 3; ++i )
       vertices[i] = makeFieldVector( cell->vertex( cgalIdx[i] )->point() );

      return vertices;
    }

    template< int codim, typename Enable = std::enable_if_t< codim == 1 > >
    static inline std::array<FVector, 2> getVertices ( const typename GridImp::template MMeshInterfaceEntity<1>& hostEntity )
    {
      const auto& cell = hostEntity.first;
      const auto& i = hostEntity.second;
      const auto& j = hostEntity.third;

      std::array<FVector, 2> vertices;
      vertices[0] = makeFieldVector( cell->vertex( i )->point() );
      vertices[1] = makeFieldVector( cell->vertex( j )->point() );

      if ( cell->vertex( i )->info().id > cell->vertex( j )->info().id )
        std::swap( vertices[0], vertices[1] );

      return vertices;
    }
  };

  //! The local geometry (2D)
  template<int md, class GridImp>
  class MMeshInterfaceGridGeometry<md, 1, GridImp> :
    public AffineGeometry <typename GridImp::ctype, md, 1>
  {
  public:
    static constexpr int mydim = md;
    static constexpr int coorddim = 1;
    typedef AffineGeometry <typename GridImp::ctype, mydim, coorddim> BaseType;
    typedef FieldVector<typename GridImp::ctype, coorddim> FVector;

    enum { dimension = GridImp::dimension };
    enum { dimensionworld = GridImp::dimensionworld };
    enum { coorddimension = coorddim };
    enum { mydimension = mydim };

    //! Constructor for local geometry of intersection from intersection index
    MMeshInterfaceGridGeometry(int i)
     : BaseType( GeometryTypes::simplex(mydim), std::array<FVector, 1>( { i } ) )
    {}

    //! Constructor from vertex array
    MMeshInterfaceGridGeometry(const std::array<FVector, 2>& points)
     : BaseType( GeometryTypes::simplex(mydim), points )
    {}

  };

}  // namespace Dune

#endif
