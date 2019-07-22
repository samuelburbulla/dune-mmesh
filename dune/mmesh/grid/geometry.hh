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

// CGAL includes
#include <CGAL/Kernel/global_functions.h>

// local includes
#include <dune/mmesh/grid/pointfieldvector.hh>

namespace Dune
{

  /** \brief Geometry
   */

  template<int mydim, int coorddim, class GridImp>
  class MMeshGeometry {};

  //! 2D Geometry

  template<int mydim, class GridImp>
  class MMeshGeometry<mydim, 2, GridImp> :
    public AffineGeometry <typename GridImp::ctype, mydim, 2>
  {
    static constexpr int coorddim = 2;
    typedef AffineGeometry <typename GridImp::ctype, mydim, coorddim> BaseType;
    typedef FieldVector<typename GridImp::ctype, coorddim> FVector;

  public:
    //! Constructor from host geometry with codim 0
    MMeshGeometry(const typename GridImp::template HostGridEntity<0>& hostEntity)
     : BaseType( GeometryTypes::simplex(mydim),
         std::array<FVector, 3>( {
           makeFieldVector( hostEntity->vertex(0)->point() ),
           makeFieldVector( hostEntity->vertex(1)->point() ),
           makeFieldVector( hostEntity->vertex(2)->point() )
         } ) )
    {
        // obtain circumcenter by CGAL
        circumcenter_ = makeFieldVector(
            CGAL::circumcenter(
                hostEntity->vertex(0)->point(),
                hostEntity->vertex(1)->point(),
                hostEntity->vertex(2)->point()
            )
        );
    }

    //! Constructor from host geometry with codim 0
    MMeshGeometry(const typename GridImp::HostGridType::Face& face)
     : BaseType( GeometryTypes::simplex(mydim),
         std::array<FVector, 3>( {
           makeFieldVector( face.vertex(0)->point() ),
           makeFieldVector( face.vertex(1)->point() ),
           makeFieldVector( face.vertex(2)->point() )
         } ) )
    {
        // obtain circumcenter by CGAL
        circumcenter_ = makeFieldVector(
            CGAL::circumcenter(
                face.vertex(0)->point(),
                face.vertex(1)->point(),
                face.vertex(2)->point()
            )
        );
    }

    //! Constructor from host geometry with codim 1
    MMeshGeometry(const typename GridImp::template HostGridEntity<1>& hostEntity)
     : BaseType( GeometryTypes::simplex(mydim),
         std::array<FVector, 2>( {
           makeFieldVector( hostEntity.first->vertex( (hostEntity.second+1)%3 )->point() ),
           makeFieldVector( hostEntity.first->vertex( (hostEntity.second+2)%3 )->point() )
       } ) )
    {
        circumcenter_ = this->corner(0);
        circumcenter_ += this->corner(1);
        circumcenter_ *= 0.5;
    }

    //! Constructor from host geometry with codim 2
    MMeshGeometry(const typename GridImp::template HostGridEntity<2>& hostEntity)
     : BaseType( GeometryTypes::simplex(mydim), std::array<FVector, 1>( { makeFieldVector( hostEntity->point() ) } ) ),
       circumcenter_( this->corner(0) )
    {}

    /** \brief Obtain the circumcenter */
    const FVector circumcenter () const
    {
      return circumcenter_;
    }

  private:
    FVector circumcenter_;
  };

  //! 3D Geometry

  template<int mydim, class GridImp>
  class MMeshGeometry<mydim, 3, GridImp> :
    public AffineGeometry <typename GridImp::ctype, mydim, 3>
  {
    static constexpr int coorddim = 3;
    typedef AffineGeometry <typename GridImp::ctype, mydim, coorddim> BaseType;
    typedef typename GridImp::ctype ctype;
    typedef FieldVector<ctype, coorddim> FVector;

  public:
    //! Constructor from host geometry with codim 0
    MMeshGeometry(const typename GridImp::template HostGridEntity<0>& hostEntity)
     : BaseType( GeometryTypes::simplex(mydim),
         std::array<FVector, 4>( {
           makeFieldVector( hostEntity->vertex(0)->point() ),
           makeFieldVector( hostEntity->vertex(1)->point() ),
           makeFieldVector( hostEntity->vertex(2)->point() ),
           makeFieldVector( hostEntity->vertex(3)->point() )
         } ) )
    {
        // obtain circumcenter by CGAL
        circumcenter_ = makeFieldVector(
            CGAL::circumcenter(
                hostEntity->vertex(0)->point(),
                hostEntity->vertex(1)->point(),
                hostEntity->vertex(2)->point(),
                hostEntity->vertex(3)->point()
            )
        );
    }

    //! Constructor from host geometry with codim 0
    MMeshGeometry(const typename GridImp::HostGridType::Cell& cell)
     : BaseType( GeometryTypes::simplex(mydim),
         std::array<FVector, 4>( {
           makeFieldVector( cell.vertex(0)->point() ),
           makeFieldVector( cell.vertex(1)->point() ),
           makeFieldVector( cell.vertex(2)->point() ),
           makeFieldVector( cell.vertex(3)->point() )
         } ) )
    {
        // obtain circumcenter by CGAL
        circumcenter_ = makeFieldVector(
            CGAL::circumcenter(
                cell.vertex(0)->point(),
                cell.vertex(1)->point(),
                cell.vertex(2)->point(),
                cell.vertex(3)->point()
            )
        );
    }

    //! Constructor from host geometry with codim 1
    MMeshGeometry(const typename GridImp::template HostGridEntity<1>& hostEntity)
     : BaseType( GeometryTypes::simplex(mydim), getVertices<1>(hostEntity) )
    {
      // obtain circumcenter
      const auto& cell     = hostEntity.first;
      const auto& facetIdx = hostEntity.second;

      // use the CGAL index convention to obtain the vertices
      circumcenter_ = makeFieldVector(
        CGAL::circumcenter(
          cell->vertex( (facetIdx + 1) % 4 )->point(),
          cell->vertex( (facetIdx + 2) % 4 )->point(),
          cell->vertex( (facetIdx + 3) % 4 )->point()
        )
      );
    }

    //! Constructor from host geometry with codim 1
    MMeshGeometry(const typename GridImp::template HostGridEntity<2>& hostEntity)
     : BaseType( GeometryTypes::simplex(mydim), getVertices<2>( hostEntity ) )
    {
      // obtain circumcenter
      circumcenter_ = this->corner(0);
      circumcenter_ += this->corner(1);
      circumcenter_ *= 0.5;
    }

    //! Constructor from host geometry with codim 3
    MMeshGeometry(const typename GridImp::template HostGridEntity<3>& hostEntity)
     : BaseType( GeometryTypes::simplex(mydim), std::array<FVector, 1>( { makeFieldVector( hostEntity->point() ) } ) ),
                 circumcenter_( this->corner(0) )
    {}

    /** \brief Obtain the circumcenter */
    const FVector circumcenter () const
    {
      return circumcenter_;
    }

  private:

     template< int codim, typename Enable = std::enable_if_t< codim == 1 > >
     static inline std::array<FVector, 3> getVertices ( const typename GridImp::template HostGridEntity<1>& hostEntity )
     {
       std::array<FVector, 3> vertices;

       const auto& cell     = hostEntity.first;
       const auto& facetIdx = hostEntity.second;

       // use the CGAL index convention to obtain the vertices
       for ( int i = 0; i < 3; ++i )
         vertices[i] = makeFieldVector( cell->vertex( (facetIdx + i + 1) % 4 )->point() );

       return vertices;
     }

     template< int codim, typename Enable = std::enable_if_t< codim == 2 > >
     static inline std::array<FVector, 2> getVertices ( const typename GridImp::template HostGridEntity<2>& hostEntity )
     {
       std::array<FVector, 2> vertices;

       const auto& cell       = hostEntity.first;
       const auto& vertexIdx1 = hostEntity.second;
       const auto& vertexIdx2 = hostEntity.third;

       vertices[0] = makeFieldVector( cell->vertex( vertexIdx1 )->point() );
       vertices[1] = makeFieldVector( cell->vertex( vertexIdx2 )->point() );

       return vertices;
     }

    FVector circumcenter_;
  };

}  // namespace Dune

#endif
