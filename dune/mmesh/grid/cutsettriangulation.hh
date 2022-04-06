// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_GRID_CUTSETTRIANGULATION_HH
#define DUNE_MMESH_GRID_CUTSETTRIANGULATION_HH

/** \file
 * \brief The CutSetTriangulation class
 * This class computes the overlapping polytope of to intersecting triangles
 * and returns a retriangulation of this domain as temporary entities.
 */

#include <dune/mmesh/grid/polygoncutting.hh>
#include <dune/mmesh/grid/pointfieldvector.hh>

namespace Dune
{

  namespace MMeshImpl
  {

    template<class Entity>
    class CutSetTriangulation
    {
      static constexpr int dim = Entity::dimension;
      using ctype = typename Entity::Geometry::ctype;
      using GlobalCoordinate = typename Entity::Geometry::GlobalCoordinate;

      using EntityImpl = typename Entity::Implementation;
      using EntityList = std::vector< Entity >;
      using CachingEntity = MMeshCachingEntity< 0, dim, const typename EntityImpl::Grid >;

    public:
      CutSetTriangulation(const CachingEntity& caching, const Entity& element)
      {
        static_assert( dim == 2 );

        const auto& cgeo = caching.geometry();
        const auto& host = element.impl().hostEntity();

        std::array<GlobalCoordinate, 3> c, e;

        for ( int i = 0; i < 3; ++i )
        {
          c[i] = cgeo.corner(i);
          // obtain vertices in ccw order
          e[i] = makeFieldVector( host->vertex(i)->point() );
        }

        // check orientation of caching entity
        int o = (c[1][1]-c[0][1])*(c[2][0]-c[1][0])-(c[1][0]-c[0][0])*(c[2][1]-c[1][1]);
        if (o > 0) // clock wise
          std::swap(c[1], c[2]);

        using PC = Dune::PolygonCutting<ctype, GlobalCoordinate>;
        auto points = PC::sutherlandHodgman(c, e);

        if (points.size() < 3)
          return;

        // we know the intersection polygon of two triangles is convex
        for (std::size_t i = 1; i < points.size()-1; ++i)
        {
          EntityImpl entity ( &element.impl().grid(), { points[0], points[i], points[i+1] } );
          if ( entity.geometry().volume() > 1e-8 )
            triangles_.emplace_back( entity );
        }
      }

      const EntityList& triangles() const
      {
        return triangles_;
      }

      EntityList& triangles()
      {
        return triangles_;
      }

    private:
      EntityList triangles_;

    }; // end of CutSetTriangulation

  } // namespace MMeshImpl

} // namespace Dune

#endif
