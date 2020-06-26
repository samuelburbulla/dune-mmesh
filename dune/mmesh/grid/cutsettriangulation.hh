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

    public:
      CutSetTriangulation(const Entity& first, const Entity& second)
      {
        static_assert( dim == 2 );

        const auto& geo1 = first.geometry();
        const auto& geo2 = second.geometry();

        std::array<GlobalCoordinate, 3> firstPoints, secondPoints;

        for ( int i = 0; i < 3; ++i )
        {
          firstPoints[i] = geo1.corner(i);
          secondPoints[i] = geo2.corner(i);
        }

        using PC = Dune::PolygonCutting<ctype, GlobalCoordinate>;
        auto points = PC::sutherlandHodgman(firstPoints, secondPoints);

        if (points.size() < 3)
          return;

        for (int i = 1; i < points.size()-1; ++i)
        {
          EntityImpl entity ( { points[0], points[i], points[i+1] } );
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
