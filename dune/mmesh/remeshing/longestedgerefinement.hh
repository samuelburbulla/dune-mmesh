// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*!
 * \file
 * \ingroup MMesh Adaptive
 * \brief   Class defining a longest edge refinement strategy.
 */

#ifndef DUNE_MMESH_REMESHING_LONGESTEDGEREFINEMENT_HH
#define DUNE_MMESH_REMESHING_LONGESTEDGEREFINEMENT_HH

#include <memory>
#include <dune/common/exceptions.hh>
#include <dune/grid/common/partitionset.hh>

namespace Dune
{

/*!
 * \ingroup MMesh Adaptive
 * \brief   Class defining a longest edge refinement strategy.
 */
template<class Grid>
class LongestEdgeRefinement
{
    static constexpr int dim = Grid::dimension;
    static constexpr int edgeCodim = dim - 1;
    static constexpr int vertexCodim = dim;
    using ctype = typename Grid::ctype;
    using Element = typename Grid::Traits::template Codim<0>::Entity;
    using Vertex = typename Grid::Traits::template Codim<vertexCodim>::Entity;

public:
    /*!
     * \brief Returns the refinement/coarsening point for each grid cell.
     */

    /*!
     * \brief return refinement point (center of longest edge)
     *
     * \param element A grid element
     */
    template<class Element>
    static auto refinement (const Element& element)
    {
      ctype longest = -1e100;
      int longestIdx = -1;

      for( std::size_t i = 0; i < element.subEntities(edgeCodim); ++i )
      {
        const ctype edgeLength = element.template subEntity<edgeCodim>(i).geometry().volume();
        if (edgeLength > longest)
        {
          longest = edgeLength;
          longestIdx = i;
        }
      }

      const auto& longestEdge = element.template subEntity<edgeCodim>(longestIdx);
      return std::make_pair( longestEdge, longestEdge.geometry().center() );
    }

    /*!
     * \brief return coarsening vertex (vertex of shortest edge)
     *
     * \param element A grid element
     */
    template<class Element>
    static Vertex coarsening (const Element& element)
    {
      // for interface any vertex is fine
      if constexpr ( Grid::dimension != Grid::dimensionworld )
      {
        for( std::size_t i = 0; i < element.subEntities(vertexCodim); ++i )
        {
          const Vertex& v = element.template subEntity<vertexCodim>(i);
          if ( isRemoveable( v ) )
            return v;
        }
        return Vertex();
      }

      // return some non-interface vertex at the shortest edge
      ctype shortest = 1e100;
      Vertex vertex;

      for( std::size_t i = 0; i < element.subEntities(edgeCodim); ++i )
      {
        const auto& edge = element.template subEntity<edgeCodim>(i);
        const ctype edgeLength = edge.geometry().volume();

        if (edgeLength < shortest)
        {
          const auto& v = edge.impl().template subEntity<vertexCodim>(0);
          if( !v.impl().isInterface() && !atBoundary(v) ) // give prio to interior points
          {
            vertex = v;
            shortest = edgeLength;
          }
          else
          {
            const auto& v2 = edge.impl().template subEntity<vertexCodim>(1);
            if( !v2.impl().isInterface() && boundaryFlag(v2) == 0 )
            {
              vertex = v2;
              shortest = edgeLength;
            }
          }
        }
      }

      return vertex;
    }

    //! return if vertex is incident to infinite vertex
    static inline bool atBoundary( const Vertex& v )
    {
      const auto& hostgrid = v.impl().grid().getHostGrid();
      for ( const auto& vertex : incidentVertices( v, true ) )
        if ( hostgrid.is_infinite( vertex.impl().hostEntity() ) )
          return true;
      return false;
    }

    //! return if vertex is removable at the boundary
    static inline int boundaryFlag( const Vertex& v )
    {
      // check if we have to initialize the boundary flag
      if ( v.impl().boundaryFlag() == -1 )
      {
        v.impl().hostEntity()->info().boundaryFlag = 0;
        if ( atBoundary( v ) )
        {
          // check that all incident vertices are in one plane
          std::vector< Vertex > incidentAtBoundary;
          for ( const auto& vertex : incidentVertices( v ) )
            if ( atBoundary( vertex ) )
              incidentAtBoundary.push_back( vertex );

          // could be that there are more than two incident boundary vertices, e.g. in grid corners
          if( incidentAtBoundary.size() > 2 )
            return 0;

          assert( incidentAtBoundary.size() == 2 );
          auto d1 = v.geometry().center();
          auto d2 = d1;
          d1 -= incidentAtBoundary[0].geometry().center();
          d2 -= incidentAtBoundary[1].geometry().center();
          d1 /= d1.two_norm();
          d2 /= d2.two_norm();
          if ( ( d1 + d2 ).two_norm() > 1e-8 && ( d1 - d2 ).two_norm() > 1e-8 )
            v.impl().hostEntity()->info().boundaryFlag = 1;
        }
      }

      return v.impl().boundaryFlag();
    }

    //! return if interface vertex is neither a tip nor a junction
    template< class Vertex >
    static inline bool isRemoveable ( const Vertex& vertex )
    {
      int count = 0;
      for ( const auto& edge : incidentInterfaceElements( vertex ) )
      {
        std::ignore = edge;
        count++;
      }

      return ( count == 2 );
    }
};

} // end namespace Dune

#endif
