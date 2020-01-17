// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*!
 * \file
 * \ingroup MMesh Adaptive
 * \brief   Class defining a longest edge refinement strategy.
 */

#ifndef DUNE_MMESH_REMESHING_INTERFACEREFINEMENT_HH
#define DUNE_MMESH_REMESHING_INTERFACEREFINEMENT_HH

#include <memory>
#include <dune/common/exceptions.hh>
#include <dune/grid/common/partitionset.hh>

namespace Dune
{

/*!
 * \ingroup MMesh Adaptive
 * \brief   Class defining a interface refinement strategy.
 */
template<class Grid>
class InterfaceRefinement
{
    static constexpr int dim = Grid::dimensionworld;
    static constexpr int edgecd = dim - 1;
    static constexpr int vertexCodim = dim;
    using ctype = typename Grid::ctype;
    using Element = typename Grid::Traits::template Codim<0>::Entity;

public:
    /*!
     * \brief Returns the refinement/coarsening point for each grid cell.
     */

    /*!
     * \brief return refinement point (center of longest interface edge)
     *
     * \param element A grid element
     */
    template<class Element>
    static auto refinement (const Element& element)
    {
      ctype longest = -1e100;
      int longestIdx = -1;

      for( int i = 0; i < element.subEntities(edgecd); ++i )
      {
        const auto& edge = element.template subEntity<edgecd>(i);

        if( element.impl().grid().isInterface( edge ) )
        {
          ctype edgelength = edge.geometry().volume();
          if (edgelength > longest)
          {
            longest = edgelength;
            longestIdx = i;
          }
        }
      }
      assert( longestIdx != -1 );

      const auto& longestEdge = element.template subEntity<edgecd>(longestIdx);
      return std::make_pair( longestEdge, longestEdge.geometry().center() );
    }

    /*!
     * \brief return coarsening vertex (vertex of shortest edge with highest insertion index)
     *
     * \param element A grid element
     */
    template<class Element>
    static auto coarsening (const Element& element, const ctype minH)
    {
      ctype shortest = 1e100;
      int shortestIdx = -1;

      for( int i = 0; i < element.subEntities(edgecd); ++i )
      {
        const auto& edge = element.template subEntity<edgecd>(i);

        if( element.impl().grid().isInterface( edge ) )
        {
          const ctype edgelength = edge.geometry().volume();
          if (edgelength < shortest)
          {
            shortest = edgelength;
            shortestIdx = i;
          }
        }
      }

      if( shortestIdx == -1 || shortest >= minH )
      {
        for( int i = 0; i < element.subEntities(vertexCodim); ++i )
        {
          const auto& vertex = element.template subEntity<vertexCodim>(i);
          if( !vertex.impl().isInterface() )
            return vertex;
        }
        assert(false);
      }

      const auto& shortestEdge = element.template subEntity<edgecd>(shortestIdx);

      const auto& v0 = shortestEdge.impl().template subEntity<vertexCodim>(0);
      const auto& v1 = shortestEdge.impl().template subEntity<vertexCodim>(1);

      // try to return the vertex which is not interface
      if ( v0.impl().isInterface() ^ v1.impl().isInterface() )
        return v0.impl().isInterface() ? v1 : v0;

      // otherwise, return vertex with higher insertion level
      return ( v1.impl().insertionLevel() > v0.impl().insertionLevel() ) ? v1 : v0;
    }
};

} // end namespace Dumux

#endif
