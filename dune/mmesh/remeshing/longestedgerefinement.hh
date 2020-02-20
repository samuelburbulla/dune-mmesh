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
      if ( Grid::dimension != Grid::dimensionworld )
        return element.template subEntity<vertexCodim>(0);

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
          if( !v.impl().isInterface() && !atBoundary(v) )
          {
            vertex = v;
            shortest = edgeLength;
          }
          else
          {
            const auto& v2 = edge.impl().template subEntity<vertexCodim>(1);
            if( !v2.impl().isInterface() && !atBoundary(v2) )
            {
              vertex = v2;
              shortest = edgeLength;
            }
          }
        }
      }

      return vertex;
    }

  private:
    //! return if vertex is at the boundary
    // TODO improve this by marking vertices
    static inline bool atBoundary( const Vertex& v )
    {
      const auto& hostgrid = v.impl().grid().getHostGrid();
      for ( const auto& vertex : incidentVertices( v, true ) )
        if ( hostgrid.is_infinite( vertex.impl().hostEntity() ) )
          return true;
      return false;
    }
};

} // end namespace Dumux

#endif
