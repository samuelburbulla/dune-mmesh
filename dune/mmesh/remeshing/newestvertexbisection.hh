// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*!
 * \file
 * \ingroup MMesh Adaptive
 * \brief   Class defining a newest vertex bisection refinement strategy.
 */

#ifndef DUNE_MMESH_REMESHING_NEWESTVERTEXBISECTION_HH
#define DUNE_MMESH_REMESHING_NEWESTVERTEXBISECTION_HH

#include <memory>
#include <dune/common/exceptions.hh>
#include <dune/grid/common/partitionset.hh>

namespace Dune
{

/*!
 * \ingroup MMesh Adaptive
 * \brief   Class defining a newest vertex bisection refinement strategy.
 */
template<class Grid>
class NewestVertexBisection
{
    static constexpr int dim = Grid::dimension;
    static constexpr int edgeCodim = dim - 1;
    static constexpr int vertexCodim = dim;
    using ctype = typename Grid::ctype;
    using Element = typename Grid::Traits::template Codim<0>::Entity;

public:
    /*!
     * \brief Returns the refinement/coarsening point for each grid cell.
     */

    /*!
     * \brief return refinement point (center of edge opposite to newest vertex)
     *
     * \param element A grid element
     */
    template<class Element>
    static auto refinement (const Element& element)
    {
      assert( dim <= 2 ); // TODO 3d

      if ( dim == 1 )
      {
        const auto& oppositeEdge = element.template subEntity<edgeCodim>( 0 );
        return std::make_pair( oppositeEdge, oppositeEdge.geometry().center() );
      }

      std::size_t maxInsertionLevel = 0;
      for( int i = 0; i < element.subEntities(vertexCodim); ++i )
      {
        const auto& v = element.template subEntity<vertexCodim>(i);
        const std::size_t insertionLevel = v.impl().insertionLevel();
        maxInsertionLevel = std::max( maxInsertionLevel, insertionLevel );
      }

      // if multiple edges are opposite of a vertex with maxInsertionLevel, take the longest one
      std::vector<int> edges;
      for( int i = 0; i < element.subEntities(vertexCodim); ++i )
      {
        const auto& v = element.template subEntity<vertexCodim>(i);
        const std::size_t insertionLevel = v.impl().insertionLevel();
        if (insertionLevel == maxInsertionLevel)
          edges.push_back( dim - i );
      }

      std::sort( edges.begin(), edges.end(),
        [&element]( int i, int j )
        {
          return element.template subEntity<edgeCodim>( i ).geometry().volume()
            > element.template subEntity<edgeCodim>( j ).geometry().volume();
        }
      );

      const auto& oppositeEdge = element.template subEntity<edgeCodim>( edges[0] );
      return std::make_pair( oppositeEdge, oppositeEdge.geometry().center() );
    }

    /*!
     * \brief return coarsening vertex (vertex with highest insertion index)
     *
     * \param element A grid element
     */
    template<class Element>
    static auto coarsening (const Element& element)
    {
      std::size_t maxInsertionLevel = 0;
      int vertexIdx = -1;

      for( int i = 0; i < element.subEntities(vertexCodim); ++i )
      {
        const auto& v = element.template subEntity<vertexCodim>(i);
        if( !v.impl().isInterface() )
        {
          const std::size_t insertionLevel = v.impl().insertionLevel();
          if (insertionLevel >= maxInsertionLevel)
          {
            maxInsertionLevel = insertionLevel;
            vertexIdx = i;
          }
        }
      }

      if ( vertexIdx >= 0 )
        return element.template subEntity<vertexCodim>(vertexIdx);
      else
      {
        DUNE_THROW(GridError, "No vertex could be used for coarsening as they are all part of the interface or boundary.");
        return element.template subEntity<vertexCodim>(0); // dummy!
      }
    }
};

} // end namespace Dumux

#endif
