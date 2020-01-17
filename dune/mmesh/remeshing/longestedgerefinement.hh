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
     * \brief return coarsening vertex (vertex of shortest edge with highest insertion index)
     *
     * \param element A grid element
     */
    template<class Element>
    static auto coarsening (const Element& element)
    {
      for( std::size_t i = 0; i < element.subEntities(vertexCodim); ++i )
      {
        const auto& v = element.template subEntity<vertexCodim>(i);
        if( !v.impl().isInterface() )
          return v;
      }

      DUNE_THROW(GridError, "No vertex could be used for coarsening as they are all part of the interface or boundary.");
      return element.template subEntity<vertexCodim>(0); // dummy
    }
};

} // end namespace Dumux

#endif
