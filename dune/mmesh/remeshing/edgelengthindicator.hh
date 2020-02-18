// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*!
 * \file
 * \ingroup MMesh Adaptive
 * \brief   Class defining an indicator for grid remeshing regarding a maximal and minimal edge length.
 */

#ifndef DUNE_MMESH_REMESHING_EDGELENGTHINDICATOR_HH
#define DUNE_MMESH_REMESHING_EDGELENGTHINDICATOR_HH

#include <memory>
#include <dune/common/exceptions.hh>
#include <dune/grid/common/partitionset.hh>

namespace Dune
{

/*!
 * \ingroup MMesh Adaptive

 * \brief   Class defining an indicator for grid adaptation regarding a maximal/minimal edge.
 */
template<class Grid>
class EdgeLengthIndicator
{
  static constexpr int dim = Grid::dimensionworld;
  static constexpr int vertexCodim = dim;
  static constexpr int edgeCodim = dim - 1;
  using ctype = typename Grid::ctype;
  using Element = typename Grid::Traits::template Codim<0>::Entity;

public:
    /*!
     * \brief Calculates the indicator for each grid cell.
     *
     * \param minH     The minimal allowed edge length
     * \param maxH     The maximal allowed edge length
     */
    EdgeLengthIndicator(const Grid& grid, ctype minH = 1e100, ctype maxH = -1e100)
     : minH_(minH), maxH_(maxH)
    {
      if ( maxH == -1e100 || minH == 1e100 )
        for ( const auto& edge : edges( grid.leafGridView() ) )
          if ( grid.isInterface( edge ) )
          {
            if ( maxH == -1e100 )
              maxH_ = std::max( maxH_, 2.1 * edge.geometry().volume() );
            if ( minH == 1e100 )
              minH_ = std::min( minH_, 0.5 * edge.geometry().volume() );
          }
    }

    /*!
     * \brief function call operator to return mark
     *
     * \return  1 if an element should be refined
     *         -1 if an element should be coarsened
     *          0 otherwise
     *
     * \param element A grid element
     */
    int operator() (const Element& element) const
    {
      int refine = 0;
      int coarse = 0;

      for( int i = 0; i < element.subEntities(edgeCodim); ++i )
      {
        const auto& edge = element.template subEntity<edgeCodim>(i);

        // skip edges where all vertices are interface
        int count = 0;
        for( std::size_t j = 0; j < edge.subEntities(vertexCodim); ++j )
          if( edge.impl().template subEntity<vertexCodim>(j).impl().isInterface() )
            count++;
        if ( count == edge.subEntities(vertexCodim) )
          continue;

        const ctype edgeLength = edge.geometry().volume();

        if (edgeLength < minH_)
          coarse++;

        if (edgeLength > maxH_)
          refine++;
      }

      if (refine == 0 && coarse == 0)
        return 0;

      // first coarsen
      if (coarse > 0)
        return -1;

      // second refine
      if (refine > 0)
        return 1;

      return 0;
    }

    ctype maxH () const
    {
      return maxH_;
    }

    ctype minH () const
    {
      return minH_;
    }

  private:
    ctype minH_;
    ctype maxH_;
};

} // end namespace Dumux

#endif
