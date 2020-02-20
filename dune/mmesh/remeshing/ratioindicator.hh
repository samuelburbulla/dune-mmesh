// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*!
 * \file
 * \ingroup MMesh Adaptive
 * \brief   Class defining an indicator for grid remeshing regarding the ratio of outer to inner radius.
 */

#ifndef DUNE_MMESH_REMESHING_RATIOINDICATOR_HH
#define DUNE_MMESH_REMESHING_RATIOINDICATOR_HH

#include <memory>
#include <dune/common/exceptions.hh>
#include <dune/grid/common/partitionset.hh>

namespace Dune
{

/*!
 * \ingroup MMesh Adaptive
 * \brief   Class defining an indicator for grid remeshing regarding the ratio of outer to inner sphere radius (normalized by factor 1/dimension).
 *          By default, we take 2x length of the longest edge contained in the interface as maximal edge length and 0.5x length of the shortest edge as minimal edge length.
 */
template<class Grid>
class RatioIndicator
{
  static constexpr int dim = Grid::dimensionworld;
  using ctype = typename Grid::ctype;
  using GlobalCoordinate = FieldVector<ctype, dim>;

public:
    /*!
     * \brief Calculates the indicator for each grid cell.
     *
     * \param h        The objective edge length (aims at edge length in [h/4, 2*h])
     */
    RatioIndicator( ctype h = 0.0 )
     : edgeRatio_( 4. ),                  // the edge ratio. Decreases minH!
       K_( 2. ),                          // the factor for the maximal edge length. Defaults to 2, because we use edge bisection.
       maxH_( K_ * h ),
       k_( K_ / ( 2. * edgeRatio_ ) ),    // ensure that a triangle with two edges longer than maxH splits up to a triangle where all edges are longer than minH
       minH_( k_ * h ),
       radiusRatio_( 30. )               // additional radius ratio to avoid very ugly cells
    {}

    /*!
     * \brief Calculates minH and maxH for the current interface edge length.
     *
     * \param grid     The grid implementation
     */
    void init(const Grid& grid)
    {
      for ( const auto& edge : edges( grid.interfaceGrid().leafGridView() ) )
      {
        const ctype h = edge.geometry().volume();
        maxH_ = std::max( maxH_, K_ * h );
        minH_ = std::min( minH_, k_ * h );
      }
    };

    /*!
     * \brief function call operator to return mark
     *
     * \return  1 if an element should be refined
     *         -1 if an element should be coarsened
     *          0 otherwise
     *
     * \param element A grid element
     */
    template< class Element >
    int operator() (const Element& element) const
    {
      static constexpr int edgeCodim = Element::dimension - 1;

      int refine = 0;
      int coarse = 0;

      ctype minE =  1e100;
      ctype maxE = -1e100;

      ctype sumE = 0.0;

      // edge length
      for( int i = 0; i < element.subEntities(edgeCodim); ++i )
      {
        const auto& edge = element.template subEntity<edgeCodim>(i);
        const ctype len = edge.geometry().volume();

        if (len < minH_)
          coarse++;

        if (len > maxH_)
          refine++;

        minE = std::min( minE, len );
        maxE = std::max( maxE, len );

        sumE += len;
      }

      // edge ratio criterion
      const ctype edgeRatio = maxE / minE;
      if (edgeRatio > edgeRatio_)
        coarse++;

      // radius ratio criterion
      const auto& geo = element.geometry();
      const ctype innerRadius = dim * geo.volume() / sumE;
      const ctype outerRadius = (geo.impl().circumcenter() - geo.corner(0)).two_norm();
      const ctype radiusRatio = outerRadius / (innerRadius * dim);

      if (radiusRatio > radiusRatio_)
        coarse++;

      // priority on coarse
      if (coarse > 0)
        if ( dim != 3 ) // disable coarsening in 3d until it is implemented
          return -1;

      // then refine
      if (refine > 0)
        return 1;

      // nothing to do
      return 0;
    }

    ctype edgeRatio () const
    {
      return edgeRatio_;
    }

    ctype maxH () const
    {
      return maxH_;
    }

    ctype minH () const
    {
      return minH_;
    }

    ctype radiusRatio () const
    {
      return radiusRatio_;
    }

  private:
    const ctype edgeRatio_;
    const ctype K_;
    ctype maxH_;
    const ctype k_;
    ctype minH_;
    const ctype radiusRatio_;
};

} // end namespace Dumux

#endif
