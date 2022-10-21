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
#include <dune/mmesh/remeshing/distance.hh>

namespace Dune
{

/*!
 * \ingroup MMesh Adaptive
 * \brief   Class defining an indicator for grid remeshing regarding the edge length ratio.
 *          By default, we take 2x length of the longest edge contained in the interface as maximal edge length and 0.5x length of the shortest edge as minimal edge length.
 */
template<class Grid>
class RatioIndicator
{
  static constexpr int dim = Grid::dimensionworld;
  using ctype = typename Grid::ctype;
  using GlobalCoordinate = FieldVector<ctype, dim>;

public:
  using DistanceType = Distance<Grid>;

    /*!
     * \brief Calculates the indicator for each grid cell.
     *
     * \param h                The objective edge length (aims at edge length in [h/4, 2*h]).
     * \param distProportion   Cells with distance to interface of value greater than distProportion * max(dist) are refined to ...
     * \param factor           ... edge length in [factor * minH, factor * maxH].
     */
    RatioIndicator( ctype h = 0.0, ctype distProportion = 1.0, ctype factor = 1.0 )
     : edgeRatio_( 4. ),                  // the edge ratio. Decreases minH!
       K_( 2. ),                          // the factor for the maximal edge length. Defaults to 2, because we use edge bisection.
       maxH_( K_ * h ),
       k_( K_ / ( 2. * edgeRatio_ ) ),    // ensure that a triangle with two edges longer than maxH splits up to a triangle where all edges are longer than minH
       minH_( k_ * h ),
       radiusRatio_( 30. ),               // additional radius ratio to avoid very ugly cells
       distProportion_( distProportion ), // cells with distance to interface of value greater than distProportion_ * max(dist) are refined to ...
       factor_( factor )                  // ... edge length in [factor * minH_, factor * maxH_]
    {}

    //! Calculates minH_ and maxH_ for the current interface edge length and sets factor_ to maxh / minh.
    void init(const Grid& grid)
    {
      maxH_ = 0.0;
      minH_ = 1e100;

      for ( const auto& edge : edges( grid.interfaceGrid().leafGridView() ) )
      {
        const ctype h = edge.geometry().volume();
        maxH_ = std::max( maxH_, h );
        minH_ = std::min( minH_, h );
      }

      maxH_ = K_ * maxH_;
      minH_ = k_ * minH_;

      ctype maxh = 0.0;
      ctype minh = 1e100;
      for ( const auto& edge : edges( grid.leafGridView() ) )
      {
        const ctype h = edge.geometry().volume();
        maxh = std::max( maxh, h );
        minh = std::min( minh, h );
      }

      // fallback if there is no interface
      if ( grid.interfaceGrid().size(1) == 0 )
      {
        maxH_ = K_ * maxh;
        minH_ = k_ * minh;
      }

      factor_ = maxh / minh;
      distance_ = DistanceType(grid);
    };

    //! Update the distances of all vertices
    void update()
    {
      distance_.update();
      maxDist_ = distProportion_ * distance_.maximum();
    }

    /*!
     * \brief Function call operator to return mark
     *
     * \return  1 if an element should be refined,
     *         -1 if an element should be coarsened,
     *          0 otherwise.
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
      ctype dist = std::min( maxDist_, distance_( element ) );
      const ctype l = dist / maxDist_;
      const ctype minH = (1. - l) * minH_ + l * factor_ * minH_;
      const ctype maxH = (1. - l) * maxH_ + l * factor_ * maxH_;

      for( std::size_t i = 0; i < element.subEntities(edgeCodim); ++i )
      {
        const auto& edge = element.template subEntity<edgeCodim>(i);
        const ctype len = edge.geometry().volume();

        if (len < minH)
          coarse++;

        if (len > maxH)
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

    //! Returns maxH
    ctype maxH () const
    {
      return maxH_;
    }

    //! Returns reference to maxH
    ctype& maxH ()
    {
      return maxH_;
    }

    //! Returns minH
    ctype minH () const
    {
      return minH_;
    }

    //! Returns reference to minH
    ctype& minH ()
    {
      return minH_;
    }

    ctype radiusRatio () const
    {
      return radiusRatio_;
    }

    //! Returns reference to distProportion
    ctype& distProportion ()
    {
      return distProportion_;
    }

    //! Returns reference to factor
    ctype& factor ()
    {
      return factor_;
    }

    //! Returns distance object
    const DistanceType& distance () const
    {
      if (!distance_.initialized())
        distance_.update();
      return distance_;
    }

  private:
    const ctype edgeRatio_;
    const ctype K_;
    ctype maxH_;
    const ctype k_;
    ctype minH_;
    const ctype radiusRatio_;
    ctype distProportion_;
    ctype factor_;
    ctype maxDist_;
    mutable DistanceType distance_;
};

} // end namespace Dune

#endif
