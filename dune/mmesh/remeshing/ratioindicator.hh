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
 *          By default, we take 2x length of the longest edge contained in the interface as maximal edge length.
 */
template<class Grid>
class RatioIndicator
{
  static constexpr int dim = Grid::dimensionworld;
  using ctype = typename Grid::ctype;
  using GlobalCoordinate = FieldVector<ctype, dim>;
  using Element = typename Grid::Traits::template Codim<0>::Entity;

public:
    /*!
     * \brief Calculates the indicator for each grid cell.
     *
     * \param grid     The grid implementation
     */
    RatioIndicator(const Grid& grid, ctype ratio = 5.0, ctype maxH = -1.0, ctype minH = 1e100)
     : maxRatio_( ratio ), maxH_( maxH )
    {
      if ( maxH < 0.0 )
        for ( const auto& edge : edges( grid.leafGridView() ) )
          if ( grid.isInterface( edge ) )
            maxH_ = std::max( maxH_, 2 * edge.geometry().volume() );

      if ( minH == 1e100 )
        for ( const auto& edge : edges( grid.leafGridView() ) )
          if ( grid.isInterface( edge ) )
            minH_ = std::min( minH_, 0.5 * edge.geometry().volume() );
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
    int operator() (const Element& element) const
    {
      static constexpr int vertexCodim = Element::dimension;
      static constexpr int edgeCodim = Element::dimension - 1;

      bool allInterfaceOrBoundary = true;
      for( std::size_t i = 0; i < element.subEntities(vertexCodim); ++i )
      {
        const auto& v = element.template subEntity<vertexCodim>(i);
        allInterfaceOrBoundary &= ( v.impl().hostEntity()->info().isInterface );
      }

      if (!allInterfaceOrBoundary)
      {
        const auto& geo = element.geometry();

        // ratio criterion for coarsening
        ctype sumE = 0.0;
        for( std::size_t i = 0; i < element.subEntities(1); ++i )
          sumE += element.template subEntity<1>(i).geometry().volume();

        const ctype innerRadius = dim * geo.volume() / sumE;
        const ctype outerRadius = (geo.impl().circumcenter() - geo.corner(0)).two_norm();
        const ctype ratio = outerRadius / (innerRadius * dim);

        if (ratio > maxRatio_)
          return -1;
      }

      // refine by edge length criterion
      for( std::size_t i = 0; i < element.subEntities(edgeCodim); ++i )
      {
        const auto& edge = element.impl().template subEntity<edgeCodim>(i);

        const ctype edgeLength = edge.geometry().volume();
        if ( edgeLength > maxH_ )
          return 1;
      }

      return 0;
    }

    ctype maxRatio () const
    {
      return maxRatio_;
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
    const ctype maxRatio_;
    ctype maxH_ = -1.0;
    ctype minH_ = 1e100;
};

} // end namespace Dumux

#endif
