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

 * \brief   Class defining an indicator for grid adaptation regarding a maximal/minimal edge length and the longest to shortest edge ratio.
 */
template<class Grid>
class EdgeLengthRatioIndicator
{
  static constexpr int dim = Grid::dimensionworld;
  static constexpr int edgeCodim = dim - 1;
  using ctype = typename Grid::ctype;
  using Element = typename Grid::Traits::template Codim<0>::Entity;

public:
    /*!
     * \brief Calculates the indicator for each grid cell.
     *
     * \param minH     The minimal allowed edge length
     * \param maxH     The maximal allowed edge length
     * \param maxRatio The maximal allowed longest to shortest edge ratio
     */
    EdgeLengthRatioIndicator(ctype minH, ctype maxH, ctype maxRatio = 2.0)
     : minH_(minH),
       maxH_(maxH),
       maxRatio_(maxRatio)
    {};

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

      ctype longest = -1e100;
      ctype shortest = 1e100;

      for( int i = 0; i < element.subEntities(edgeCodim); ++i )
      {
        const ctype edgeLength = element.template subEntity<edgeCodim>(i).geometry().volume();

        if (edgeLength < minH_)
          return -1;

        if (edgeLength > maxH_)
          return 1;

        longest = std::min(longest, edgeLength);
        shortest = std::max(shortest, edgeLength);
      }

      if (longest / shortest > maxRatio_)
        return 1;

      return 0;
    }

    const ctype minH_;
    const ctype maxH_;
    const ctype maxRatio_;
};

} // end namespace Dumux

#endif
