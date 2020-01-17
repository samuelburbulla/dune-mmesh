// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*!
 * \file
 * \ingroup MMesh Adaptive
 * \brief   Class defining an indicator for grid remeshing regarding the cell volume.
 */

#ifndef DUNE_MMESH_REMESHING_VOLUMEINDICATOR_HH
#define DUNE_MMESH_REMESHING_VOLUMEINDICATOR_HH

#include <memory>
#include <dune/common/exceptions.hh>
#include <dune/grid/common/partitionset.hh>

namespace Dune
{

/*!
 * \ingroup MMesh Adaptive
 * \brief   Class defining an indicator for grid adaptation regarding the cell volume.
 */
template<class Grid>
class VolumeIndicator
{
  static constexpr int dim = Grid::dimensionworld;
  using ctype = typename Grid::ctype;
  using GlobalCoordinate = FieldVector<ctype, dim>;
  using Element = typename Grid::Traits::template Codim<0>::Entity;

public:
    /*!
     * \brief Calculates the indicator for each grid cell.
     *
     * \param minVol     The minimal allowed volume
     * \param maxH       The maximal allowed edge length
     * \param maxRation  The maximal allowed ratio between circumradius and inradius
     */
    VolumeIndicator(ctype minVol,
                    const std::function<ctype(GlobalCoordinate)>& maxH,
                    ctype maxRatio = 100.)
     : minVol_(minVol),
       maxH_(maxH),
       maxRatio_(maxRatio)
    {
      // assure that maxRatio is greater than minimal possible ratio
      assert(maxRatio > dim);
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

      bool allInterface = true;
      for( int i = 0; i < element.subEntities(vertexCodim); ++i )
        allInterface &= element.template subEntity<vertexCodim>(i).impl().hostEntity()->info().isInterface;

      if (!allInterface)
      {
        const auto& geo = element.geometry();
        // coarsen
        const ctype volume = geo.volume();
        if (volume < minVol_)
          return -1;

        // ratio criterion
        ctype sumE = 0.0;
        for( int i = 0; i < element.subEntities(1); ++i )
          sumE += element.template subEntity<1>(i).geometry().volume();

        const ctype innerRadius = dim * volume / sumE;
        const ctype outerRadius = (geo.impl().circumcenter() - geo.corner(0)).two_norm();
        const ctype ratio = outerRadius / innerRadius;

        if (ratio > maxRatio_)
          return -1;
      }

      // refine
      for( int i = 0; i < element.subEntities(edgeCodim); ++i )
      {
        const auto& geo = element.impl().template subEntity<edgeCodim>(i).geometry();
        const ctype edgeLength = geo.volume();
        if ( edgeLength > maxH_(geo.center()) )
          return 1;
      }

      return 0;
    }

    const ctype minVol_;
    const std::function<ctype(GlobalCoordinate)> maxH_;
    const ctype maxRatio_;
};

} // end namespace Dumux

#endif
