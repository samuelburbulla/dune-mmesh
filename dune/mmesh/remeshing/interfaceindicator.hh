// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*!
 * \file
 * \ingroup MMesh Adaptive
 * \brief   Class defining an indicator for grid remeshing regarding the cell volume.
 */

#ifndef DUNE_MMESH_REMESHING_INTERFACEINDICATOR_HH
#define DUNE_MMESH_REMESHING_INTERFACEINDICATOR_HH

#include <memory>
#include <dune/common/exceptions.hh>
#include <dune/grid/common/partitionset.hh>

namespace Dune
{

/*!
 * \ingroup MMesh Adaptive
 * \brief   Class defining an indicator for grid adaptation regarding the interface.
 */
template<class Grid>
class InterfaceIndicator
{
  static constexpr int dim = Grid::dimensionworld;
  using ctype = typename Grid::ctype;
  using GlobalCoordinate = FieldVector<ctype, dim>;
  using Element = typename Grid::Traits::template Codim<0>::Entity;

public:
    /*!
     * \brief Calculates the indicator for each grid cell.
     *
     * \param minH     The minimal allowed interface length
     * \param maxH     The maximal allowed interface length
     */
    InterfaceIndicator(const Grid& grid, ctype minH, ctype maxH)
     : grid_(grid),
       minH_(minH),
       maxH_(maxH)
    {}

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
      static constexpr int edgecd = Element::dimension - 1;

      const auto& vol = element.geometry().volume();

      for( int i = 0; i < element.subEntities(edgecd); ++i )
      {
        const auto& edge = element.template subEntity<edgecd>(i);

        if( grid_.isInterface( edge ) )
        {
          ctype edgelength = edge.geometry().volume();

          if ( edgelength > maxH_ )
            return 1;

          if ( edgelength < minH_ )
            return -1;

          const ctype cellheight = vol / edgelength;
          if( cellheight < 0.5 * edgelength )
            return -1;
        }
      }

      return 0;
    }

    ctype minH() const
    {
      return minH_;
    }

    const Grid& grid_;
    const ctype minH_, maxH_;
};

} // end namespace Dumux

#endif
