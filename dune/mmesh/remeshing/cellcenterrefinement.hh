// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*!
 * \file
 * \ingroup MMesh Adaptive
 * \brief   Class defining a longest edge refinement strategy.
 */

#ifndef DUNE_MMESH_REMESHING_CELLCENTERREFINEMENT_HH
#define DUNE_MMESH_REMESHING_CELLCENTERREFINEMENT_HH

#include <memory>
#include <dune/common/exceptions.hh>
#include <dune/grid/common/partitionset.hh>

namespace Dune
{

/*!
 * \ingroup MMesh Adaptive
 * \brief   Class defining a cell center refinement strategy.
 */
template<class Grid>
class CellCenterRefinement
{
    using Element = typename Grid::Traits::template Codim<0>::Entity;
    using Vertex = typename Grid::Vertex;

public:
    /*!
     * \brief Returns the refinement/coarsening point for each grid cell.
     */

    /*!
     * \brief return refinement point (center of cell)
     *
     * \param element A grid element
     */
    template<class Element>
    static auto refinement (const Element& element)
    {
      return std::make_pair( element, element.geometry().center() );
    }

    /*!
     * \brief return coarsening vertex (vertex of shortest edge with highest insertion index)
     *
     * \param element A grid element
     */
    template<class Element>
    static Vertex coarsening (const Element& element)
    {
      DUNE_THROW(NotImplemented, "Coarsening of interface!");
    }
};

} // end namespace Dumux

#endif
