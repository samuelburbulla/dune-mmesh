// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*!
 * \file
 * \ingroup MMesh Remeshing
 * \brief   Class for computing the distance to the interface.
 */

#ifndef DUNE_MMESH_REMESHING_DISTANCE_HH
#define DUNE_MMESH_REMESHING_DISTANCE_HH

#include <memory>
#include <dune/common/exceptions.hh>
#include <dune/common/timer.hh>
#include <dune/grid/common/partitionset.hh>

namespace Dune
{

/*!
 * \ingroup MMesh Remeshing
 * \brief   Class for computing the distance to the interface.
 */
template<class Grid>
class Distance
{
  static constexpr int dim = Grid::dimensionworld;
  using ctype = typename Grid::ctype;
  using GlobalCoordinate = FieldVector<ctype, dim>;
  using Vertex = typename Grid::Vertex;
  using Element = typename Grid::template Codim<0>::Entity;

public:
    /*!
     * \brief Computes the distance to the interface for each vertex
     *
     */
    Distance(const Grid& grid)
     : grid_( grid ), indexSet_( grid_.leafIndexSet() )
    {
      update();
    }

    //! Update the distances of all vertices
    void update()
    {
      // Timer timer;
      distances_.resize( indexSet_.size(dim) );
      std::fill( distances_.begin(), distances_.end(), 1e100 );

      for ( const auto& vertex : vertices( grid_.leafGridView() ) )
      {
        const auto& idx = indexSet_.index( vertex );

        if ( vertex.impl().isInterface() )
          distances_[ idx ] = 0.0;

        for ( const auto& ivertex : vertices( grid_.interfaceGrid().leafGridView() ) )
        {
          ctype dist = (ivertex.geometry().center() - vertex.geometry().center()).two_norm();
          distances_[ idx ] = std::min( distances_[ idx ], dist );
        }
      }
      // std::cout << "Distance::update() took " << timer.stop() << "s." << std::endl;
    };

    /*!
     * \brief function call operator to return distance of vertex
     *
     * \param vertex    A grid vertex
     */
    ctype operator() (const Vertex& vertex) const
    {
      return distances_[ indexSet_.index( vertex ) ];
    }

    /*!
     * \brief function call operator to return distance of element (average of the vertex distances)
     *
     * \param element    A grid element
     */
    ctype operator() (const Element& element) const
    {
      ctype dist = 0.0;
      for ( std::size_t i = 0; i < dim+1; ++i )
      {
        const auto& vertex = element.template subEntity<dim>(i);
        dist += distances_[ indexSet_.index( vertex ) ];
      }
      dist /= dim+1;
      return dist;
    }

    /*!
     * \brief function call operator to return distance
     *
     * \param index    Index of a vertex
     */
    template< class Index >
    ctype operator[] (const Index& index) const
    {
      return distances_[ index ];
    }

    //! return size of distances vector
    std::size_t size() const
    {
      return distances_.size();
    }

  private:
    std::vector<ctype> distances_;
    const Grid& grid_;
    const typename Grid::LeafIndexSet& indexSet_;
};

} // end namespace Dumux

#endif
