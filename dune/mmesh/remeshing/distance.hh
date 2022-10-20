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
  using ThisType = Distance<Grid>;
  static constexpr int dim = Grid::dimensionworld;
  using ctype = typename Grid::ctype;
  using GlobalCoordinate = FieldVector<ctype, dim>;
  using Vertex = typename Grid::Vertex;
  using Element = typename Grid::template Codim<0>::Entity;
  using Facet = typename Grid::template Codim<1>::Entity;
  using InterfaceElement = typename Grid::InterfaceGrid::template Codim<0>::Entity;

public:
    //! Default constructor
    Distance() {}

    //! Constructor with grid reference
    Distance(const Grid& grid)
     : grid_( &grid ), initialized_(false)
    {}

    //! Update the distances of all vertices
    void update()
    {
      // Resize distance_ and set to high default value
      distances_.resize( indexSet().size(dim) );
      std::fill(distances_.begin(), distances_.end(), 1e100);

      // Set all interface vertices to zero and initialize queue

      for ( const InterfaceElement& ielement : elements( grid_->interfaceGrid().leafGridView(), Partitions::interior ) )
      {
        // Convert to bulk facet
        const Facet facet = grid_->entity( ielement.impl().hostEntity() );

        // Compute vertex distances to this facet
        handleFacet(facet);
      }

      initialized_ = true;
    };

    //! Return if distance has been initialized
    bool initialized() const
    {
      return initialized_;
    }

    /*!
     * \brief function call operator to return distance of vertex
     *
     * \param vertex    A grid vertex
     */
    ctype operator() (const Vertex& vertex) const
    {
      assert(initialized_);
      assert( indexSet().index( vertex ) < size() );
      return distances_[ indexSet().index( vertex ) ];
    }

    /*!
     * \brief Set distance of vertex
     *
     * \param vertex    A grid vertex
     * \param value     Distance value
     */
    void set(const Vertex& vertex, ctype value)
    {
      assert( indexSet().index( vertex ) < size() );
      distances_[ indexSet().index( vertex ) ] = value;
    }

    /*!
     * \brief function call operator to return distance of element (average of the vertex distances)
     *
     * \param element    A grid element
     */
    ctype operator() (const Element& element) const
    {
      assert(initialized_);
      ctype dist = 0.0;
      for ( std::size_t i = 0; i < dim+1; ++i )
      {
        const auto& vertex = element.template subEntity<dim>(i);
        dist += distances_[ indexSet().index( vertex ) ];
      }
      dist /= dim+1;
      return dist;
    }

    //! Interface element
    ctype operator() (const InterfaceElement& element) const
    {
      return 0.0;
    }

    /*!
     * \brief function call operator to return distance
     *
     * \param index    Index of a vertex
     */
    template< class Index >
    ctype operator[] (const Index& index) const
    {
      assert(initialized_);
      return distances_[ index ];
    }

    //! return maximum distance
    ctype maximum() const
    {
      assert(initialized_);
      double maximum = 0.0;
      for ( const auto& d : distances_ )
        maximum = std::max( d, maximum );
      return maximum;
    }

    //! return size of distances vector
    std::size_t size() const
    {
      assert(initialized_);
      return distances_.size();
    }

  private:
    //! Handle facet: Compute all vertex distances
    void handleFacet(const Facet& facet)
    {
      for (const auto& v : vertices(grid_->leafGridView(), Partitions::interior))
      {
        ctype dist = computeDistance(v, facet);
        ctype &currentDist = distances_[ indexSet().index( v ) ];
        if (dist < currentDist)
          currentDist = dist;
      }
    }

    //! Compute vertex distance value
    ctype computeDistance(const Vertex& v, const Facet &facet) const
    {
      GlobalCoordinate vp = v.geometry().center();
      const auto& geo = facet.geometry();

      if constexpr (dim == 2)
      {
        for (std::size_t i = 0; i < dim; ++i)
        {
          const GlobalCoordinate& vi = geo.corner(i);
          const GlobalCoordinate& vj = geo.corner(1-i);
          ctype sgn = (vi - vj) * (vp - vi);
          if (sgn >= 0)
            return (vp - vi).two_norm();
        }

        Dune::Plane<GlobalCoordinate> plane (geo.corner(0), geo.corner(1));
        return std::abs(plane.signedDistance(vp));
      }
      else // dim == 3
      {
        // We only use an estimate for 3d.
        ctype dist = 1e100;
        for (std::size_t i = 0; i < dim; ++i)
        {
          const GlobalCoordinate& vi = geo.corner(i);
          dist = std::min(dist, (vp - vi).two_norm2());

          const GlobalCoordinate& vj = geo.corner((i+1)%dim);
          dist = std::min(dist, (vp - 0.5*(vi+vj)).two_norm2());
        }

        const GlobalCoordinate& vc = geo.center();
        dist = std::min(dist, (vp - vc).two_norm2());

        return std::sqrt(dist);
      }

      DUNE_THROW(InvalidStateException, "Compute distance failed");
      return -1.;
    }

    const typename Grid::LeafIndexSet& indexSet() const
    {
      return grid_->leafIndexSet();
    }

    std::vector<ctype> distances_;
    const Grid* grid_;
    bool initialized_;
};

} // end namespace Dune

#endif
