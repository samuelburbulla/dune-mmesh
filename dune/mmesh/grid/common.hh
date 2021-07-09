// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_GRID_COMMON_HH
#define DUNE_MMESH_GRID_COMMON_HH

/** \file
 * \brief Some common helper methods
 */

namespace Dune
{

  //! Hash a UInt vector
  struct HashUIntVector {
    std::size_t operator() (const std::vector<std::size_t>& a) const
    {
      std::size_t hash = std::hash<std::size_t>{}(a[0]);
      for ( std::size_t i = 1; i < a.size(); ++i )
        hash = hash ^ (std::hash<std::size_t>{}(a[i]) << i);
      return hash;
    }
  };

  //! Hash a UInt array
  struct HashUIntArray {
    template< std::size_t dim >
    std::size_t operator() (const std::array<std::size_t, dim>& a) const
    {
      std::size_t hash = std::hash<std::size_t>{}(a[0]);
      for ( std::size_t i = 1; i < a.size(); ++i )
        hash = hash ^ (std::hash<std::size_t>{}(a[i]) << i);
      return hash;
    }
  };

  namespace MMeshImpl
  {

    template < int dim >
    class sort_indices
    {
      public:
        sort_indices(std::array< std::size_t, dim > ids)
         : ids_(ids) {}

        bool operator()(std::size_t i, std::size_t j) const
        {
          return ids_[i] < ids_[j];
        }

        private:
          std::array< std::size_t, dim > ids_;
    };

    //! return reference element
    template< int dim, typename ctype = double >
    static inline auto& ref()
    {
      return ReferenceElements<ctype, dim>::simplex();
    }

    //! Return list of indices sorted by id
    template < typename HostEntity, int dim >
    static inline auto computeCGALIndices( const HostEntity& hostEntity )
    {
      std::array< std::size_t, dim+1 > ids, indices;
      for ( std::size_t i = 0; i < dim+1; ++i )
      {
        ids[i] = hostEntity->vertex( i )->info().id;
        indices[i] = i;
      }
      std::sort(indices.begin(), indices.end(), sort_indices<dim+1>(ids));

      return indices;
    }

    // for a given dune facet index compute corresponding CGAL .second value
    template< std::size_t dim >
    static inline std::size_t duneFacetToCgalSecond ( const std::size_t duneFacet, const std::array< std::size_t, dim+1 >& cgalIndex )
    {
      std::size_t sum = 0;
      for ( std::size_t k = 0; k < dim; ++k )
        sum += cgalIndex[ ref<dim>().subEntity(duneFacet, 1, k, dim) ];

      static constexpr int max = (dim == 2) ? 3 : 6;
      return max - sum;
    }

    // for a given CGAL .second value compute corresponding dune facet index
    template< std::size_t dim, typename HostFacet >
    static inline std::size_t cgalFacetToDuneFacet ( const HostFacet& facet )
    {
      const auto& i = facet.second;
      const auto& cgalIndex = facet.first->info().cgalIndex;

      // invert cgalIndex
      auto duneIndex = cgalIndex;
      for ( std::size_t k = 0; k < dim+1; ++k )
        duneIndex[ cgalIndex[k] ] = k;

      std::size_t sum = 0;
      for ( std::size_t k = 0; k < dim; ++k )
        sum += duneIndex[ (i+k+1)%(dim+1) ];

      static const int thr = (dim == 2) ? 1 : 3;
      return sum - thr;
    }

    template< std::size_t dim, typename HostEdge >
    static inline std::size_t cgalEdgeToDuneEdge( const HostEdge& cgalEdge )
    {
      const auto& c = cgalEdge.first;
      const auto& i = cgalEdge.second;
      const auto& j = cgalEdge.third;

      const auto& cgalIndex = c->info().cgalIndex;

      // invert cgalIndex
      auto duneIndex = cgalIndex;
      for ( std::size_t k = 0; k < dim+1; ++k )
        duneIndex[ cgalIndex[k] ] = k;

      auto i0 = duneIndex[ i ];
      auto j0 = duneIndex[ j ];

      if ( i0 > j0 )
        std::swap(i0, j0);

      if (j0 != 3)
        return i0+j0-1;
      else
        return 3+i0;
    }

  }

}

#endif
