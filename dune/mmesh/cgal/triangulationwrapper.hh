// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_CGAL_TRIANGULATIONWRAPPER_HH
#define DUNE_MMESH_CGAL_TRIANGULATIONWRAPPER_HH

/** \file
 * \brief A CGAL triangulation wrapper class
 */
#include <dune/grid/common/exceptions.hh>
#include <dune/mmesh/grid/pointfieldvector.hh>

namespace Dune
{

  template< class Coordinate >
  class Plane
  {
  public:
    Plane( const Coordinate& a, const Coordinate& b )
    {
      n_ = a;
      n_ -= b;
      std::swap( n_[0], n_[1] );
      n_[0] *= -1.0;
      n_ /= n_.two_norm();

      p_ = a;
    }

    template< class Vertex >
    Plane( const Vertex& a, const Vertex& b )
     : Plane( makeFieldVector( a->point() ), makeFieldVector( b->point() ) )
    {}

    auto signedDistance( const Coordinate& x ) const
    {
      return n_ * (x - p_);
    }

    template< class Vertex >
    auto signedDistance( const Vertex& vh ) const
    {
      return signedDistance( makeFieldVector( vh->point() ) );
    }

  private:
    Coordinate p_, n_;
  };

  //! TriangulationWrapper<2>
  template<>
  class TriangulationWrapper<2> :
    public MMeshDefaults::Triangulation<2>::type
  {
      using ThisType = TriangulationWrapper<2>;
      using BaseType = typename MMeshDefaults::Triangulation<2>::type;

      using CellHandle = typename BaseType::Face_handle;
      using FacetHandle = std::pair<typename BaseType::Face_handle, int>;
      using VertexHandle = typename BaseType::Vertex_handle;

      using Vertex_pair = std::pair< VertexHandle, VertexHandle >;
      using Vertex_pair_Facet_map = std::map< Vertex_pair, FacetHandle >;
      using Vertex_handle_unique_hash_map = std::map< VertexHandle, VertexHandle >;

  public:
    template< class VertexHandle, class Elements >
    void removeAndGiveNewElements (const VertexHandle& vh, Elements& elements)
    {
      std::list<FacetHandle> hole;
      this->make_hole(vh, hole);
      this->fill_hole_delaunay(hole, std::back_inserter(elements));
      this->delete_vertex(vh);
    }

    template< class VertexHandle, class Elements >
    void remeshHoleConstrained(const VertexHandle& vh, std::list<FacetHandle>& hole, Elements& elements, const std::vector<VertexHandle> constraint)
    {
      assert( constraint.size() == 2 );

      std::list<FacetHandle> hole1;
      std::list<FacetHandle> hole2;

      Plane< FieldVector<double, 2> > plane ( constraint[0], constraint[1] );

      VertexHandle right;
      VertexHandle left;

      auto helperface = create_face();
      FacetHandle helperfacet (helperface, 0);

      // collect facets corresponding to their side of the constraint
      for( const FacetHandle& facet : hole )
      {
        const auto& vh = facet.first->vertex( (facet.second+1)%3 );

        if ( vh == constraint[0] || vh == constraint[1] )
        {
          const auto& vh2 = facet.first->vertex( (facet.second+2)%3 );
          if ( plane.signedDistance( vh2 ) > 0 )
          {
            hole1.push_back( facet );
            left = vh;
            hole1.push_back( helperfacet );
          }
          else
          {
            hole2.push_back( facet );
            right = vh;
            hole2.push_back( helperfacet );
          }
        }
        else
        {
          if ( plane.signedDistance( vh ) > 0 )
            hole1.push_back( facet );
          else
            hole2.push_back( facet );
        }
      }

      // fill hole1
      helperface->set_vertices( VertexHandle(), right, left );
      this->fill_hole(vh, hole1, std::back_inserter(elements));

      // fill hole2
      helperface->set_vertices( VertexHandle(), left, right );
      this->fill_hole(vh, hole2, std::back_inserter(elements));

      // glue the facets at helperface together
      Face_handle f1; int i1;
      for ( auto f : elements )
        for ( int i = 0; i < 3; ++i )
          if (f->neighbor(i) == helperface)
          {
            if (f1 == Face_handle())
            {
              f1 = f;
              i1 = i;
            }
            else
            {
              f->set_neighbor(i, f1);
              f1->set_neighbor(i1, f);

              delete_vertex(vh);
              delete_face(helperface);
              return;
            }
          }

      DUNE_THROW( GridError, "Could not remesh hole with interface." );
    }

  private:
    Vertex_pair make_vertex_pair ( FacetHandle f ) const
    {
      Vertex_pair vp;
      vp.first = f.first->vertex( (f.second+1)%3 );
      vp.second = f.first->vertex( (f.second+2)%3 );

      if ( vp.first > vp.second )
        std::swap( vp.first, vp.second );

      return vp;
    }
  };


  //! TriangulationWrapper<3>
  template<>
  class TriangulationWrapper<3> :
    public MMeshDefaults::Triangulation<3>::type
  {
    using ThisType = TriangulationWrapper<3>;
    using BaseType = typename MMeshDefaults::Triangulation<3>::type;

  public:
    template< class VertexHandle, class Elements >
    void removeAndGiveNewElements (const VertexHandle& v, Elements fit)
    {
      DUNE_THROW( GridError, "Removal is not supported in 3D!" );
    }
  }; // end of TriangulationWrapper<3>


  //! DelaunayTriangulationWrapper
  template<int dim>
  class DelaunayTriangulationWrapper :
    public MMeshDefaults::Delaunay<dim>::type
  {
      using ThisType = DelaunayTriangulationWrapper<dim>;
      using BaseType = typename MMeshDefaults::Delaunay<dim>::type;
  };

} // namespace Dune

#endif
