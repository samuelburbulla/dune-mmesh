// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_CGAL_TRIANGULATIONWRAPPER_HH
#define DUNE_MMESH_CGAL_TRIANGULATIONWRAPPER_HH

/** \file
 * \brief A CGAL triangulation wrapper class
 */
#include <dune/mmesh/mmesh.hh>

// CGAL includes
#include <CGAL/Triangulation_2.h>
#include <CGAL/Triangulation_3.h>

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

    auto signedDistance( const Coordinate& x ) const
    {
      return n_ * (x - p_);
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
    void removeAndGiveNewElements (const VertexHandle& vh, const Elements& elements)
    {
      std::list<FacetHandle> hole;
      this->make_hole(vh, hole);
      this->fill_hole_delaunay(hole, elements);
      this->delete_vertex(vh);
    }

    template< class VertexHandle, class Elements >
    void remeshHoleConstrained(const VertexHandle& vh, std::list<FacetHandle>& hole, Elements elements, const std::vector<VertexHandle> constraint)
    {
      assert( constraint.size() == 2 );

      Plane< FieldVector<double, 2> > plane (
        makeFieldVector( constraint[0]->point() ),
        makeFieldVector( constraint[1]->point() )
      );

      // collect vertices corresponding to their side of the constraint
      std::vector<VertexHandle> insideVhs, outsideVhs;
      for( const FacetHandle& facet : hole )
      {
        for ( std::size_t i = 1; i <= 2; ++i )
        {
          const auto& vh = facet.first->vertex( (facet.second+i)%3 );

          if ( vh == constraint[0] || vh == constraint[1] )
            continue;

          if ( plane.signedDistance( makeFieldVector( vh->point() ) ) > 0 )
            insideVhs.push_back( vh );
          else
            outsideVhs.push_back( vh );
        }
      }

      // build outer map
      Vertex_pair_Facet_map outer_map;
      for (const auto& f : hole)
      {
        Vertex_pair vt = make_vertex_pair(f);
        outer_map[vt] = f;

        for (int i = 0; i < 3; i++)
          if ( i != f.second )
          {
            const auto& vh = f.first->vertex(i);
            vh->set_face(f.first);

            if ( vh != constraint[0] && vh != constraint[1] )
            {
              if ( plane.signedDistance( makeFieldVector( vh->point() ) ) > 0 )
                insideVhs.push_back( vh );
              else
                outsideVhs.push_back( vh );
            }
          }
      }

      // unify the points
      std::sort( insideVhs.begin(), insideVhs.end() );
      auto last = std::unique( insideVhs.begin(), insideVhs.end() );
      insideVhs.erase( last, insideVhs.end() );

      Vertex_pair_Facet_map outer_map2;
      // remesh the insideVhs
      {
        Vertex_pair_Facet_map inner_map;
        Vertex_handle_unique_hash_map vmap;

        ThisType del;
        // add all vertices on one side to a triangulation
        for ( const auto& vertex : insideVhs )
        {
          auto vh = del.insert(vertex->point());
          vmap[vh] = vertex;
        }

        // add the constraint vertices
        for ( const auto& vertex : constraint )
        {
          auto vh = del.insert(vertex->point());
          vmap[vh] = vertex;
        }

        CGAL_triangulation_assertion(del.dimension() == 2);

        for(auto it = del.finite_faces_begin(); it != del.finite_faces_end(); ++it)
          for(int i = 0; i < 3; i++)
          {
            FacetHandle f (it, i);
            Vertex_pair vp_aux = make_vertex_pair(f);
            Vertex_pair vp( vmap[ vp_aux.first ], vmap[ vp_aux.second ] );
            if ( vp.first > vp.second )
              std::swap( vp.first, vp.second );
            inner_map[vp] = f;
          }

        Vertex_pair vtc ( constraint[0], constraint[1] );
        if ( vtc.first > vtc.second )
          std::swap( vtc.first, vtc.second );

        while(!outer_map.empty())
        {
          auto oit = outer_map.begin();
          auto o_vt_f_pair = *oit;
          auto o_ch = o_vt_f_pair.second.first;
          unsigned int o_i = o_vt_f_pair.second.second;

          auto iit = inner_map.find( o_vt_f_pair.first );

          if(iit == inner_map.end())
          {
            outer_map2[ o_vt_f_pair.first ] = o_vt_f_pair.second;
            outer_map.erase(oit);
            continue;
          }

          auto i_vt_f_pair = *iit;
          auto i_ch = i_vt_f_pair.second.first;
          unsigned int i_i = i_vt_f_pair.second.second;

          //! check if we got the same face handle
          if ( vmap[i_ch->vertex(0)] == o_ch->vertex(0)
            && vmap[i_ch->vertex(1)] == o_ch->vertex(1)
            && vmap[i_ch->vertex(2)] == o_ch->vertex(2))
          {
            auto tmpch = i_ch;
            i_ch = tmpch->neighbor(i_i);
            i_i = i_ch->index(tmpch);
          }

          // create a new cell and glue it to the outer surface
          auto new_ch = create_face();
          *elements++ = new_ch;

          new_ch->set_vertices(
            vmap[i_ch->vertex(0)],
            vmap[i_ch->vertex(1)],
            vmap[i_ch->vertex(2)]
          );

          o_ch->set_neighbor(o_i, new_ch);
          new_ch->set_neighbor(i_i, o_ch);

          // for the other faces check, if they can also be glued
          for(int i = 0; i < 3; i++)
          {
            if(i != i_i)
            {
              FacetHandle f (new_ch, i);
              Vertex_pair vt = make_vertex_pair(f);
              auto oit2 = outer_map.find(vt);
              if(oit2 == outer_map.end())
              {
                if ( vt == vtc )
                {
                  outer_map2[vt] = f;
                }
                else
                {
                  outer_map[vt] = f;
                }
              }
              else
              {
                // glue the faces
                auto o_vt_f_pair2 = *oit2;
                auto o_ch2 = o_vt_f_pair2.second.first;
                int o_i2 = o_vt_f_pair2.second.second;
                o_ch2->set_neighbor(o_i2, new_ch);
                new_ch->set_neighbor(i, o_ch2);
                outer_map.erase(oit2);
              }
            }
          }
          outer_map.erase(oit);
        }
      }

      // unify the points
      std::sort( outsideVhs.begin(), outsideVhs.end() );
      last = std::unique( outsideVhs.begin(), outsideVhs.end() );
      outsideVhs.erase( last, outsideVhs.end() );

      // remesh the outsideVhs
      {
        Vertex_pair_Facet_map inner_map;
        Vertex_handle_unique_hash_map vmap;

        ThisType del;
        // add all vertices on one side to a triangulation
        for ( const auto& vertex : outsideVhs )
        {
          auto vh = del.insert(vertex->point());
          vmap[vh] = vertex;
        }

        // add the constraint vertices
        for ( const auto& vertex : constraint )
        {
          auto vh = del.insert(vertex->point());
          vmap[vh] = vertex;
        }

        CGAL_triangulation_assertion(del.dimension() == 2);

        for(auto it = del.finite_faces_begin(); it != del.finite_faces_end(); ++it)
          for(int i = 0; i < 3; i++)
          {
            FacetHandle f (it, i);
            Vertex_pair vp_aux = make_vertex_pair(f);
            Vertex_pair vp( vmap[ vp_aux.first ], vmap[ vp_aux.second ] );
            if ( vp.first > vp.second )
              std::swap( vp.first, vp.second );
            inner_map[vp] = f;
          }

        while(!outer_map2.empty())
        {
          auto oit = outer_map2.begin();
          auto o_vt_f_pair = *oit;
          auto o_ch = o_vt_f_pair.second.first;
          unsigned int o_i = o_vt_f_pair.second.second;

          auto iit = inner_map.find( o_vt_f_pair.first );

          if(iit == inner_map.end())
            DUNE_THROW( InvalidStateException, "Not found in inner map.");

          auto i_vt_f_pair = *iit;
          auto i_ch = i_vt_f_pair.second.first;
          unsigned int i_i = i_vt_f_pair.second.second;

          //! check if we got the same face handle
          if ( vmap[i_ch->vertex(0)] == o_ch->vertex(0)
            && vmap[i_ch->vertex(1)] == o_ch->vertex(1)
            && vmap[i_ch->vertex(2)] == o_ch->vertex(2))
          {
            auto tmpch = i_ch;
            i_ch = tmpch->neighbor(i_i);
            i_i = i_ch->index(tmpch);
          }

          // create a new cell and glue it to the outer surface
          auto new_ch = create_face();
          *elements++ = new_ch;

          new_ch->set_vertices(
            vmap[i_ch->vertex(0)],
            vmap[i_ch->vertex(1)],
            vmap[i_ch->vertex(2)]
          );

          o_ch->set_neighbor(o_i, new_ch);
          new_ch->set_neighbor(i_i, o_ch);

          // for the other faces check, if they can also be glued
          for(int i = 0; i < 3; i++)
          {
            if(i != i_i)
            {
              FacetHandle f (new_ch, i);
              Vertex_pair vt = make_vertex_pair(f);
              auto oit2 = outer_map2.find(vt);
              if(oit2 == outer_map2.end())
              {
                outer_map2[vt] = f;
              }
              else
              {
                // glue the faces
                auto o_vt_f_pair2 = *oit2;
                auto o_ch2 = o_vt_f_pair2.second.first;
                int o_i2 = o_vt_f_pair2.second.second;
                o_ch2->set_neighbor(o_i2, new_ch);
                new_ch->set_neighbor(i, o_ch2);
                outer_map2.erase(oit2);
              }
            }
          }
          outer_map2.erase(oit);
        }
      }

      delete_vertex(vh);
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
