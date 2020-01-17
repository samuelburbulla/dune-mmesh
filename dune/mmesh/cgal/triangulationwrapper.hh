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
  //! TriangulationWrapper<2>
  template<>
  class TriangulationWrapper<2> :
    public MMeshDefaults::Triangulation<2>::type
  {
      using ThisType = TriangulationWrapper<2>;
      using BaseType = typename MMeshDefaults::Triangulation<2>::type;
      using EdgeHandle = std::pair<typename BaseType::Face_handle, int>;

  public:
    template< class VertexHandle, class Elements >
    void removeAndGiveNewElements (const VertexHandle& vh, const Elements& elements)
    {
      // this->remove_and_give_new_faces(vh, cells);
      std::list<EdgeHandle> hole;
      this->make_hole(vh, hole);
      this->fill_hole_delaunay(hole, elements);
      this->delete_vertex(vh);
    }
  };

  template< class Triangulation >
  class VertexRemover
  {
  public:
    typedef CGAL::Nullptr_t Hidden_points_iterator;

    VertexRemover(Triangulation &tmp_) : tmp(tmp_) {}

    Triangulation &tmp;

    using Cell_handle = typename Triangulation::Cell_handle;
    using Point = typename Triangulation::Point;

    template< class T >
    void add_hidden_points(T) {}
    Hidden_points_iterator hidden_points_begin() { return NULL; }
    Hidden_points_iterator hidden_points_end() { return NULL; }

    CGAL::Bounded_side side_of_bounded_circle(const Point &p, const Point &q,
      const Point &r, const Point &s, bool perturb = false) const
    {
      assert(false);
      return CGAL::Bounded_side();
    }
  };

  //! TriangulationWrapper<3>
  template<>
  class TriangulationWrapper<3> :
    public MMeshDefaults::Triangulation<3>::type
  {
    using ThisType = TriangulationWrapper<3>;
    using BaseType = typename MMeshDefaults::Triangulation<3>::type;
    using Remover = VertexRemover<ThisType>;
    /*
    using Delaunay = typename MMeshDefaults::Delaunay<3>::type;
    using Remover = VertexRemover<Delaunay>;

    typedef CGAL::Triple<Vertex_handle, Vertex_handle, Vertex_handle> Vertex_triple;
    typedef boost::unordered_map<Vertex_triple, Facet> Vertex_triple_Facet_map;
    using Delaunay_cell_handle = typename Delaunay::Cell_handle;
    using Delaunay_vertex_handle = typename Delaunay::Vertex_handle;
    using Delaunay_facet = typename Delaunay::Facet;
    typedef CGAL::Unique_hash_map<Delaunay_vertex_handle, Vertex_handle, CGAL::Handle_hash_function> Vertex_handle_unique_hash_map;
    typedef CGAL::Triple<Delaunay_vertex_handle, Delaunay_vertex_handle, Delaunay_vertex_handle> Delaunay_vertex_triple;
    typedef boost::unordered_map<Vertex_triple, Delaunay_facet> Delaunay_vertex_triple_Facet_map;

    Vertex_triple make_vertex_triple(const Facet& f) const
    {
      auto ch = f.first;
      int i = f.second;

      return Vertex_triple(ch->vertex(vertex_triple_index(i,0)),
          ch->vertex(vertex_triple_index(i,1)),
          ch->vertex(vertex_triple_index(i,2)));
    }

    Delaunay_vertex_triple make_vertex_triple(const Delaunay_facet& f) const
    {
      auto ch = f.first;
      int i = f.second;

      return Delaunay_vertex_triple(ch->vertex(vertex_triple_index(i,0)),
          ch->vertex(vertex_triple_index(i,1)),
          ch->vertex(vertex_triple_index(i,2)));
    }

    void make_canonical(auto& t) const
    {
      int i = (t.first < t.second) ? 0 : 1;
      if(i==0) {
        i = (t.first < t.third) ? 0 : 2;
      } else {
        i = (t.second < t.third) ? 1 : 2;
      }
      auto tmp = t.first;
      switch(i){
      case 0: return;
      case 1:
        tmp = t.first;
        t.first = t.second;
        t.second = t.third;
        t.third = tmp;
        return;
      default:
        tmp = t.first;
        t.first = t.third;
        t.third = t.second;
        t.second = tmp;
      }
    }

    void make_hole_3D( Vertex_handle v, Vertex_triple_Facet_map& outer_map, std::vector<Cell_handle> & hole)
    {
      CGAL_triangulation_expensive_precondition( ! test_dim_down(v) );

      incident_cells(v, std::back_inserter(hole));

      for (typename std::vector<Cell_handle>::iterator cit = hole.begin(),
           end = hole.end(); cit != end; ++cit) {
        int indv = (*cit)->index(v);
        Cell_handle opp_cit = (*cit)->neighbor( indv );
        Facet f(opp_cit, opp_cit->index(*cit));
        Vertex_triple vt = make_vertex_triple(f);
        make_canonical(vt);
        outer_map[vt] = f;
        for (int i=0; i<4; i++)
          if ( i != indv )
            (*cit)->vertex(i)->set_cell(opp_cit);
      }
    }
  public:
    using Perturbation_order = typename BaseType::Perturbation_order;

    CGAL::Orientation coplanar_orientation(Point p, Point q, Point r)
    {
      return this->coplanar_orientation( p, q, r );
    }
    */

  public:
    template< class VertexHandle, class Elements >
    void removeAndGiveNewElements (const VertexHandle& v, Elements fit)
    {
      // [CGAL version]
      ThisType tmp;
      VertexRemover<ThisType> remover(tmp);
      remove_and_give_new_cells(v, remover, fit);
      return;

      /*
      // [Delaunay version]
      // Compare to 'remove_and_give_new_cells(vh, remover, elements)'
      // (taken from CGAL/Triangulation_3.hh:5178-5316)
      // We use a Delaunay_triangulation_3 in the hole and glue it
      // with a Triangulation_3 around.
      // -------------------------------------------------------------

      using namespace CGAL;
      CGAL_triangulation_precondition(dimension() == 3);

      std::vector<Cell_handle> hole;
      hole.reserve(64);

      // Construct the set of vertex triples on the boundary
      // with the facet just behind
      Vertex_triple_Facet_map outer_map;

      make_hole_3D(v, outer_map, hole);

      // collect all vertices on the boundary
      std::vector<Vertex_handle> vertices;
      vertices.reserve(64);

      adjacent_vertices(v, std::back_inserter(vertices));

      std::sort(vertices.begin(), vertices.end());

      remeshHole( v, outer_map, vertices, fit );

      tds().delete_vertex(v);
      tds().delete_cells(hole.begin(), hole.end());
      */
    }

    /*
    void remeshHole( const auto& v, auto& outer_map, const auto& vertices, auto fit )
    {
      // create a Delaunay triangulation of the points on the boundary
      // and make a map from the vertices in remover.tmp towards the vertices
      // in *this
      Delaunay tmp;
      Remover remover(tmp);
      Delaunay_vertex_triple_Facet_map inner_map;

      bool inf = false;
      unsigned int i;

      Vertex_handle_unique_hash_map vmap;
      Delaunay_cell_handle ch = Delaunay_cell_handle();
      for(i=0; i < vertices.size(); i++){
        if(! is_infinite(vertices[i])){
          Delaunay_vertex_handle vh = remover.tmp.insert(vertices[i]->point(), ch);
          ch = vh->cell();
          vmap[vh] = vertices[i];
        }else {
          inf = true;
        }
      }

      if(remover.tmp.dimension()==2){
        Delaunay_vertex_handle fake_inf = remover.tmp.insert(v->point());
        vmap[fake_inf] = infinite_vertex();
      } else {
        vmap[remover.tmp.infinite_vertex()] = infinite_vertex();
      }

      CGAL_triangulation_assertion(remover.tmp.dimension() == 3);

      // Construct the set of vertex triples of remover.tmp
      // We reorient the vertex triple so that it matches those from outer_map
      // Also note that we use the vertices of *this, not of remover.tmp
      if(inf){
        for(auto it = remover.tmp.all_cells_begin(),
              end = remover.tmp.all_cells_end(); it != end; ++it)
        {
          for(i=0; i < 4; i++){
            typename Delaunay::Facet f = std::pair<Delaunay_cell_handle,int>(it,i);
            Delaunay_vertex_triple vt_aux = make_vertex_triple(f);
            Vertex_triple vt(vmap[vt_aux.first],vmap[vt_aux.third],vmap[vt_aux.second]);
            make_canonical(vt);
            inner_map[vt] = f;
          }
        }
      } else {
        for(auto it = remover.tmp.finite_cells_begin(),
            end = remover.tmp.finite_cells_end(); it != end; ++it)
        {
          for(i=0; i < 4; i++){
            typename Delaunay::Facet f = std::pair<Delaunay_cell_handle,int>(it,i);
            Delaunay_vertex_triple vt_aux = make_vertex_triple(f);
            Vertex_triple vt(vmap[vt_aux.first],vmap[vt_aux.third],vmap[vt_aux.second]);
            make_canonical(vt);
            inner_map[vt] = f;
          }
        }
      }

      // Grow inside the hole, by extending the surface
      while(! outer_map.empty()){
        typename Vertex_triple_Facet_map::iterator oit = outer_map.begin();
        while(is_infinite(oit->first.first) ||
              is_infinite(oit->first.second) ||
              is_infinite(oit->first.third)){
          ++oit;
          // otherwise the lookup in the inner_map fails
          // because the infinite vertices are different
        }
        typename Vertex_triple_Facet_map::value_type o_vt_f_pair = *oit;
        Cell_handle o_ch = o_vt_f_pair.second.first;
        unsigned int o_i = o_vt_f_pair.second.second;

        typename Delaunay_vertex_triple_Facet_map::iterator iit = inner_map.find( o_vt_f_pair.first );

        if(iit == inner_map.end())
          DUNE_THROW( GridError, "During remeshing an outer facet could not be matched to an inner facet!" );

        typename Delaunay_vertex_triple_Facet_map::value_type i_vt_f_pair = *iit;
        Delaunay_cell_handle i_ch = i_vt_f_pair.second.first;
        unsigned int i_i = i_vt_f_pair.second.second;

        // create a new cell and glue it to the outer surface
        Cell_handle new_ch = tds().create_cell();
        *fit++ = new_ch;

        new_ch->set_vertices(vmap[i_ch->vertex(0)], vmap[i_ch->vertex(1)],
                             vmap[i_ch->vertex(2)], vmap[i_ch->vertex(3)]);

        o_ch->set_neighbor(o_i,new_ch);
        new_ch->set_neighbor(i_i, o_ch);

        // for the other faces check, if they can also be glued
        for(i = 0; i < 4; i++){
          if(i != i_i){
            Facet f = std::pair<Cell_handle,int>(new_ch,i);
            Vertex_triple vt = make_vertex_triple(f);
            make_canonical(vt);
            std::swap(vt.second,vt.third);
            typename Vertex_triple_Facet_map::iterator oit2 = outer_map.find(vt);
            if(oit2 == outer_map.end()){
              std::swap(vt.second,vt.third);
              outer_map[vt]= f;
            } else {
              // glue the faces
              typename Vertex_triple_Facet_map::value_type o_vt_f_pair2 = *oit2;
              Cell_handle o_ch2 = o_vt_f_pair2.second.first;
              int o_i2 = o_vt_f_pair2.second.second;
              o_ch2->set_neighbor(o_i2,new_ch);
              new_ch->set_neighbor(i, o_ch2);
              outer_map.erase(oit2);
            }
          }
        }
        outer_map.erase(oit);
      }
    }
    */

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
