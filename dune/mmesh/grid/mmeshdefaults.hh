// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_GRID_MMESHDEFAULTS_HH
#define DUNE_MMESH_GRID_MMESHDEFAULTS_HH

#include <dune/mmesh/mmesh.hh>

#define CGAL_NO_POSTCONDITIONS

// CGAL includes
#include <CGAL/exceptions.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// 2D
#include <CGAL/Triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_hierarchy_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_2.h>

// 3D
#include <CGAL/Triangulation_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>

/** \file
 * \brief Some defaults to generate common CGAL triangulations
 */

namespace Dune
{
  namespace MMeshDefaults
  {

    template< int dim >
    class Delaunay;

    /*!
     * \brief A Delaunay triangulation in 2D
     */
    template<>
    class Delaunay<2>
    {
    public:
      struct FaceInfo {
        std::size_t id;
        bool idWasSet = false;
        std::size_t index;
      };

      struct VertexInfo {
        std::size_t id;
        bool idWasSet = false;
        std::size_t index;
        std::size_t boundaryid = 0;
        std::size_t insertionLevel = 0;
      };

    private:
      typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

      typedef CGAL::Delaunay_mesh_vertex_base_2<K> Vbbb;
      typedef CGAL::Triangulation_vertex_base_with_info_2<VertexInfo, K, Vbbb> Vbb;
      typedef CGAL::Triangulation_hierarchy_vertex_base_2<Vbb> Vb;

      typedef CGAL::Delaunay_mesh_face_base_2<K> Fbbb;
      typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo, K, Fbbb> Fbb;
      typedef CGAL::Triangulation_face_base_2<K,Fbb> Fb;

      typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;

    public:
      typedef CGAL::Delaunay_triangulation_2<K, Tds> Triangulation;
    };


    /*!
     * \brief A Delaunay triangulation in 3D
     */
    template<>
    class Delaunay<3>
    {
    public:
      struct CellInfo {
        std::size_t id;
        bool idWasSet = false;
        std::size_t index;
      };

      struct VertexInfo {
        std::size_t id;
        bool idWasSet = false;
        std::size_t index;
        std::size_t boundaryid = 0;
      };

    private:
      typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

      typedef CGAL::Triangulation_vertex_base_with_info_3<VertexInfo, K> Vb;

      typedef CGAL::Delaunay_triangulation_cell_base_3<K> Fbb;
      typedef CGAL::Triangulation_cell_base_with_info_3<CellInfo, K, Fbb> Fb;

      typedef CGAL::Triangulation_data_structure_3<Vb, Fb> Tds;

    public:
      typedef CGAL::Delaunay_triangulation_3<K, Tds> Triangulation;
    };

  } // end namespace MMeshDefaults

} // end namespace Dune

#endif
