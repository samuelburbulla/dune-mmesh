// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_MMESHPOINTFIELDVECTOR_HH
#define DUNE_GRID_MMESHPOINTFIELDVECTOR_HH

// Dune includes
#include <dune/common/fvector.hh>

/** \file
 * \brief Helpers for conversion from CGAL::Point_x to DUNE::FieldVector
 */

namespace Dune
{

  /** \brief Helper function to create DUNE::FieldVector from CGAL::Point_2
   *
   * \tparam Kernel     the kernel of CGAL::Point_2
   */
  template< class Kernel >
  static inline FieldVector<typename Kernel::RT, 2> makeFieldVector( const CGAL::Point_2<Kernel>& p )
  {
      return FieldVector<typename Kernel::RT, 2>( { p.x(), p.y() } );
  }

  /** \brief Helper function to create DUNE::FieldVector from CGAL::Vector_2
   *
   * \tparam Kernel     the kernel of CGAL::Vector_2
   */
  template< class Kernel >
  static inline FieldVector<typename Kernel::RT, 2> makeFieldVector( const CGAL::Vector_2<Kernel>& p )
  {
      return FieldVector<typename Kernel::RT, 2>( { p.x(), p.y() } );
  }

  /** \brief Helper function to create DUNE::FieldVector from CGAL::Point_3
   *
   * \tparam Kernel     the kernel of CGAL::Point_3
   */
  template< class Kernel >
  static inline FieldVector<typename Kernel::RT, 3> makeFieldVector( const CGAL::Point_3<Kernel>& p )
  {
      return FieldVector<typename Kernel::RT, 3>( { p.x(), p.y(), p.z() } );
  }

  /** \brief Helper function to create DUNE::FieldVector from CGAL::Vector_3
   *
   * \tparam Kernel     the kernel of CGAL::Vector_3
   */
  template< class Kernel >
  static inline FieldVector<typename Kernel::RT, 3> makeFieldVector( const CGAL::Vector_3<Kernel>& p )
  {
      return FieldVector<typename Kernel::RT, 3>( { p.x(), p.y(), p.z() } );
  }


  /** \brief Convert FieldVector to CGAL Point
   */
  typedef CGAL::Exact_predicates_inexact_constructions_kernel PointKernel;

  /** \brief Convert FieldVector to CGAL Point
   *  \ingroup 2D
   */
  template< typename ctype >
  static inline auto makePoint( const Dune::FieldVector<ctype, 2>& v )
  {
    return CGAL::Point_2<PointKernel> ( v[ 0 ], v[ 1 ] );
  }

  /** \brief Convert FieldVector to CGAL Point
   *  \ingroup 3D
   */
  template< typename ctype >
  static inline auto makePoint( const Dune::FieldVector<ctype, 3>& v )
  {
    return CGAL::Point_3<PointKernel> ( v[ 0 ], v[ 1 ], v[ 2 ] );
  }

  /** \brief Compute circumcenter
   */
  template <class Geometry>
  auto computeCircumcenter(const Geometry& geo)
  {
    static constexpr int mydim = Geometry::mydim;

    if constexpr (mydim == 3)
    {
      // obtain circumcenter by CGAL
      return makeFieldVector(
                             CGAL::circumcenter(
                                                makePoint(geo.corner(0)),
                                                makePoint(geo.corner(1)),
                                                makePoint(geo.corner(2)),
                                                makePoint(geo.corner(3))
                                                )
                             );
    }

    if constexpr (mydim == 2)
    {
      // obtain circumcenter by CGAL
      return makeFieldVector(
                             CGAL::circumcenter(
                                                makePoint(geo.corner(0)),
                                                makePoint(geo.corner(1)),
                                                makePoint(geo.corner(2))
                                                )
                             );
    }

    if constexpr (mydim == 1)
      return 0.5 * (geo.corner(0) + geo.corner(1));

    if constexpr (mydim == 0)
      return geo.corner(0);
  }

}  // namespace Dune

#endif
