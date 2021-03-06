// Copyright (c) 2018  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL: https://github.com/CGAL/cgal/blob/v5.2.2/Intersections_3/include/CGAL/Intersections_3/Point_3_Ray_3.h $
// $Id: Point_3_Ray_3.h 0f3305f 2020-01-16T17:20:13+01:00 Maxime Gimeno
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Maxime Gimeno

#ifndef CGAL_INTERSECTIONS_3_POINT_3_RAY_3_H
#define CGAL_INTERSECTIONS_3_POINT_3_RAY_3_H

#include <CGAL/Ray_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/Intersection_traits_3.h>

namespace CGAL {

namespace Intersections {

namespace internal {

template <class K>
inline
bool
do_intersect(const typename K::Point_3 &pt,
             const typename K::Ray_3 &ray,
             const K& k)
{
  return k.has_on_3_object()(ray,pt);
}


template <class K>
inline
bool
do_intersect(const typename K::Ray_3 &ray,
             const typename K::Point_3 &pt,
             const K& k)
{
  return k.has_on_3_object()(ray,pt);
}


template <class K>
typename CGAL::Intersection_traits
<K, typename K::Point_3, typename K::Ray_3>::result_type
intersection(const typename K::Point_3 &pt,
             const typename K::Ray_3 &ray,
             const K& k)
{
  if (do_intersect(pt,ray, k)) {
    return intersection_return<typename K::Intersect_3, typename K::Point_3, typename K::Ray_3>(pt);
  }
  return intersection_return<typename K::Intersect_3, typename K::Point_3, typename K::Ray_3>();
}

template <class K>
typename CGAL::Intersection_traits
<K, typename K::Ray_3, typename K::Point_3>::result_type
intersection(const typename K::Ray_3 &ray,
             const typename K::Point_3 &pt,
             const K& k)
{
  return internal::intersection(pt, ray, k);
}

} // namespace internal
} // namespace Intersections

CGAL_INTERSECTION_FUNCTION(Point_3, Ray_3, 3)
CGAL_DO_INTERSECT_FUNCTION(Point_3, Ray_3, 3)


} //namespace CGAL
#endif // CGAL_INTERSECTIONS_3_POINT_3_RAY_3_H
