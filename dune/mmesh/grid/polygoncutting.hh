#ifndef DUNE_MMESH_GRID_POLYGONCUTTING_HH
#define DUNE_MMESH_GRID_POLYGONCUTTING_HH

#include <iostream>
#include <algorithm>

namespace Dune
{
  template<class Scalar, class Point>
  class PolygonCutting
  {
  public:

    //implements the Sutherland-Hodgman algorithm to determine the polygon that
    //emerges from clipping a subject polygon with respect to the edges of a
    //clip polygon, subject polygon and clip polygon as well as the resulting
    //polygon are represented by lists of their corners points with the points
    //arranged in a consecutive counterclockwise order
    template< typename Polygon >
    static std::vector<Point> sutherlandHodgman(const Polygon& subjectPolygon,
      const Polygon& clipPolygonInput)
    {
      std::vector<Point>
        outputPolygon(subjectPolygon.begin(), subjectPolygon.end());
      std::vector<Point>
        clipPolygon(clipPolygonInput.begin(), clipPolygonInput.end());

      const int clipSize = clipPolygon.size();

      //iterate over the edges of the clipping polygon represented by
      //consecutive points (indices clipIdx and clipIdxNext) in clipPolygon
      for (int clipIdx = 0; clipIdx < clipSize; clipIdx++)
      {
        const int clipIdxNext = (clipIdx + 1) % clipSize;

        std::vector<Point> inputPolygon(outputPolygon);
        outputPolygon.clear();
        const int inputSize = inputPolygon.size();

        //iterate over the edges of the subject polygon represented by
        //consectutive points (indices subjIdx and subjIdxNext) and clip each
        //subject edge with respect to current clipping edge (indices clipIdx
        //and clipIdxNext)
        for (int inputIdx = 0; inputIdx < inputSize; inputIdx++)
        {
          const int inputIdxNext = (inputIdx + 1) % inputSize;

          //Case differentation: Determine the positions of the corner points
          //of the subject edge p1, p2 with respect to the current clipping
          //edge given as the line from q1 to q2.
          //This is done by evaluating the cross product
          //d := (q1 - q2) x (p_i - q_1) = [(q1 - q2) x p_i] + q2 x q1
          //for i = 1, 2.
          //d < 0 => p_i inside, d > 0 => p_i outside,
          //d = 0 => p_i on clipping edge

          Point deltaQ = clipPolygon[clipIdxNext] - clipPolygon[clipIdx];

          Scalar q2_cross_q1 =
            clipPolygon[clipIdx][1] * clipPolygon[clipIdxNext][0] -
            clipPolygon[clipIdx][0] * clipPolygon[clipIdxNext][1];

          Scalar firstPos = - deltaQ[0] * inputPolygon[inputIdx][1]
            + deltaQ[1] * inputPolygon[inputIdx][0] + q2_cross_q1;
          Scalar secondPos = - deltaQ[0] * inputPolygon[inputIdxNext][1]
            + deltaQ[1] * inputPolygon[inputIdxNext][0] + q2_cross_q1;

          static constexpr double EPSILON = 1e-14;
          //case 1: first point outside, second point inside
          if ( firstPos > EPSILON && secondPos < -EPSILON )
          {
            //add intersection point and second point to output polygon
            outputPolygon.push_back( lineIntersectionPoint(
              inputPolygon[inputIdx], inputPolygon[inputIdxNext],
              clipPolygon[clipIdx], clipPolygon[clipIdxNext] )
            );
            outputPolygon.push_back( inputPolygon[inputIdxNext] );
          }

          //case 2: first point inside, second point outside
          else if ( firstPos < -EPSILON && secondPos > EPSILON )
          {
            //add intersection point to output polygon
            outputPolygon.push_back( lineIntersectionPoint(
              inputPolygon[inputIdx], inputPolygon[inputIdxNext],
              clipPolygon[clipIdx], clipPolygon[clipIdxNext] )
            );
          }

          //case 3: first point outside or on the edge, scond point outside
          else if ( firstPos >= -EPSILON  &&  secondPos > EPSILON )
          {
            //do nothing
            continue;
          }

          //case 4: first point inside, second point inside or on the edge
          else
          {
            //add second point to output polygon
            outputPolygon.push_back( inputPolygon[inputIdxNext] );
          }
        }
      }

      return outputPolygon;
    }

    //computes the point of intersection in 2d between the line given by the
    //points p1, p2 and the line defined by the points q1, q2
    static Point lineIntersectionPoint(const Point& p1, const Point& p2,
      const Point& q1, const Point& q2)
    {
      Point deltaP = p1 - p2;
      Point deltaQ = q1 - q2;

      Scalar scaleFactor = deltaP[0] * deltaQ[1] - deltaQ[0] * deltaP[1];

      deltaP *= q1[0] * q2[1] - q1[1] * q2[0];
      deltaQ *= p1[0] * p2[1] - p1[1] * p2[0];

      Point intersection = deltaQ - deltaP;
      intersection /= scaleFactor;

      return intersection;
    }


    //computes the area of a convex polygon by evaluating the shoelace formula,
    //the polygon is represented by a corner point list with the points arranged
    //consecutively in counterclockwise order (for clockwise order (-1)*area
    //will be returned)
    template< typename Points >
    static Scalar polygonArea (const Points& polygonPoints)
    {
      const int polygonSize = polygonPoints.size();

      if (polygonSize < 3)
      {
        return 0.0;
      }

      Scalar area(0.0);

      for (int i = 0; i < polygonSize; i++)
      {
        const int j = (i+1) % polygonSize;
        area += polygonPoints[i][0] * polygonPoints[j][1];
        area -= polygonPoints[j][0] * polygonPoints[i][1];
      }

      return 0.5 * area;
    }

#if CGAL_INTERSECTION
    template< typename Triang >
    static Scalar intersectionVolumeCGAL( const Triang& subjectTriangle,
      const Triang& clipTriangle )
    {
      using Epeck = CGAL::Exact_predicates_exact_constructions_kernel;
      using PointCGAL = CGAL::Point_2<Epeck>;
      using Triangle = CGAL::Triangle_2<Epeck>;

      std::array<PointCGAL, 3> v;
      for ( int i = 0; i < 3; ++i )
      {
        const auto& p = makePoint( subjectTriangle[i] );
        v[i] = PointCGAL( p.x(), p.y() );
      }
      const Triangle triangle1 ( v[0], v[1], v[2] );

      for ( int i = 0; i < 3; ++i )
      {
        const auto& p = makePoint( clipTriangle[i] );
        v[i] = PointCGAL( p.x(), p.y() );
      }
      const Triangle triangle2 ( v[0], v[1], v[2] );

      const auto result = CGAL::intersection( triangle1, triangle2 );

      using PointList = std::vector<PointCGAL>;
      using Polygon = CGAL::Polygon_2<Epeck>;
      if (result)
      {
        // intersection is triangle
        if (const Triangle* t = boost::get<Triangle>(&*result))
          return std::abs( CGAL::to_double(
            CGAL::area( t->vertex(0), t->vertex(1), t->vertex(2) ) ) );

        // intersection is polygon
        else if ( const PointList* points = boost::get<PointList>(&*result))
        {
          Polygon polygon( points->begin(), points->end() );
          return std::abs( CGAL::to_double( polygon.area() ) );
        }
      }

      // there is no intersection
      return 0;
    }
#endif

  };
}

#endif
