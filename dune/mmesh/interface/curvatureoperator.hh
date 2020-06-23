// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*!
 * \file
 * \brief   Class defining an operator that assigns a curvature to each element
 *          or each vertex of the interface grid.
 */

#ifndef DUNE_MMESH_INTERFACE_CURVATUREOPERATOR_HH
#define DUNE_MMESH_INTERFACE_CURVATUREOPERATOR_HH

#include <algorithm>
#include <dune/common/dynvector.hh>
#include <dune/common/dynmatrix.hh>

namespace Dune
{

/*!
 * \brief   used as template parameter that determines whether the curvature is
 *          approximated at the elemets or the vertices of the grid
 */
enum CurvatureLayout
{
  Vertex,
  Element
};

/*!
 * \brief   Class defining an operator that assigns a curvature to each element
 *          or each vertex of the interface grid. The curvature is approximated
 *          by interpolation with a sphere (through the center points of an
 *          interface element and its neighbors). The enum CurvatureLayout
 *          determines whether the curvature is approximated at the vertices
 *          or elements of the grid.
 *          NOTE that the algorithm IN 3D might provide poor approximations for
 *          non-spherical interfaces and IN 3D the algorithm in general does not
 *          converge under grid refinement
 */
template<class IGridView, class IGridMapper, CurvatureLayout CL = Vertex>
class CurvatureOperator
{
public:
  static constexpr int dim = IGridView::dimensionworld;
  using Scalar = typename IGridView::ctype;
  using Vector = Dune::FieldVector<Scalar, dim>;

  /*!
   * \brief Constructor.
   *
   * \param iGridView The grid view of the interface grid
   * \param iGridMapper Element or Vertex mapper for the interface grid
   *        depending on the specified curvature layout
   */
  CurvatureOperator(const IGridView& iGridView, const IGridMapper& iGridMapper)
   : iGridView_(iGridView), iGridMapper_(iGridMapper) {};

  /*!
   * \brief Operator that assigns a curvature to each element or each vertex
   *        of the interface grid. The curvature is approximated by
   *        interpolation with a sphere (osculating circle/sphere through
   *        incident interface vertices).
   *
   * \param curvatures A container to store the curvature of each element or
   *        each vertex of the interface grid
   * \param centers A container to store the centers of curvatures (sphere
   *        center points) corresponding to the approximated curvature of each
   *        element or each vertex of the interface grid
   */
  template<class Curvatures, class Centers, CurvatureLayout Layout = CL>
  std::enable_if_t<Layout == Vertex, void>
  operator() (Curvatures& curvatures, Centers& centers) const
  {
    //storage for incident interface vertices
    std::vector<Vector> vertexPoints;
    Vector center;

    for (const auto& vertex : vertices(iGridView_))
    {
      const int iVertexIdx = iGridMapper_.index(vertex);
      vertexPoints.clear();
      vertexPoints.push_back(vertex.geometry().center());

      for (const auto& incidentVertex : incidentInterfaceVertices(vertex))
      {
        vertexPoints.push_back(incidentVertex.geometry().center());
      }

      //determine curvature associated with iVertexMapper (approximated by a
      //sphere with center "center")
      curvatures[iVertexIdx] = 1.0 / getRadius(vertexPoints, center);
      centers[iVertexIdx] = center;
    }
  }

  template<class Curvatures, class Centers, CurvatureLayout Layout = CL>
  std::enable_if_t<Layout == Element, void>
  operator() (Curvatures& curvatures, Centers& centers) const
  {
   //storage for vertices of incident interface elements
   std::vector<Vector> vertexPoints;
   Vector center;

   for (const auto& iElem : elements(iGridView_))
   {
     int iElemIdx = iGridMapper_.index(iElem);
     vertexPoints.clear();

     for (const auto& intersct : intersections(iGridView_, iElem))
     {
       const auto& neighborGeo = intersct.outside().geometry();

       for (int i = 0; i < neighborGeo.corners(); i++)
       {
         if (std::find(vertexPoints.begin(), vertexPoints.end(),
           neighborGeo.corner(i)) == vertexPoints.end())
         {
           vertexPoints.push_back(neighborGeo.corner(i));
         }
       }
     }

     //determine curvature associated with iElem (approximated by a sphere
     //with center "center")
     curvatures[iElemIdx] = 1.0 / getRadius(vertexPoints, center);
     centers[iElemIdx] = center;
   }
 }

private:
  /*!
   * \brief Determines the radius and center of the sphere that is uniquely
   *        defined by dim+1 points
   *
   * \retrun The radius of the sphere (or infinity if all points lie on a
   *         hyperplane)
   *
   * \param points Array with dim+1 points on the sphere
   * \param center storage for the center of the sphere
   */
  template< class Points >
  double getRadius (const Points& points, Vector& center) const
  {
    using LargeVector = Dune::FieldVector<Scalar, dim+1>;
    using Matrix = Dune::FieldMatrix<Scalar, dim+1, dim+1>;
    using LargeDynVector = Dune::DynamicVector<Scalar>;
    using DynMatrix = Dune::DynamicMatrix<Scalar>;

    LargeDynVector b(points.size(), 0.0);
    DynMatrix A(points.size(), dim+1, 0.0);
    LargeVector x(0.0);
    LargeVector ATb(0.0); //A.transpose()*b
    Matrix ATA(0.0); //A.transpose()*A

    for (unsigned int i = 0; i < points.size(); i++)
    {
      b[i] = -points[i].two_norm2();
      A[i][0] = 1.0;

      for (int j = 0; j < dim; j++)
      {
        A[i][j+1] = points[i][j];
      }
    }

    A.mtv(b, ATb);

    for (int i = 0; i <= dim; i++)
      for (int j = 0; j <= dim; j++)
        for (std::size_t k = 0; k < points.size(); k ++)
          ATA[i][j] += A[k][i] * A[k][j];

    if (std::abs( ATA.determinant() ) > epsilon)
    {
      ATA.solve(x, ATb);

      double r = 0.0;
      for (int j = 1; j <= dim; j++)
      {
        r += x[j]*x[j];
        center[j-1] = -0.5*x[j];
      }

      r = sqrt(std::abs(0.25*r - x[0]));

      return r;
    }

    return std::numeric_limits<double>::infinity();
  }

  //The grid view of the interface grid
  const IGridView& iGridView_;
  //Element mapper for the interface grid
  const IGridMapper& iGridMapper_;

  //floating point accuracy
  static constexpr double epsilon = std::numeric_limits<double>::epsilon();
};

} // end namespace Dune

#endif
