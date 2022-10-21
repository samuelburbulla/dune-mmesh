#include "config.h"

#include <iostream>
#include <cmath>
#include <string>

#include <dune/common/exceptions.hh>

#include <dune/grid/common/mcmgmapper.hh>

#include <dune/mmesh/mmesh.hh>
#include <dune/mmesh/interface/curvatureoperator.hh>


constexpr int dim = GRIDDIM;
using Vector = Dune::FieldVector<double, dim>;

double exactCurvature (Vector pos)
{
  if constexpr (dim == 2)
  {
    //ellipse with semiaxis a = 0.3 and b = 0.15
    double a4 = 8.1e-3;
    double b4 = 5.0625e-4;
    return (a4 * b4) / pow(a4 * pos[1] * pos[1] + b4 * pos[0] * pos[0], 1.5);
  }
  else
  {
    //sphere with radius r = 0.15
    double curvature = 1.0 / 0.15;
    return curvature;
  }
}

Vector exactCenterOfCurvature (Vector pos, double a = 0.3, double b = 0.15)
{
  if constexpr (dim == 2)
  {
    //ellipse with semiaxis a = 0.3 and b = 0.15
    double a2 = 9e-2;
    double b2 = 2.25e-2;
    double e2 = a2 - b2;

    Vector center;
    center[0] = e2 * pow(pos[0], 3) / (a2 * a2);
    center[1] = -e2 * pow(pos[1], 3) / (b2 * b2);

    return center;
  }
  else
  {
    DUNE_THROW(Dune::Exception, "Not implemented");
  }
}

template< class Lengths, class GridFiles, class NumOfIElems, class ErrorsCurvature, class ErrorsCenter >
void output(const Lengths& characteristicLengths,
            const GridFiles& gridfiles,
            const NumOfIElems& numOfIElems,
            const ErrorsCurvature& errorsCurvature,
            const ErrorsCenter& errorsCenter)
{
  int numOfGridFiles = numOfIElems.size();

  for (int i = 0; i < numOfGridFiles; i++)
  {
    std::cout << "characteristic length: " << gridfiles[i] <<
      "\t#interface elements: " << numOfIElems[i] << " \terror curvature: "
      << errorsCurvature[i];
    if constexpr (dim == 2)
    {
      std::cout << "\terror center: " << errorsCenter[i];
    }
    std::cout << std::endl;
  }

  std::cout << "\nEOC curvature: " << std::endl;
  for (int i = 1; i < numOfGridFiles; i++)
  {
    const double EOC = log(errorsCurvature[i-1] / errorsCurvature[i])
    / log(characteristicLengths[i-1] / characteristicLengths[i]);
    std::cout << "EOC: " << EOC << "\t(interface elements " <<
      numOfIElems[i-1] << " -> " << numOfIElems[i] << ")\n";

    if constexpr (dim == 2)
    {
      assert(EOC > 1.6);
    }
  }

  if constexpr (dim == 2)
  {
    std::cout << "\nEOC center of curvature: " << std::endl;
    for (int i = 1; i < numOfGridFiles; i++)
    {
      const double EOC = log(errorsCenter[i-1] / errorsCenter[i])
      / log(characteristicLengths[i-1] / characteristicLengths[i]);
      std::cout << "EOC: " << EOC << "\t(interface elements " <<
        numOfIElems[i-1] << " -> " << numOfIElems[i] << ")\n";
      assert(EOC > 1.6);
    }
  }
}

int main(int argc, char** argv)
{
  try
  {
    Dune::MPIHelper::instance(argc, argv);
    using Grid = Dune::MovingMesh<dim>;
    using GridFactory = Dune::GmshGridFactory<Grid, /*implicit=*/false>;
    using IGridView = typename Grid::InterfaceGrid::LeafGridView;
    using IMapper =
      typename Dune::MultipleCodimMultipleGeomTypeMapper<IGridView>;
    using CurvatureOperatorElement =
      Dune::CurvatureOperator<IGridView, IMapper, Dune::CurvatureLayout::Element>;
    using CurvatureOperatorVertex =
      Dune::CurvatureOperator<IGridView, IMapper, Dune::CurvatureLayout::Vertex>;

    const std::string gridType = (dim == 2) ? "ellipse" : "sphere";

    //get interface grid view
    GridFactory gridFactory("grids/" + gridType + std::to_string(dim) + "d.msh");
    const Grid& grid = *gridFactory.grid();
    const IGridView& iGridView = grid.interfaceGrid().leafGridView();

    //define element mapper
    const IMapper iElemMapper(iGridView, Dune::mcmgElementLayout());
    const IMapper iVertexMapper(iGridView, Dune::mcmgVertexLayout());

    //storage for errors
    double errorCurvatureElement = 0;
    double errorCenterElement = 0;
    double errorCurvatureVertex = 0;
    double errorCenterVertex = 0;

    //storage curvatures on the interface grid
    std::vector<double> curvaturesElement(iElemMapper.size());
    std::vector<Vector> centersElement(iElemMapper.size());
    std::vector<double> exactCurvaturesElement(iElemMapper.size());
    std::vector<double> curvaturesVertex(iVertexMapper.size());
    std::vector<Vector> centersVertex(iVertexMapper.size());
    std::vector<double> exactCurvaturesVertex(iVertexMapper.size());

    //determine curvatures and sphere center points
    CurvatureOperatorElement curvOpElem(iGridView, iElemMapper);
    curvOpElem(curvaturesElement, centersElement);

    CurvatureOperatorVertex curvOpVertex(iGridView, iVertexMapper);
    curvOpVertex(curvaturesVertex, centersVertex);

    //calculate L2 error for element layout
    for (const auto& iElem: elements(iGridView))
    {
      const int iElemIdx = iElemMapper.index(iElem);
      const auto iGeo = iElem.geometry();
      const auto& center = iGeo.center();
      const double iVolume = iGeo.volume();
      double exCurv = exactCurvature(center);

      errorCurvatureElement += pow(curvaturesElement[iElemIdx] - exCurv, 2)*iVolume;

      if constexpr (dim == 2)
        errorCenterElement += (centersElement[iElemIdx] - exactCenterOfCurvature(center)).two_norm2()*iVolume;

      exactCurvaturesElement[iElemIdx] = exCurv;
    }

    //calculate maximum error for vertex layout
    for (const auto& vertex : vertices(iGridView))
    {
      const int vertexIdx = iVertexMapper.index(vertex);
      const auto& center = vertex.geometry().center();

      double exCurv =  exactCurvature(center);

      errorCurvatureVertex = std::max(errorCurvatureVertex, std::abs(curvaturesVertex[vertexIdx] - exCurv));

      if constexpr (dim == 2)
        errorCenterVertex = std::max(errorCenterVertex,
          (centersVertex[vertexIdx] - exactCenterOfCurvature(center)).infinity_norm()
        );

      exactCurvaturesVertex[vertexIdx] = exCurv;
    }

    Dune::VTKWriter<IGridView> vtkWriter(iGridView);
    vtkWriter.addCellData(curvaturesElement, "curvatureElement");
    vtkWriter.addCellData(exactCurvaturesElement, "exactCurvatureElement");
    vtkWriter.addVertexData(curvaturesVertex, "curvatureVertex");
    vtkWriter.addVertexData(exactCurvaturesVertex, "exactCurvatureVertex");
    vtkWriter.write("curvature");

    std::cout << "ErrorCurvatureElement: " << errorCurvatureElement << std::endl;
    std::cout << "ErrorCenterElement: " << errorCenterElement << std::endl;
    std::cout << "ErrorCurvatureVertex: " << errorCurvatureVertex << std::endl;
    std::cout << "ErrorCenterVertex: " << errorCenterVertex << std::endl;

    assert( errorCurvatureElement < 1.0 );
    assert( errorCenterElement < 1e-3 );
    assert( errorCurvatureVertex < 1.0 );
    assert( errorCenterVertex < 1e-2 );
  }
  catch (Dune::Exception& e)
  {
    std::cerr << "Dune reported error: " << e.what() << std::endl;
  }
  catch (std::exception& e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
  }
  catch (...)
  {
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
