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
    using Grid = Dune::MovingMesh<dim>;
    using GridFactory = Dune::GmshGridFactory<Grid, /*implicit=*/false>;
    using IGridView = typename Grid::InterfaceGrid::LeafGridView;
    using IMapper =
      typename Dune::MultipleCodimMultipleGeomTypeMapper<IGridView>;
    using CurvatureOperatorElement =
      Dune::CurvatureOperator<IGridView, IMapper, Dune::CurvatureLayout::Element>;
    using CurvatureOperatorVertex =
      Dune::CurvatureOperator<IGridView, IMapper, Dune::CurvatureLayout::Vertex>;
    using StringList = std::vector<std::string>;
    using NumList = std::vector<double>;

    constexpr int numOfGridFiles = (dim == 2) ? 9 : 3;
    const StringList gridfiles = (dim == 2) ? StringList({"5e-2", "3e-2",
      "2e-2", "1e-2", "9e-3", "8e-3", "7e-3", "6e-3", "5e-3"})
      : StringList({"9e-2", "7e-2", "5e-2"});
    const NumList characteristicLengths = (dim == 2)
      ? NumList({5e-2, 3e-2, 2e-2, 1e-2, 9e-3, 8e-3, 7e-3, 6e-3, 5e-3})
      :  NumList({9e-2, 7e-2, 5e-2});

    std::array<double, numOfGridFiles> errorsCurvatureElement;
    std::array<double, numOfGridFiles> errorsCenterElement;
    std::array<double, numOfGridFiles> errorsCurvatureVertex;
    std::array<double, numOfGridFiles> errorsCenterVertex;
    std::array<double, numOfGridFiles> numOfIElems;

    for (int i = 0; i < numOfGridFiles; i++)
    {
      errorsCurvatureElement[i] = 0.0;
      errorsCenterElement[i] = 0.0;
      errorsCurvatureVertex[i] = 0.0;
      errorsCenterVertex[i] = 0.0;

      const std::string gridType = (dim == 2) ? "ellipse" : "sphere";

      //get interface grid view
      GridFactory gridFactory("grids/" + gridType + std::to_string(dim) + "d_"
        + gridfiles[i] + ".msh");
      const Grid& grid = *gridFactory.grid();
      const IGridView& iGridView = grid.interfaceGrid().leafGridView();

      //get number of interface segments
      numOfIElems[i] = iGridView.size(0);

      //define element mapper
      const IMapper iElemMapper(iGridView, Dune::mcmgElementLayout());
      const IMapper iVertexMapper(iGridView, Dune::mcmgVertexLayout());

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

        errorsCurvatureElement[i] +=
          pow(curvaturesElement[iElemIdx] - exCurv, 2)*iVolume;

        if constexpr (dim == 2)
        {
          errorsCenterElement[i] += (centersElement[iElemIdx]
            - exactCenterOfCurvature(center)).two_norm2()*iVolume;
        }

        exactCurvaturesElement[iElemIdx] = exCurv;
      }

      //calculate maximum error for vertex layout
      for (const auto& vertex : vertices(iGridView))
      {
        const int vertexIdx = iVertexMapper.index(vertex);
        const auto& center = vertex.geometry().center();

        double exCurv =  exactCurvature(center);

        errorsCurvatureVertex[i] = std::max(errorsCurvatureVertex[i],
          std::abs(curvaturesVertex[vertexIdx] - exCurv));

        if constexpr (dim == 2)
        {
          errorsCenterVertex[i] = std::max(errorsCenterVertex[i],
            (centersVertex[vertexIdx] -
              exactCenterOfCurvature(center)).infinity_norm());
        }

        exactCurvaturesVertex[vertexIdx] = exCurv;
      }

      if (i == numOfGridFiles - 1)
      { //write output
        Dune::VTKWriter<IGridView> vtkWriter(iGridView);
        vtkWriter.addCellData(curvaturesElement, "curvatureElement");
        vtkWriter.addCellData(exactCurvaturesElement, "exactCurvatureElement");
        vtkWriter.addVertexData(curvaturesVertex, "curvatureVertex");
        vtkWriter.addVertexData(exactCurvaturesVertex, "exactCurvatureVertex");
        vtkWriter.write("curvature");
      }

      errorsCurvatureElement[i] = sqrt(errorsCurvatureElement[i]);

      if constexpr (dim == 2)
      {
        errorsCenterElement[i] = sqrt(errorsCenterElement[i]);
      }

      if constexpr (dim == 3)
      {
        assert(errorsCurvatureElement[i] < 5e-8);
        assert(errorsCurvatureVertex[i] < 5e-8);
      }
    }

    std::cout << "Element Layout:\n\n";
    output(characteristicLengths, gridfiles, numOfIElems,
      errorsCurvatureElement, errorsCenterElement);

    std::cout << "\n\n\nVertex Layout:\n\n";
    output(characteristicLengths, gridfiles, numOfIElems, errorsCurvatureVertex,
      errorsCenterVertex);

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
