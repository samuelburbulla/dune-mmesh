#include <iostream>
#include <cmath>
#include <string>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>

#include <dune/grid/common/mcmgmapper.hh>

#include <dune/mmesh/mmesh.hh>
#include <dune/mmesh/interface/curvatureoperator.hh>

static constexpr int dim = GRIDDIM;
using Vector = Dune::FieldVector<double, dim>;

#if GRIDDIM == 2
double exactRadiusOfCurvature (Vector pos, double a = 0.3, double b = 0.15)
{
  double a4 = pow(a, 4);
  double b4 = pow(b, 4);
  return pow(a4 * pow(pos[1], 2) + b4 * pow(pos[0], 2), 1.5) / (a4*b4);
}
#elif GRIDDIM == 3
double exactRadiusOfCurvature (Vector pos, double a = 0.35, double b = 0.15,
  double c = 0.25)
{
  const double a2 = a * a;
  const double b2 = b * b;
  const double c2 = c * c;
  const double x2 = pos[0] * pos[0];
  const double y2 = pos[1] * pos[1];
  const double z2 = pos[2] * pos[2];
  const double R2 = x2 + y2 + z2;
  const double H = 1.0 / sqrt(x2 / (a2 * a2) + y2 / (b2 * b2) + z2 / (c2 * c2));
  const double H2 = H * H;

  const double meanCurvature =
    0.5 * H2 * H * (a2 + b2 + c2 - R2) / (a2 * b2 * c2);
  const double gaussianCurvature = H2 * H2 / (a2 * b2 * c2);

  return meanCurvature;
}
#endif

#if GRIDDIM == 2
Vector exactCenterOfCurvature (Vector pos, double a = 0.3, double b = 0.15)
{
  double e2 = pow(a, 2) - pow(b, 2);
  Vector center;
  center[0] = e2 * pow(pos[0], 3) / pow(a, 4);
  center[1] = -e2 * pow(pos[1], 3) / pow(b, 4);
  return center;
}
#endif

int main(int argc, char** argv)
{
  try
  {
    using Grid = Dune::MovingMesh<dim>;
    using GridFactory = Dune::GmshGridFactory<Grid, /*implicit=*/false>;
    using IGridView = typename Grid::InterfaceGrid::LeafGridView;
    using IMapper =
      typename Dune::MultipleCodimMultipleGeomTypeMapper<IGridView>;
    using CurvatureOp = Dune::CurvatureOperator<IGridView, IMapper>;

    Dune::MPIHelper::instance(argc, argv);

    #if GRIDDIM == 2

    static constexpr int numOfGridFiles = 11;
    const std::array<std::string, numOfGridFiles>
      gridfiles({"9e-2", "7e-2", "5e-2", "3e-2", "2e-2", "1e-2", "9e-3",
      "8e-3", "7e-3", "6e-3", "5e-3"});
    const std::array<double, numOfGridFiles> characteristicLengths({9e-2, 7e-2,
      5e-2, 3e-2, 2e-2, 1e-2, 9e-3, 8e-3, 7e-3, 6e-3, 5e-3});

    #elif GRIDDIM == 3

    static constexpr int numOfGridFiles = 5;
    const std::array<std::string, numOfGridFiles>
      gridfiles({"9e-2", "7e-2", "5e-2", "3e-2", "2e-2"});
    const std::array<double, numOfGridFiles>
      characteristicLengths({9e-2, 7e-2, 5e-2, 3e-2, 2e-2});

    #endif

    std::array<double, numOfGridFiles> errorsCurvature;
    std::array<double, numOfGridFiles> errorsCenter;
    std::array<double, numOfGridFiles> numOfIElems;

    for (int i = 0; i < numOfGridFiles; i++)
    {
      errorsCurvature[i] = 0.0;
      errorsCenter[i] = 0.0;

      //get interface grid view
      GridFactory gridFactory("grids/ellipse" + std::to_string(dim) + "d_"
        + gridfiles[i] + ".msh");
      const IGridView& iGridView =
        gridFactory.grid()->interfaceGrid().leafGridView();

      //get number of interface segments
      numOfIElems[i] = iGridView.size(0);

      //define element mapper
      IMapper iElemMapper(iGridView, Dune::mcmgElementLayout());

      //storage curvatures on the interface grid
      std::vector<double> curvatures(iElemMapper.size());
      std::vector<Vector> centers(iElemMapper.size());

      std::vector<double> exactCurvatures(iElemMapper.size());

      //determine curvatures and sphere center points
      CurvatureOp curvatureOp(iGridView, iElemMapper);
      curvatureOp(curvatures, centers);

      //calculate L2 error
      for (const auto& iElem: elements(iGridView))
      {
        const int iElemIdx = iElemMapper.index(iElem);
        const auto iGeo = iElem.geometry();
        const auto& center = iGeo.center();
        const double iVolume = iGeo.volume();
        double exactCurvature =  1.0 / exactRadiusOfCurvature(center);

        errorsCurvature[i] +=
          pow(curvatures[iElemIdx] - exactCurvature, 2)*iVolume;

        #if GRIDDIM == 2
        errorsCenter[i] += (centers[iElemIdx]
          - exactCenterOfCurvature(center)).two_norm2()*iVolume;
        #endif

        exactCurvatures[iElemIdx] = exactCurvature;
      }

      if (i == numOfGridFiles - 1)
      { //write output
        Dune::VTKWriter<IGridView> vtkWriter(iGridView);
        vtkWriter.addCellData(curvatures, "curvature");
        vtkWriter.addCellData(exactCurvatures, "exactCurvature");
        vtkWriter.write("curvature");
      }

      errorsCurvature[i] = sqrt(errorsCurvature[i]);

      #if GRIDDIM == 2
      errorsCenter[i] = sqrt(errorsCenter[i]);
      #endif
    }

    for (int i = 0; i < numOfGridFiles; i++)
    {
      std::cout << "characteristic length: " << gridfiles[i] <<
        "\t#interface elements: " << numOfIElems[i] << " \terror curvature: "
        << errorsCurvature[i]
        #if GRIDDIM == 2
        << "\terror center: " << errorsCenter[i]
        #endif
        << std::endl;
    }

    std::cout << "\nEOC curvature: " << std::endl;
    for (int i = 1; i < numOfGridFiles; i++)
    {
      const double EOC = log(errorsCurvature[i-1] / errorsCurvature[i])
      / log(characteristicLengths[i-1] / characteristicLengths[i]);
      std::cout << "EOC: " << EOC << "\t(interface elements " <<
        numOfIElems[i-1] << " -> " << numOfIElems[i] << ")\n";
    }

    #if GRIDDIM == 2
    std::cout << "\nEOC center of curvature: " << std::endl;
    for (int i = 1; i < numOfGridFiles; i++)
    {
      const double EOC = log(errorsCenter[i-1] / errorsCenter[i])
      / log(characteristicLengths[i-1] / characteristicLengths[i]);
      std::cout << "EOC: " << EOC << "\t(interface elements " <<
        numOfIElems[i-1] << " -> " << numOfIElems[i] << ")\n";
    }
    #endif

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
