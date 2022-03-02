 // -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_GRID_GMSHGRIDFACTORY_HH
#define DUNE_MMESH_GRID_GMSHGRIDFACTORY_HH

// System includes
#include <memory>
#include <utility>
#include <vector>

// Dune includes
#include <dune/grid/common/intersection.hh>

// MMesh includes
#include <dune/mmesh/grid/explicitgridfactory.hh>
#include <dune/mmesh/grid/implicitgridfactory.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------
  template< class GridImp, class IntersectionImp >
  class Intersection;

  // GmshGridFactory for MMesh
  // ------------------------
  template< class Grid, bool useImplicitGridFactory = false >
  struct GmshGridFactory
  {
    const static int dimension = Grid::dimension;
    typedef MPIHelper::MPICommunicator MPICommunicatorType;
    typedef typename Grid::template Codim<0>::Entity Element;
    typedef typename Grid::template Codim<dimension>::Entity Vertex;
    typedef typename Grid::FieldType FieldType;
    typedef FieldVector< FieldType, dimension > GlobalPosition;

    // Choose explicit or implicit grid factory
    typedef std::conditional_t< useImplicitGridFactory,
      MMeshImplicitGridFactory<Grid>,
      MMeshExplicitGridFactory<Grid>
    > GridFactory;

    //! Constructor using istream
    explicit GmshGridFactory ( std::istream &input )
    {
      input.clear();
      input.seekg( 0 );
      if( !input )
        DUNE_THROW(DGFException, "Error resetting input stream." );
      generate( input );
    }

    //! Constructor using filename
    explicit GmshGridFactory ( const std::string &filename )
    {
      factory_ = new GridFactory();
      std::ifstream input( filename.c_str() );
      if( !input )
        DUNE_THROW( DGFException, "Macrofile " << filename << " not found." );
      if( !generate( input ) )
        DUNE_THROW( DGFException, "Could not generate MMesh from macrofile " << filename );
      input.close();
    }

    explicit GmshGridFactory ( Dune::GridFactory< Grid >& gridFactory, const std::string &filename )
    {
      factory_ = &gridFactory;
      std::ifstream input( filename.c_str() );
      if( !input )
        DUNE_THROW( DGFException, "Macrofile " << filename << " not found." );
      if( !generate( input ) )
        DUNE_THROW( DGFException, "Could not generate MMesh from macrofile " << filename );
      input.close();
    }

    //! return grid pointer
    std::shared_ptr<Grid> grid() const
    {
      if ( !grid_ )
        grid_ = factory_->createGrid();

      return grid_;
    }

    //! returns if intersection was inserted
    template< class Intersection >
    bool wasInserted ( const Intersection &intersection ) const
    {
      return factory_->wasInserted( intersection );
    }

  private:
    bool generate( std::istream &gridFile )
    {
      // read file until we get to the list of nodes
      std::string line;
      std::getline(gridFile, line);
      while (line.find("$Nodes") == std::string::npos)
      {
        if (!std::getline(gridFile, line))
          DUNE_THROW(IOError, "Gmsh file not valid!");
      }

      // read all vertices
      std::getline(gridFile, line);
      const auto numVertices = convertString<std::size_t>(line);

      std::getline(gridFile, line);
      std::size_t vertexCount = 0;
      while (line.find("$EndNodes")  == std::string::npos)
      {
          // drop first entry in line (vertex index) and read coordinates
          std::istringstream stream(line);
          std::string buf; stream >> buf;
          GlobalPosition v;
          for (auto& coord : v)
          {
              stream >> coord;
              if (stream.fail()) DUNE_THROW(Dune::IOError, "Could not read vertex coordinate");
          }

          // insert vertex and move to next line
          factory_->insertVertex( v );
          if (!std::getline(gridFile, line))
            DUNE_THROW(IOError, "Gmsh file not valid!");
          vertexCount++;
      }

      // we should always find as many vertices as the mesh file states
      if (vertexCount != numVertices)
          DUNE_THROW(Dune::InvalidStateException, "Couldn't find as many vertices as stated in the .msh file");

      // read file until we get to the list of elements
      while(line.find("$Elements") == std::string::npos)
        if (!std::getline(gridFile, line))
          DUNE_THROW(IOError, "Gmsh file not valid!");

      // read elements
      std::getline(gridFile, line);
      const auto numElements = convertString<std::size_t>(line);

      std::size_t elemCount = 0;
      std::getline(gridFile, line);
      while (line.find("$EndElements") == std::string::npos)
      {
          // pass all indices into vector
          std::istringstream stream(line);
          std::string buf;
          std::vector<std::size_t> lineData;
          while (stream >> buf) lineData.push_back(convertString<std::size_t>(buf));
          assert(lineData.size() >= 4 && "Grid format erroneous or unsupported");

          // obtain geometry type
          const auto gt = obtainGeometryType( lineData[1] );

          std::vector< unsigned int > cornersIndices;
          auto it = lineData.begin()+2+lineData[2]+1;
          for (; it != lineData.end(); ++it)
              cornersIndices.push_back(*it-1); // gmsh indices start from 1

          // insert element
          if ( gt.dim() == dimension )
              factory_->insertElement( gt, cornersIndices, lineData[2+lineData[2]] );

          // insert interface/boundary segments
          if ( gt.dim() == dimension-1 )
          {
              factory_->insertInterface( cornersIndices, lineData[2+lineData[2]] );
              factory_->insertBoundarySegment( cornersIndices );
          }

          // insert interface's boundary segments
          if ( gt.dim() == dimension-2 )
              factory_->insertInterfaceBoundarySegment( cornersIndices );

          // get next line
          if (!std::getline(gridFile, line))
            DUNE_THROW(IOError, "Gmsh file not valid!");
          elemCount++;
      }

      // make sure we read all elements
      if (elemCount != numElements)
          DUNE_THROW(Dune::InvalidStateException, "Didn't read as many elements as stated in the .msh file");

      return true;
    }

  private:
    //! converts a value contained in a string
    template<class T>
    T convertString(const std::string& string) const
    {
        T value;
        std::istringstream stream(string);
        stream >> value;
        if (stream.fail())
            DUNE_THROW(Dune::InvalidStateException, "Couldn't convert string: " << string << "to type: " << typeid(T).name());
        return value;
    }

    //! obtain Dune::GeometryType from a given gmsh element type
    Dune::GeometryType obtainGeometryType(std::size_t gmshElemType) const
    {
        switch (gmshElemType)
        {
            case 15: return Dune::GeometryTypes::vertex;        // points
            case 1:  return Dune::GeometryTypes::line;          // lines
            case 2:  return Dune::GeometryTypes::triangle;      // triangle
            case 3:  return Dune::GeometryTypes::quadrilateral; // quadrilateral
            case 4:  return Dune::GeometryTypes::tetrahedron;   // tetrahedron
            case 5:  return Dune::GeometryTypes::hexahedron;    // hexahedron
            default:
                DUNE_THROW(Dune::NotImplemented, "FacetCoupling gmsh reader for gmsh element type " << gmshElemType);
        }
    }

    mutable std::shared_ptr<Grid> grid_;
    GridFactory* factory_;
  };

} // end namespace Dune

#endif
