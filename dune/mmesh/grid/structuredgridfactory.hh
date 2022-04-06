 // -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_STRUCTURED_GRID_FACTORY_HH
#define DUNE_MMESH_STRUCTURED_GRID_FACTORY_HH

// System includes
#include <memory>
#include <utility>
#include <vector>

// Dune includes
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/common/intersection.hh>
#include <dune/grid/utility/multiindex.hh>

// MMesh includes
#include <dune/mmesh/grid/implicitgridfactory.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------
  template< class GridImp, class IntersectionImp >
  class Intersection;


  // StructuredGridFactory for MMesh
  // -------------------------------
  template< class Grid >
  struct MMeshStructuredGridFactory
  {
    const static int dimworld = Grid::dimension;
    typedef MPIHelper::MPICommunicator MPICommunicatorType;
    typedef typename Grid::template Codim<0>::Entity Element;
    typedef typename Grid::template Codim<dimworld>::Entity Vertex;
    typedef typename Grid::FieldType FieldType;
    typedef FieldVector< FieldType, dimworld > GlobalPosition;
    typedef Dune::MMeshImplicitGridFactory<Grid> GridFactory;

    //! Constructor
    template <long unsigned int dim>
    explicit MMeshStructuredGridFactory (const GlobalPosition& lowerLeft,
                                         const GlobalPosition& upperRight,
                                         const std::array<unsigned int,dim>& elements,
                                         MPICommunicatorType comm = MPIHelper::getCommunicator() )
    {
        grid_ = createStructuredGrid (lowerLeft, upperRight, elements);
    }

    //! return grid pointer
    const std::shared_ptr<Grid> grid() const
    {
      return grid_;
    }

    //! returns if intersection was inserted
    template< class Intersection >
    bool wasInserted ( const Intersection &intersection ) const
    {
      return factory_.wasInserted( intersection );
    }

    //! return boundary id
    template< class Intersection >
    int boundaryId ( const Intersection &intersection ) const
    {
      DUNE_THROW( NotImplemented, "BoundaryId in structured grid factory" );
    }

public:

    /** \brief Create a structured simplicial mmesh grid

        If the grid dimension is less than the world dimension, the coefficients (dim+1,...,dimworld) in
        the vertex coordinates are set to the corresponding values of the lowerLeft input argument.

        \param lowerLeft Lower left corner of the grid
        \param upperRight Upper right corner of the grid
        \param elements Number of elements in each coordinate direction
     */
    template <long unsigned int dim>
    static typename Grid::GridPtrType createStructuredGrid(const GlobalPosition& lowerLeft,
                                                          const GlobalPosition& upperRight,
                                                          const std::array<unsigned int,dim>& elements)
    {
      // The grid factory
      GridFactory factory;

      // Insert uniformly spaced vertices
      std::array<unsigned int,dim> vertices = elements;
      for( size_t i = 0; i < vertices.size(); ++i )
        vertices[i]++;

      FactoryUtilities::MultiIndex<dim> index(vertices);

      // Compute the total number of vertices to be created
      int numVertices = index.cycle();

      // Create vertices
      for (int i=0; i<numVertices; i++, ++index)
      {
        // scale the multiindex to obtain a world position
        GlobalPosition pos(0);
        for (std::size_t j=0; j<dim; j++)
          pos[j] = lowerLeft[j] + index[j] * (upperRight[j]-lowerLeft[j])/(vertices[j]-1);
        for (std::size_t j=dim; j<dimworld; j++)
          pos[j] = lowerLeft[j];

        factory.insertVertex(pos);
      }

      // Create the grid and hand it to the calling method
      return factory.createGrid();

    }

private:
    std::shared_ptr<Grid> grid_;
    GridFactory factory_;
  };

} // end namespace Dune

#endif // #ifndef DUNE_MMESH_STRUCTURED_GRID_FACTORY_HH
