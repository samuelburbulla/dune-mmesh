// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_MMESH_GRID_IMPLICITGRIDFACTORY_HH
#define DUNE_MMESH_GRID_IMPLICITGRIDFACTORY_HH

// System includes
#include <algorithm>
#include <array>
#include <limits>
#include <map>
#include <memory>
#include <unordered_map>
#include <unordered_set>

// Dune includes
#include <dune/common/to_unique_ptr.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/common/boundarysegment.hh>

// MMesh includes
#include <dune/mmesh/mmesh.hh>
#include <dune/mmesh/grid/common.hh>
#include <dune/mmesh/grid/pointfieldvector.hh>

namespace Dune
{

  /** \brief specialization of the implicit GridFactory for MMesh
   *
   *  \ingroup GridFactory
   */

  //! The implicit grid factory for MMesh
  template< class Grid >
  class MMeshImplicitGridFactory
    : public GridFactoryInterface< Grid >
  {
    typedef MMeshImplicitGridFactory< Grid > This;

  public:
    //! type of (scalar) coordinates
    typedef typename Grid::ctype ctype;
    //! type of the hostgrid
    typedef typename Grid::HostGridType HostGrid;

    //! dimension of the grid
    static const int dimension = Grid::dimension;
    //! dimension of the world
    static const int dimensionworld = Grid::dimensionworld;

    //! type of vector for world coordinates
    typedef FieldVector< ctype, dimensionworld > WorldVector;
    //! type of matrix from world coordinates to world coordinates
    typedef FieldMatrix< ctype, dimensionworld, dimensionworld > WorldMatrix;

    //! type of a Dune boundary segment
    typedef Dune::BoundarySegment< dimension, dimensionworld > BoundarySegment;
    //! type of the boundary segment id map
    typedef std::unordered_map< std::vector< std::size_t >, std::size_t, HashUIntVector > BoundarySegments;

    //! type of the interface segment set
    typedef std::unordered_set< std::vector< std::size_t >, HashUIntVector > InterfaceSegments;

    template< int codim >
    struct Codim
    {
      typedef typename Grid::template Codim< codim >::Entity Entity;
    };

  private:
    typedef typename HostGrid::Point Point;
    typedef typename HostGrid::Vertex_handle Vertex_handle;
    typedef typename Grid::template HostGridEntity<0> Element_handle;
    typedef typename Grid::template HostGridEntity<1> Face_handle;

  public:
    //! are boundary ids supported by this factory?
    static const bool supportsBoundaryIds = true;
    //! the factory is not able to create periodic meshes
    static const bool supportPeriodicity = false;

    /** default constructor */
    MMeshImplicitGridFactory ()
    {}

    /** \brief insert an element into the macro grid
     *
     *  \param[in]  type      GeometryType of the new element
     *  \param[in]  v         indices of the element vertices (starting with 0)
     *
     *  \note The implicit grid factory ignores the insertion of specific elements.
     *  The CGAL backend constructs the elements by inserting vertices implicitly.
     */
    void insertElement ( const GeometryType &type,
                         const std::vector< unsigned int > &v )
    {
      Element_handle fh;
      if( isElement( v, fh ) )
      {
        fh->info().id = countElements;
        fh->info().idWasSet = true;
      }

      // Increase element count in each case
      countElements++;
    };

    /** \brief Returns if there is a face with the given vertices in the triangulation
     *  \ingroup 2D
     *
     *  \param[in]  v        indices of the element vertices
     */
    template< int d = dimension >
    std::enable_if_t< d == 2, bool >
    isElement( const std::vector< unsigned int >& v, Element_handle& fh ) const
    {
        assert( v.size() == dimension+1 );
        return tr_.is_face( vhs_[v[0]], vhs_[v[1]], vhs_[v[2]], fh );
    }

    /** \brief Returns if there is a cell with the given vertices in the triangulation
     *  \ingroup 3D
     *
     *  \param[in]  v        indices of the element vertices
     */
    template< int d = dimension >
    std::enable_if_t< d == 3, bool >
    isElement( const std::vector< unsigned int >& v, Element_handle& fh ) const
    {
        assert( v.size() == dimension+1 );
        return tr_.is_cell( vhs_[v[0]], vhs_[v[1]], vhs_[v[2]], vhs_[v[3]], fh );
    }

    /** \brief Insert a vertex into the macro grid
     *
     *  \param[in]  pos  position of the vertex (in world coordinates)
     *  \note This method assumes that the vertices are inserted consecutively
     *        with respect to their index.
     */
    void insertVertex ( const WorldVector &pos )
    {
      // Insert vertex
      Vertex_handle vh = tr_.insert( makePoint( pos ) );

      // Store insertion index
      vh->info().id = countVertices;
      vh->info().idWasSet = true;

      // Increase vertex counter
      countVertices++;

      // Store vertex handle for later use
      vhs_.push_back( vh );
    }

    /** \brief insert an interface into the macro grid
     *
     *  \param[in]  v         indices of the interface vertices (starting with 0)
     */
    void insertInterface ( const std::vector< unsigned int > &v )
    {
      assert( v.size() == dimension );

      // mark all vertices as interface
      for( std::size_t i = 0; i < v.size(); ++i )
        vhs_[v[i]]->info().isInterface = true;

      std::vector< std::size_t > sorted_vertices( v.begin(), v.end() );
      std::sort(sorted_vertices.begin(), sorted_vertices.end());
      interfaceSegments_.insert( sorted_vertices );
    }

    /** \brief return index of inserted vertex within the macro grid
     *
     *  \param[in]  pos  position of the vertex (in world coordinates)
     */
    unsigned int insertionIndex ( const typename Codim<0>::Entity &entity ) const
    {
      return entity.impl().hostEntity()->info().id;
    }

    /** \brief return insertion index of entity
     *
     *  \param[in]  entity  Entity
     */
    unsigned int insertionIndex ( const typename Codim< dimension >::Entity &entity ) const
    {
      return entity.impl().hostEntity()->info().id;
    }

    /** \brief insert boundary segment
     *
     *  \param[in]  vertices  Vertices
     */
    void insertBoundarySegment(const std::vector<unsigned int>& vertices)
    {
      std::vector< std::size_t > sorted_vertices;
      for ( const auto& v : vertices )
        sorted_vertices.push_back( vhs_[v]->info().id );  // vertices give the number of the vertex
      std::sort(sorted_vertices.begin(), sorted_vertices.end());

      if( boundarySegments_.find( sorted_vertices ) != boundarySegments_.end() )
        DUNE_THROW( GridError, "A boundary segment was inserted twice." );

      boundarySegments_.insert( std::make_pair( sorted_vertices, countBoundarySegments++ ) );
    }

    void insertBoundarySegment(const std::vector<unsigned int>& vertices,
                               const std::shared_ptr<Dune::BoundarySegment<dimension,dimension>>& boundarySegment)
    {
      DUNE_THROW( NotImplemented, "insertBoundarySegments with Dune::BoundarySegment" );
    }

    /** \brief finalize grid creation and hand over the grid
     *
     *  This version of createGrid is original to the MMesh grid factroy,
     *  allowing to specity a grid name.
     *
     *  \returns a pointer to the newly created grid
     *
     *  \note The caller takes responsibility of freeing the memory allocated
     *        for the grid.
     *  \note MMesh's grid factory provides a static method for freeing the
     *        grid (destroyGrid).
     */
    typename Grid::GridPtrType createGrid ()
    {
      /* Make interface segments gabriel
       *  \note This removes all conflicting vertices, i.e. all vertices that lie inside a circumscribing sphere of an interface segment.
       *  This throws an exception if a vertex in conflict is also an interface point.
       *  Remark: Do not delete Vertex_handles from vhs_ (but from tr_) to keep the vertex indices consistent.
       */
      int count = 0;
      for ( const auto& v : interfaceSegments_ )
      {
        const auto sphere = makeSphere_( v );

        // Walk over vertices
        for ( auto vh = tr_.finite_vertices_begin(); vh != tr_.finite_vertices_end(); vh++ )
          if ( sphere.has_on_bounded_side( vh->point() ) )
          {
            // exlude the vertices of v
            if ( std::any_of( v.begin(), v.end(), [this, vh](int i){ return vhs_[i]->info().id == vh->info().id; } ) )
              continue;

            if ( vh->info().isInterface )
              DUNE_THROW( GridError, "The interface is in conflict with the interface point " << vh->point() );

            tr_.remove( vh );
            count++;
          }
      }

      // print some information about removal
      if( count > 0 )
        std::cout << "Removed " << count << " " << ( (count == 1) ? "vertex" : "vertices") << " to make interface gabriel." << std::endl;

      // Return pointer to grid
      #if DUNE_VERSION_NEWER(DUNE_GRID, 2, 7)
        return makeToUnique<Grid>( std::move(tr_), std::move(boundarySegments_), std::move(interfaceSegments_) );
      #else
        return new Grid ( std::move(tr_), std::move(boundarySegments_), std::move(interfaceSegments_) );
      #endif
    }

    /** \brief destroy a grid previously obtained from this factory
     *
     *  \param[in]  grid  pointer to the grid to destroy
     */
    static void destroyGrid ( Grid *grid )
    {
      delete grid;
    }

  private:
    template<int d = dimension>
    inline std::enable_if_t< d == 2, std::vector< Face_handle > >
    getFaces_( const std::vector< std::size_t >& v )
    {
      Element_handle f; int i;
      tr_.is_edge( vhs_[v[0]], vhs_[v[1]], f, i );

      std::vector< Face_handle > faces;

      if ( !tr_.is_infinite( f ) )
        faces.push_back( { f, i } );

      if ( !tr_.is_infinite( f->neighbor( i ) ) )
        faces.push_back( { f->neighbor( i ), tr_.mirror_index( f, i ) } );

      return faces;
    }

    template<int d = dimension>
    inline std::enable_if_t< d == 3, std::vector< Face_handle > >
    getFaces_( const std::vector< std::size_t >& v )
    {
      Element_handle f;
      int i,j,k,l;
      tr_.is_facet( vhs_[v[0]], vhs_[v[1]], vhs_[v[2]], f, i, j, k );
      for( l = 0; l == i || l == j || l == k ; ++l );

      std::vector< Face_handle > faces;

      if ( !tr_.is_infinite( f ) )
        faces.push_back( { f, l } );

      if ( !tr_.is_infinite( f->neighbor( l ) ) )
        faces.push_back( { f->neighbor( l ), tr_.mirror_index( f, l ) } );

      return faces;
    }

    //! Return the circumscribing circle of an interface line segment
    template<int d = dimension>
    inline std::enable_if_t< d == 2, CGAL::Circle_2<typename HostGrid::Point::R> >
    makeSphere_( const std::vector< std::size_t >& v )
    {
      return CGAL::Circle_2<typename HostGrid::Point::R>(
        vhs_[v[0]]->point(),
        vhs_[v[1]]->point()
      );
    }

    //! Return the circumscribing sphere of an interface triangle segment
    template<int d = dimension>
    inline std::enable_if_t< d == 3, CGAL::Sphere_3<typename HostGrid::Point::R> >
    makeSphere_( const std::vector< std::size_t >& v )
    {
      return CGAL::Sphere_3<typename HostGrid::Point::R>(
        vhs_[v[0]]->point(),
        vhs_[v[1]]->point(),
        vhs_[v[2]]->point()
      );
    }

    //! Private members
    HostGrid tr_;
    std::size_t countVertices = 0, countElements = 0, countBoundarySegments = 0;
    std::vector< Vertex_handle > vhs_;
    BoundarySegments boundarySegments_;
    InterfaceSegments interfaceSegments_;
  };

} // end namespace Dune

#endif
