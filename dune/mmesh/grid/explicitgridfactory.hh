// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_MMESH_GRID_EXPLICITGRIDFACTORY_HH
#define DUNE_MMESH_GRID_EXPLICITGRIDFACTORY_HH

// System includes
#include <algorithm>
#include <array>
#include <limits>
#include <map>
#include <memory>

// Dune includes
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/common/boundarysegment.hh>

// MMesh includes
#include <dune/mmesh/grid/mmesh.hh>
#include <dune/mmesh/grid/common.hh>

namespace Dune
{

  /** \brief specialization of the explicit GridFactory for MMesh
   *
   *  \ingroup GridFactory
   */

  //! The explicit grid factory for MMesh
  template< class Grid >
  class MMeshExplicitGridFactory
    : public GridFactoryInterface< Grid >
  {
    typedef MMeshExplicitGridFactory< Grid > This;
    typedef GridFactoryInterface< Grid > Base;

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
    //! type of an id
    typedef typename Grid::IdType IdType;

    //! type of the boundary segment id map
    typedef std::unordered_map< IdType, std::size_t > BoundarySegments;
    typedef std::unordered_map< std::size_t, std::size_t > BoundaryIds;

    //! type of the interface segment set
    typedef std::unordered_map< IdType, std::size_t > InterfaceSegments;

    template< int codim >
    struct Codim
    {
      typedef typename Grid::template Codim< codim >::Entity Entity;
    };

  private:
    typedef typename HostGrid::Point Point;
    typedef typename HostGrid::Vertex_handle Vertex_handle;
    typedef typename Grid::template HostGridEntity<0> Element_handle;

    using Tuple = std::tuple<std::vector< unsigned int >, std::size_t, std::size_t>;

    using Base::insertElement;
  public:
    //! are boundary ids supported by this factory?
    static const bool supportsBoundaryIds = true;
    //! the factory is not able to create periodic meshes
    static const bool supportPeriodicity = false;

    /** default constructor */
    MMeshExplicitGridFactory ()
    {
      tr_.tds().clear();
    }

    /** \brief insert an element into the macro grid
     *
     *  \param[in]  type      GeometryType of the new element
     *  \param[in]  v         indices of the element vertices (starting with 0)
     */
    void insertElement ( const GeometryType &type,
                         const std::vector< unsigned int > &v )
    {
      insertElement( type, v, 0 );
    }

    /** \brief insert an element into the macro grid with a given domain marker
     *
     *  \param[in]  type          GeometryType of the new element
     *  \param[in]  v             indices of the element vertices (starting with 0)
     *  \param[in]  domainMarker  domain marker of element
     */
    void insertElement ( const GeometryType &type,
                         const std::vector< unsigned int > &v,
                         const size_t domainMarker )
    {
      assert( type == GeometryTypes::simplex(dimension) );
      assert( v.size() == dimension+1 );
      auto w = v;

      // Create element
      storeElement( w, countElements, domainMarker );

      // Increase element count
      countElements++;

      // Store vertices to check occurence of element later
      elementVerticesList.push_back( w );
    };

  private:
    //! Store element to sort elements before insertion
    void storeElement( std::vector< unsigned int >& v, const size_t insertionIndex, const size_t domainMarker )
    {
      elements_.push_back( std::make_tuple(v, insertionIndex, domainMarker) );
    }

    /** \brief Creates an element (face) in the underlying triangulation data structure
     *  \ingroup 2D
     *
     *  \param[in]  v     indices of the element vertices
     */
    template< int d = dimension >
    std::enable_if_t< d == 2, void >
    createElement( std::vector< unsigned int >& v, const size_t insertionIndex, const size_t domainMarker )
    {
      auto& p0 = vhs_[v[0]]->point();
      auto& p1 = vhs_[v[1]]->point();
      auto& p2 = vhs_[v[2]]->point();

      // Check orientation of vertices
      auto orientation = (p1.y() - p0.y()) * (p2.x() - p1.x()) - (p1.x() - p0.x()) * (p2.y() - p1.y());
       // if clockwise, swap two vertices
      if( orientation > 0.0 )
        std::swap(v[0], v[1]);

      auto&& v0 = vhs_[v[0]];
      auto&& v1 = vhs_[v[1]];
      auto&& v2 = vhs_[v[2]];

      // Create face with vertices v0, v1, v2
      Element_handle face = tr_.tds().create_face(v0, v1, v2);

      // Set insertion index
      face->info().insertionIndex = insertionIndex;
      face->info().domainMarker = domainMarker;

      // Set this face in vertices v0, v1, v2
      v0->set_face(face);
      v1->set_face(face);
      v2->set_face(face);

      // Add facets to neighbor map to obtain connectivity
      addFacetToMap( { v[0], v[1] }, face, 2 );
      addFacetToMap( { v[0], v[2] }, face, 1 );
      addFacetToMap( { v[1], v[2] }, face, 0 );
    }

    /** \brief Creates an element (cell) in the underlying triangulation data structure
     *  \ingroup 3D
     *
     *  \param[in]  v     indices of the element vertices
     */
    template< int d = dimension >
    std::enable_if_t< d == 3, void >
    createElement( const std::vector< unsigned int >& v, const size_t insertionIndex, const size_t domainMarker )
    {
      auto&& v0 = vhs_[v[0]];
      auto&& v1 = vhs_[v[1]];
      auto&& v2 = vhs_[v[2]];
      auto&& v3 = vhs_[v[3]];

      // Create cell with vertices v0, v1, v2, v3
      Element_handle cell = tr_.tds().create_cell(v0, v1, v2, v3);

      // Set insertion index
      cell->info().insertionIndex = insertionIndex;
      cell->info().domainMarker = domainMarker;

      // Set this cell in vertices v0, v1, v2, v3
      v0->set_cell(cell);
      v1->set_cell(cell);
      v2->set_cell(cell);
      v3->set_cell(cell);

      // Add facets to neighbor map to obtain connectivity
      addFacetToMap( { v[0], v[1], v[2] }, cell, 3 );
      addFacetToMap( { v[0], v[1], v[3] }, cell, 2 );
      addFacetToMap( { v[0], v[2], v[3] }, cell, 1 );
      addFacetToMap( { v[1], v[2], v[3] }, cell, 0 );
    }

    /** \brief Adds the facet index tuples as set into a map to find neighbor relations
     *
     *  \param[in]  v        indices of the facet vertices
     *  \param[in]  element  the element handle where facet belongs to
     *  \param[in]  fi       the index of the facet (= index of neighbor) in element
     */
    void addFacetToMap( const std::vector< unsigned int >& v, const Element_handle& element, const int fi )
    {
      assert( v.size() == dimension );

      // Make a set of the vertex indices
      std::set< unsigned int > facetIndices ( v.begin(), v.end() );

      // Try to insert this set into the neighborMap
      auto entry = neighborMap.insert( { facetIndices, std::pair<Element_handle, int>( element, fi ) } );

      // If facetIndices was already in neighborMap, connect the neighbors and remove the entry
      if( !entry.second )
      {
        auto facet = entry.first->second;
        Element_handle neighbor = facet.first;
        int ni = facet.second;

        element->set_neighbor(fi, neighbor);
        neighbor->set_neighbor(ni, element);

        neighborMap.erase( facetIndices );
      }
    }

  public:
    /** \brief Returns if there is a face with the given vertices in the triangulation
     *  \ingroup 2D
     *
     *  \param[in]  v        indices of the element vertices
     */
    template< int d = dimension >
    std::enable_if_t< d == 2, bool >
    isElement( const std::vector< unsigned int >& v ) const
    {
        assert( v.size() == dimension+1 );
        return tr_.is_face( vhs_[v[0]], vhs_[v[1]], vhs_[v[2]] );
    }

    /** \brief Returns if there is a cell with the given vertices in the triangulation
     *  \ingroup 3D
     *
     *  \param[in]  v        indices of the element vertices
     */
    template< int d = dimension >
    std::enable_if_t< d == 3, bool >
    isElement( const std::vector< unsigned int >& v ) const
    {
        assert( v.size() == dimension+1 );
        Element_handle cell;
        return tr_.is_cell( vhs_[v[0]], vhs_[v[1]], vhs_[v[2]], vhs_[v[3]], cell );
    }

    /** \brief insert a boundary segment into the macro grid
     *
     *  Only influences the ordering of the boundary segments
     *  \param[in]  vertices         vertex indices of boundary face
     */

    virtual void insertBoundarySegment ( const std::vector< unsigned int >& vertices )
    {
      assert( vertices.size() == dimension );

      std::vector< std::size_t > sorted_vertices;
      for ( const auto& v : vertices )
        sorted_vertices.push_back( vhs_[v]->info().id );
      std::sort(sorted_vertices.begin(), sorted_vertices.end());

      if( boundarySegments_.find( sorted_vertices ) != boundarySegments_.end() )
        DUNE_THROW( GridError, "A boundary segment was inserted twice." );

      boundarySegments_.insert( std::make_pair( sorted_vertices, countBoundarySegments++ ) );
    }

    void insertBoundarySegment ( const std::vector< unsigned int >& vertices,
                                 const std::shared_ptr< BoundarySegment >& boundarySegment )
    {
      DUNE_THROW( NotImplemented, "insertBoundarySegments with Dune::BoundarySegment" );
    }

    void insertInterfaceBoundarySegment ( const std::vector< unsigned int >& vertices )
    {
      assert( vertices.size() == dimension-1 );

      std::vector< std::size_t > sorted_vertices;
      for ( const auto& v : vertices )
        sorted_vertices.push_back( vhs_[v]->info().id );
      std::sort(sorted_vertices.begin(), sorted_vertices.end());

      if( boundarySegments_.find( sorted_vertices ) != boundarySegments_.end() )
          DUNE_THROW( GridError, "A boundary segment was inserted twice." );

      interfaceBoundarySegments_.insert( std::make_pair( sorted_vertices, countInterfaceBoundarySegments++ ) );
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
      auto vh = tr_.tds().create_vertex();
      vh->set_point( makePoint( pos ) );

      // Store insertion index
      vh->info().id = countVertices;
      vh->info().idWasSet = true;

      // Increase vertex counter
      countVertices++;

      // Store the inserted vertex handle for later use
      vhs_.push_back( vh );
    }

    /** \brief insert an interface into the macro grid
     *
     *  \param[in]  vertices  indices of the interface vertices (starting with 0)
     *  \param[in]  marker    marker value of the interface segment (default 1)
     */
    void insertInterface ( const std::vector< unsigned int > &vertices, const std::size_t marker = 1 )
    {
      assert( vertices.size() == dimension );

      std::vector< std::size_t > sorted_vertices;
      for( const auto& v : vertices )
        sorted_vertices.push_back( vhs_[v]->info().id );
      std::sort(sorted_vertices.begin(), sorted_vertices.end());
      interfaceSegments_.insert( std::make_pair( sorted_vertices, marker ) );
    }

    /** \brief return insertion index of entity
     *
     *  \param[in]  entity  Entity of codim 0
     */
    unsigned int insertionIndex ( const typename Codim<0>::Entity &entity ) const
    {
      return entity.impl().hostEntity()->info().insertionIndex;
    }

    /** \brief return insertion index of vertex entity
     *
     *  \param[in]  entity  Entity of codim dimension
     */
    unsigned int insertionIndex ( const typename Codim< dimension >::Entity &entity ) const
    {
      return entity.impl().hostEntity()->info().id;
    }

    /** \brief return insertion index of boundary intersection
     *
     *  \param[in]  intersection  Leaf intersection
     */
    unsigned int insertionIndex ( const typename Grid::LeafIntersection &intersection ) const
    {
      return intersection.impl().boundaryId();
    }

    //! returns the boundary segment to index map
    const BoundarySegments& boundarySegments() const
    {
      return boundarySegments_;
    }

    //! returns the boundary segment index to boundary id map
    const BoundaryIds& boundaryIds() const
    {
      return boundaryIds_;
    }

    //! add a boundary id
    void addBoundaryId( std::size_t boundarySegmentIndex, std::size_t boundaryId )
    {
      boundaryIds_.insert( std::make_pair( boundarySegmentIndex, boundaryId ) );
    }

    /** \brief finalize grid creation and hand over the grid
     *
     *  This version of createGrid is original to the MMesh grid factroy,
     *  allowing to specity a grid name.
     *
     *  \returns a pointer to the newly created grid
     */

    std::unique_ptr<Grid> createGrid ()
    {
      // Sort elements by x-coordinate of center for partitioning
      std::sort(elements_.begin(), elements_.end(), [this](const auto& a, const auto& b){
        const auto& va = std::get<0>(a);
        const auto& vb = std::get<0>(b);

        double xa = 0.0;
        double xb = 0.0;
        for (int i = 0; i < dimensionworld+1; ++i)
        {
          xa += vhs_[va[i]]->point().x();
          xb += vhs_[vb[i]]->point().x();
        }
        return xa < xb;
      });

      for (auto& t : elements_)
        createElement(std::get<0>(t), std::get<1>(t), std::get<2>(t));

      // Create the infinite cells (neighbors of boundary cells)
      createInfiniteVertex();

      // Create the infinite cells (neighbors of boundary cells)
      createInfiniteCells();

      // Remove interfaceSegments_ from boundarySegments_
      for( const auto& interfaceSeg : interfaceSegments_ )
        boundarySegments_.erase( interfaceSeg.first );

      // Mark interface vertices as isInterface
      for( const auto& interfaceSeg : interfaceSegments_ )
        for( const auto& v : interfaceSeg.first.vt() )
          vhs_[v]->info().isInterface = true;

      // Check if all inserted elements really exist in the triangulation
      checkOccurenceOfAllElements();

      // Return pointer to grid
      return std::make_unique<Grid> (
        std::move(tr_),
        std::move(boundarySegments_),
        std::move(interfaceBoundarySegments_),
        std::move(boundaryIds_),
        std::move(interfaceSegments_)
      );
    }

    //! return the vertex handles
    const std::vector< Vertex_handle >& vertexHandles () const
    {
      return vhs_;
    }

  private:
    /** \brief Create the infinite vertex
     *
     */
    void createInfiniteVertex()
    {
      if( countElements > 0 )
        tr_.tds().set_dimension(dimension);
      else
        tr_.tds().set_dimension(0);

      Vertex_handle infinite = tr_.tds().create_vertex();
      infinite->info().id = std::size_t(-1);
      tr_.set_infinite_vertex(infinite);
    }

    /** \brief Create the infinite cells (neighbors of boundary cells)
     *  \ingroup 2D
     *
     */
    template< int d = dimension >
    std::enable_if_t< d == 2, void >
    createInfiniteCells()
    {
      // Iterate over all unique facets
      for ( const auto& entry : neighborMap )
      {
        auto facet = entry.second;
        Element_handle face = facet.first;
        int fi = facet.second;

        // Remove real boundary segments from interfaceSegments_
        std::vector< std::size_t > vertices;
        vertices.push_back( face->vertex( (fi+2)%3 )->info().id );
        vertices.push_back( face->vertex( (fi+1)%3 )->info().id );
        std::sort(vertices.begin(), vertices.end());

        auto it = boundarySegments_.find( vertices );
        if( it != boundarySegments_.end() )
          interfaceSegments_.erase( vertices );


        // Create infinite face with correct orientation
        Element_handle iface = tr_.tds().create_face(
          face->vertex( (fi+2)%3 ),
          face->vertex( (fi+1)%3 ),
          tr_.infinite_vertex()
        );

        tr_.infinite_vertex()->set_face( iface );

        face->set_neighbor(fi, iface);
        iface->set_neighbor(2, face);


        // Map infinite neighbors
        for( int i = 0; i < 2; ++i )
        {
          std::set< std::size_t > vertexIndex ( { iface->vertex(i)->info().id } );
          auto entry = infiniteNeighborMap.insert( { vertexIndex, std::pair<Element_handle, int>( iface, (i+1)%2 ) } );

           // If vertexIndex was already in map, connect neighbors and remove entry
          if( !entry.second )
          {
            auto pair = entry.first->second;
            Element_handle neighbor = pair.first;
            int ni = pair.second;

            iface->set_neighbor((i+1)%2, neighbor);
            neighbor->set_neighbor(ni, iface);

            infiniteNeighborMap.erase( vertexIndex );
          }
        }
      }

      // Assert that all boundary facets are mapped
      assert( infiniteNeighborMap.size() == 0 );
    }

    /** \brief Create the infinite cells (neighbors of boundary cells)
     *  \ingroup 3D
     *
     */
    template< int d = dimension >
    std::enable_if_t< d == 3, void >
    createInfiniteCells()
    {
      // Iterate over all unique facets
      for ( const auto& entry : neighborMap )
      {
        auto facet = entry.second;
        Element_handle cell = facet.first;
        int fi = facet.second;

        // Remove real boundary segments from interfaceSegments_
        std::vector< std::size_t > vertices;
        vertices.push_back( cell->vertex( (fi%2==1) ? (fi+2)&3 : (fi+1)&3 )->info().id );
        vertices.push_back( cell->vertex( (fi%2==1) ? (fi+1)&3 : (fi+2)&3 )->info().id );
        vertices.push_back( cell->vertex( (fi+3)&3 )->info().id );
        std::sort(vertices.begin(), vertices.end());

        auto it = boundarySegments_.find( vertices );
        if( it != boundarySegments_.end() )
          interfaceSegments_.erase( vertices );


        // Create infinite cell with correct orientation
        Element_handle icell = tr_.tds().create_cell(
          cell->vertex( (fi%2==1) ? (fi+2)&3 : (fi+1)&3 ),
          cell->vertex( (fi%2==1) ? (fi+1)&3 : (fi+2)&3 ),
          cell->vertex((fi+3)&3),
          tr_.infinite_vertex()
        );

        tr_.infinite_vertex()->set_cell( icell );

        cell->set_neighbor(fi, icell);
        icell->set_neighbor(3, cell);


        // Map infinite neighbors
        for( int i = 0; i < 3; ++i )
        {
          std::set< std::size_t > edgeIndices ( { icell->vertex(i)->info().id, icell->vertex((i+1)%3)->info().id } );
          auto entry = infiniteNeighborMap.insert( { edgeIndices, std::pair<Element_handle, int>( icell, (i+2)%3 ) } );

           // If edgeIndices was already in map, connect neighbors and remove entry
          if( !entry.second )
          {
            auto pair = entry.first->second;
            Element_handle neighbor = pair.first;
            int ni = pair.second;

            icell->set_neighbor((i+2)%3, neighbor);
            neighbor->set_neighbor(ni, icell);

            infiniteNeighborMap.erase( edgeIndices );
          }
        }
      }

      // Assert that all boundary facets are mapped
      assert( infiniteNeighborMap.size() == 0 );
    }

    /** \brief Check if all inserted elements really exist in the triangulation
     */
    void checkOccurenceOfAllElements() const
    {
      for ( const auto& v : elementVerticesList )
        if( !isElement( v ) )
          DUNE_THROW( InvalidStateException, "Inserted element was not found in CGAL triangulation." );
    }

    //! Private members
    HostGrid tr_;
    std::vector< Vertex_handle > vhs_;
    std::vector<Tuple> elements_;
    BoundarySegments boundarySegments_, interfaceBoundarySegments_;
    BoundaryIds boundaryIds_;
    InterfaceSegments interfaceSegments_;
    std::vector< std::vector< unsigned int > > elementVerticesList;
    std::map< std::set< unsigned int >, std::pair<Element_handle, int> > neighborMap;
    std::map< std::set< std::size_t >, std::pair<Element_handle, int> > infiniteNeighborMap;
    std::size_t countElements = 0, countVertices = 0;
    std::size_t countBoundarySegments = 0, countInterfaceBoundarySegments = 0;
  };

} // end namespace Dune

#endif
