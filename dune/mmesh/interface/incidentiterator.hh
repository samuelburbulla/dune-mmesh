// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_INTERFACE_INCIDENTITERATOR_HH
#define DUNE_MMESH_INTERFACE_INCIDENTITERATOR_HH

/** \file
 * \brief The MMeshIncidentIterator class
 */

// Dune includes
#include <dune/grid/common/gridenums.hh>

namespace Dune
{

  /** \brief Iterator over all incident interface vertices
   *  \ingroup MMesh
   */

  //! Forward declaration
  template<class GridImp, int dim>
  class MMeshIncidentInterfaceVerticesIteratorImp;

  //!  The Incident Interface Vertices Iterator alias
  template<class Grid>
  using MMeshIncidentInterfaceVerticesIterator = EntityIterator<Grid::dimension, Grid, MMeshIncidentInterfaceVerticesIteratorImp<Grid, Grid::dimension>>;

  //! 2D
  template<class GridImp>
  class MMeshIncidentInterfaceVerticesIteratorImp<GridImp, 1>
  {
  private:
    //! The type of the requested vertex entity
    typedef typename GridImp::VertexHandle HostGridVertex;

    //! The type of the element circulator
   using Circulator = typename GridImp::HostGridType::Vertex_circulator;
   using ElementContainer = std::vector<HostGridVertex>;

  public:
    enum {codimension = 1};

    typedef typename GridImp::template Codim<1>::Entity Entity;

    explicit MMeshIncidentInterfaceVerticesIteratorImp(const GridImp* igrid, const HostGridVertex& hostEntity)
    : mMesh_(&igrid->getMMesh()),
      i_(0)
    {
      Circulator circulator = mMesh_->getHostGrid().incident_vertices(hostEntity);
      for ( std::size_t i = 0; i < CGAL::circulator_size(circulator); ++i, ++circulator )
        if (circulator->info().isInterface)
          elementContainer_.push_back( circulator );
    }

    /** \brief Constructor which creates the end iterator
     *  \param endDummy      Here only to distinguish it from the other constructor
     */
    explicit MMeshIncidentInterfaceVerticesIteratorImp(const GridImp* igrid, const HostGridVertex& hostEntity, bool endDummy)
    : mMesh_(&igrid->getMMesh()),
      i_(0)
    {
      Circulator circulator = mMesh_->getHostGrid().incident_vertices(hostEntity);
      for ( std::size_t i = 0; i < CGAL::circulator_size(circulator); ++i, ++circulator )
        if (circulator->info().isInterface)
          ++i_;
    }

    //! prefix increment
    void increment() {
      ++i_;
    }

    //! dereferencing
    Entity dereference() const {
      return Entity {{ &mMesh_->interfaceGrid(), elementContainer_[i_] }};
    }

    //! equality
    bool equals(const MMeshIncidentInterfaceVerticesIteratorImp& iter) const {
      return i_ == iter.i_;
    }

  private:
    const typename GridImp::MMeshType* mMesh_;
    ElementContainer elementContainer_;
    std::size_t i_;
  };


  //! 3D
  template<class GridImp>
  class MMeshIncidentInterfaceVerticesIteratorImp<GridImp, 2>
  {
  private:
    //! The type of the requested vertex entity
    typedef typename GridImp::VertexHandle HostGridVertex;

    //! The type of the element output
    using ElementOutput = std::list<HostGridVertex>;

    //! The type of the element container
    using ElementContainer = std::vector<HostGridVertex>;

  public:
    enum {codimension = 2};

    typedef typename GridImp::template Codim<2>::Entity Entity;

    explicit MMeshIncidentInterfaceVerticesIteratorImp(const GridImp* igrid, const HostGridVertex& hostEntity)
    : mMesh_(&igrid->getMMesh()),
      i_(0)
    {
      ElementOutput elements;
      mMesh_->getHostGrid().incident_vertices( hostEntity, std::back_inserter(elements) );

      typename ElementOutput::iterator fit;
      for(fit = elements.begin(); fit != elements.end(); fit++)
        if ((*fit)->info().isInterface)
          elementContainer_.push_back( *fit );
    }

    /** \brief Constructor which creates the end iterator
     *  \param endDummy      Here only to distinguish it from the other constructor
     */
    explicit MMeshIncidentInterfaceVerticesIteratorImp(const GridImp* igrid, const HostGridVertex& hostEntity, bool endDummy)
    : mMesh_(&igrid->getMMesh()),
      i_(0)
    {
      ElementOutput elements;
      mMesh_->getHostGrid().incident_vertices( hostEntity, std::back_inserter(elements) );

      typename ElementOutput::iterator fit;
      for(fit = elements.begin(); fit != elements.end(); fit++)
        if ((*fit)->info().isInterface)
          ++i_;
    }

    //! prefix increment
    void increment() {
      ++i_;
    }

    //! dereferencing
    Entity dereference() const {
      return Entity {{ &mMesh_->interfaceGrid(), elementContainer_[i_] }};
    }

    //! equality
    bool equals(const MMeshIncidentInterfaceVerticesIteratorImp& iter) const {
      return i_ == iter.i_;
    }

  private:
    const typename GridImp::MMeshType* mMesh_;
    ElementContainer elementContainer_;
    std::size_t i_;
  };


  /** \brief Iterator over all incident interface elements
   *  \ingroup MMesh
   */

  //! Forward declaration
  template<class GridImp, int dim>
  class MMeshIncidentInterfaceElementsIteratorImp;

  //!  The Incident Interface Elements Iterator alias
  template<class Grid>
  using MMeshIncidentInterfaceElementsIterator = EntityIterator<0, Grid, MMeshIncidentInterfaceElementsIteratorImp<Grid, Grid::dimension>>;

  //! 2D
  template<class GridImp>
  class MMeshIncidentInterfaceElementsIteratorImp<GridImp, 1>
  {
  private:
    //! The type of the requested element
    typedef typename GridImp::template MMeshInterfaceEntity<0> HostGridElement;

    //! The type of the vertex entity
    typedef typename GridImp::VertexHandle HostGridVertex;

    //! The type of the element circulator
   using Circulator = typename GridImp::HostGridType::Edge_circulator;
   using ElementContainer = std::vector<HostGridElement>;

  public:
    enum {codimension = 0};

    typedef typename GridImp::template Codim<0>::Entity Entity;

    explicit MMeshIncidentInterfaceElementsIteratorImp(const GridImp* igrid, const HostGridVertex& hostEntity)
    : igrid_(igrid),
      i_(0)
    {
      Circulator circulator = mMesh().getHostGrid().incident_edges(hostEntity);
      for ( std::size_t i = 0; i < CGAL::circulator_size(circulator); ++i, ++circulator )
        if ( mMesh().isInterface( mMesh().entity( *circulator ) ) )
          elementContainer_.push_back( *circulator );
    }

    /** \brief Constructor which creates the end iterator
     *  \param endDummy      Here only to distinguish it from the other constructor
     */
    explicit MMeshIncidentInterfaceElementsIteratorImp(const GridImp* igrid, const HostGridVertex& hostEntity, bool endDummy)
    : igrid_(igrid),
      i_(0)
    {
      Circulator circulator = mMesh().getHostGrid().incident_edges(hostEntity);
      for ( std::size_t i = 0; i < CGAL::circulator_size(circulator); ++i, ++circulator )
        if ( mMesh().isInterface( mMesh().entity( *circulator ) ) )
          ++i_;
    }

    //! prefix increment
    void increment() {
      ++i_;
    }

    //! dereferencing
    Entity dereference() const {
      return Entity {{ igrid_, elementContainer_[i_] }};
    }

    //! equality
    bool equals(const MMeshIncidentInterfaceElementsIteratorImp& iter) const {
      return i_ == iter.i_;
    }

  private:
    const typename GridImp::MMeshType& mMesh() const { return igrid_->getMMesh(); }

    const GridImp* igrid_;
    ElementContainer elementContainer_;
    std::size_t i_;
  };


  //! 3D
  template<class GridImp>
  class MMeshIncidentInterfaceElementsIteratorImp<GridImp, 2>
  {
  private:
    //! The type of the requested element
    typedef typename GridImp::template MMeshInterfaceEntity<0> HostGridElement;

    //! The type of the vertex entity
    typedef typename GridImp::VertexHandle HostGridVertex;

    //! The type of the element output
    using ElementOutput = std::list<HostGridElement>;

    //! The type of the element container
    using ElementContainer = std::vector<HostGridElement>;

  public:
    enum {codimension = 0};

    typedef typename GridImp::template Codim<0>::Entity Entity;

    explicit MMeshIncidentInterfaceElementsIteratorImp(const GridImp* igrid, const HostGridVertex& hostEntity)
    : igrid_(igrid),
      i_(0)
    {
      ElementOutput elements;
      mMesh().getHostGrid().finite_incident_facets( hostEntity, std::back_inserter(elements) );

      typename ElementOutput::iterator fit;
      for(fit = elements.begin(); fit != elements.end(); fit++)
        if ( mMesh().isInterface( mMesh().entity( *fit ) ) )
          elementContainer_.push_back( *fit );
    }

    /** \brief Constructor which creates the end iterator
     *  \param endDummy      Here only to distinguish it from the other constructor
     */
    explicit MMeshIncidentInterfaceElementsIteratorImp(const GridImp* igrid, const HostGridVertex& hostEntity, bool endDummy)
    : igrid_(igrid),
      i_(0)
    {
      ElementOutput elements;
      mMesh().getHostGrid().finite_incident_facets( hostEntity, std::back_inserter(elements) );

      typename ElementOutput::iterator fit;
      for(fit = elements.begin(); fit != elements.end(); fit++)
        if ( mMesh().isInterface( mMesh().entity( *fit ) ) )
          ++i_;
    }

    //! prefix increment
    void increment() {
      ++i_;
    }

    //! dereferencing
    Entity dereference() const {
      return Entity {{ igrid_, elementContainer_[i_] }};
    }

    //! equality
    bool equals(const MMeshIncidentInterfaceElementsIteratorImp& iter) const {
      return i_ == iter.i_;
    }

  private:
    const typename GridImp::MMeshType& mMesh() const { return igrid_->getMMesh(); }

    const GridImp* igrid_;
    ElementContainer elementContainer_;
    std::size_t i_;
  };


  //! Forward declaration
  template<class GridImp, int dim>
  class MMeshEdgeIncidentInterfaceElementsIteratorImp;

  //!  The Incident Interface Elements Iterator alias
  template<class Grid>
  using MMeshEdgeIncidentInterfaceElementsIterator = EntityIterator<0, Grid,
    MMeshEdgeIncidentInterfaceElementsIteratorImp<Grid, Grid::dimension>>;

  //! 3D
  template<class GridImp>
  class MMeshEdgeIncidentInterfaceElementsIteratorImp<GridImp, 2>
  {
  private:
    //! The type of the requested element
    typedef typename GridImp::template MMeshInterfaceEntity<0> HostGridElement;

    //! The type of the vertex entity
    typedef typename GridImp::EdgeHandle HostGridEdge;

    //! The type of the facet ciruclator
    using Circulator = typename GridImp::HostGridType::Facet_circulator;

    //! The type of the element container
    using ElementContainer = std::vector<HostGridElement>;

  public:
    enum {codimension = 0};

    typedef typename GridImp::template Codim<0>::Entity Entity;

    explicit MMeshEdgeIncidentInterfaceElementsIteratorImp(const GridImp* igrid, const HostGridEdge& hostEntity)
    : igrid_(igrid),
    i_(0)
    {
      Circulator circulator = mMesh().getHostGrid().incident_facets(hostEntity);
      for ( std::size_t i = 0; i < CGAL::circulator_size(circulator); ++i, ++circulator )
        if ( mMesh().isInterface( mMesh().entity( *circulator ) ) )
          elementContainer_.push_back( *circulator );
    }

    /** \brief Constructor which creates the end iterator
     *  \param endDummy      Here only to distinguish it from the other constructor
     */
    explicit MMeshEdgeIncidentInterfaceElementsIteratorImp(const GridImp* igrid, const HostGridEdge& hostEntity, bool endDummy)
    : igrid_(igrid),
    i_(0)
    {
      Circulator circulator = mMesh().getHostGrid().incident_facets(hostEntity);
      for ( std::size_t i = 0; i < CGAL::circulator_size(circulator); ++i, ++circulator )
        if ( mMesh().isInterface( mMesh().entity( *circulator ) ) )
          ++i_;
    }

    //! prefix increment
    void increment() {
      ++i_;
    }

    //! dereferencing
    Entity dereference() const {
      return Entity {{ igrid_, elementContainer_[i_] }};
    }

    //! equality
    bool equals(const MMeshEdgeIncidentInterfaceElementsIteratorImp& iter) const {
      return i_ == iter.i_;
    }

  private:
    const typename GridImp::MMeshType& mMesh() const { return igrid_->getMMesh(); }

    const GridImp* igrid_;
    ElementContainer elementContainer_;
    std::size_t i_;
  };

}  // namespace Dune

#endif
