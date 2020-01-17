// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_GRID_MULTIID_HH
#define DUNE_MMESH_GRID_MULTIID_HH

/** \file
 * \brief The multi id class
 */

#include <dune/mmesh/grid/common.hh>

namespace Dune
{
  namespace Impl
  {
    class MultiId
    {
    public:
      using ThisType = MultiId;
      using T = std::size_t;
      using VT = std::vector<T>;

      MultiId() {}

      MultiId( const VT& vt )
       : vt_ ( vt )
      {}

      MultiId( std::initializer_list<T> l )
       : vt_ ( l )
      {}

      MultiId( T t )
       : vt_ ( { t } )
      {}

      operator std::size_t() const
      {
        assert( vt_.size() == 1 );
        return vt_[0];
      }

      ThisType& operator= ( const ThisType& b )
      {
        if (this != &b)
          vt_ = b.vt_;
        return *this;
      }

      bool operator< ( const ThisType& b ) const
      {
        if( size() != b.size() )
          return size() < b.size();
        return vt_ < b.vt_;
      }

      bool operator== ( const ThisType& b ) const
      {
        return vt_ == b.vt_;
      }

      bool operator!= ( const ThisType& b ) const
      {
        return !operator==( b );
      }

      std::size_t size() const
      {
        return vt_.size();
      }

      const VT& vt() const
      {
        return vt_;
      }

    private:
      VT vt_;
    };

  }  // end namespace Impl
} // end namespace Dune

namespace std
{
  //! overload std::hash
  template <> struct hash<Dune::Impl::MultiId>
  {
    size_t operator()(const Dune::Impl::MultiId& x) const
    {
      static constexpr Dune::HashUIntVector hashUIntVector;
      return hashUIntVector( x.vt() );
    }
  };

  //! overload operator<<
  ostream& operator<<(ostream& os, const Dune::Impl::MultiId& multiId)
  {
    for( const auto& v : multiId.vt() )
      os << v << " ";
    return os;
  }
}

#endif
