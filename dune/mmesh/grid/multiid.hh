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
  namespace MMeshImpl
  {

    class MultiId
    {
    public:
      using ThisType = MultiId;
      using T = std::size_t;
      using VT = std::vector<T>;

      // we use an array as internal storage to make MultiId trivially copyable
      using Storage = std::array<T, 4>;

      MultiId() : size_(0), hash_(-1) {}

      MultiId( const MultiId& other )
       : vt_ ( other.vt_ ), size_( other.size_ ), hash_( other.hash_ )
      {}

      MultiId( const std::vector<T>& vt )
       : size_( vt.size() ), hash_(-1)
      {
        int i = 0;
        for (const auto& v : vt)
          vt_[ i++ ] = v;
      }

      MultiId( std::initializer_list<T> l )
       : MultiId ( std::vector<T>(l) )
      {}

      MultiId( T t )
       : MultiId ( { t } )
      {}

      ThisType& operator= ( const ThisType& b )
      {
        if (this != &b)
        {
          vt_ = b.vt_;
          size_ = b.size_;
          hash_ = b.hash_;
        }
        return *this;
      }

      bool operator< ( const ThisType& b ) const
      {
        if( size() != b.size() )
          return size() < b.size();

        for( int i = 0; i < size_; ++i )
          if( vt_[i] != b.vt_[i] )
            return vt_[i] < b.vt_[i];

        return false;
      }

      bool operator== ( const ThisType& b ) const
      {
        if( size() != b.size() )
          return false;

        for( int i = 0; i < size_; ++i )
          if( vt_[i] != b.vt_[i] )
            return false;

        return true;
      }

      bool operator<= ( const ThisType& b ) const
      {
        return !b.operator<(*this);
      }

      bool operator!= ( const ThisType& b ) const
      {
        return !operator==( b );
      }

      std::size_t size() const
      {
        return size_;
      }

      VT vt() const
      {
        std::vector<T> vec;
        vec.reserve(size_);
        for(int i = 0; i < size_; ++i)
          vec.push_back(vt_[i]);
        return vec;
      }

      //! Hash function with caching
      std::size_t hash() const
      {
        if (hash_ == std::size_t(-1))
        {
          if( size_ == 0 )
          {
            hash_ = 0;
            return hash_;
          }

          static constexpr std::hash<std::size_t> hasher;
          hash_ = hasher(vt_[0]);
          for ( std::size_t i = 1; i < size_; ++i )
            hash_ = hash_ ^ (hasher(vt_[i]) << i);
        }
        return hash_;
      }

    private:
      Storage vt_;
      std::size_t size_;
      mutable std::size_t hash_;
    };

  }  // end namespace MMeshImpl
} // end namespace Dune

namespace std
{
  //! overload std::hash
  template <> struct hash<Dune::MMeshImpl::MultiId>
  {
    size_t operator()(const Dune::MMeshImpl::MultiId& id) const
    {
      return id.hash();
    }
  };

  //! overload operator<<
  inline ostream& operator<<(ostream& os, const Dune::MMeshImpl::MultiId& multiId)
  {
    for( const auto& v : multiId.vt() )
      os << v << " ";
    return os;
  }
}

#endif
