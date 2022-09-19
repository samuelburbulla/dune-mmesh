// This implementation has been taken from Dune-ALUGrid!
// (dune/alugrid/impl/serial/serialize.h)
// Credit: Bernhard Schupp, 1997-1998 and Robert Kloefkorn 2006

#ifndef DUNE_MMESH_MISC_OBJECTSTREAM_HH
#define DUNE_MMESH_MISC_OBJECTSTREAM_HH

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <iostream>
#include <utility>

namespace Dune
{

namespace MMeshImpl
{

  // ObjectStream
  // ------------

  struct ObjectStreamTraits
  {
    template< class T >
    static void copy ( T *dest, const void *src, std::size_t n )
    {
      for( std::size_t i = 0; i < n; ++i )
        dest[ i ] = static_cast< const T * >( src )[ i ];
    }

    template< class T >
    static void copy ( void *dest, const T *src, std::size_t n )
    {
      for( std::size_t i = 0; i < n; ++i )
        static_cast< T * >( dest )[ i ] = src[ i ];
    }
  };

  class ObjectStream
  {
    using Traits = ObjectStreamTraits;
  public:
    char *_buf;
    size_t _rb, _wb, _len;
  protected:
    const size_t _bufChunk;
    mutable bool _owner;

  public :
    class EOFException
    {
    public:
      virtual std::string what () const { return "EOFException"; }
    };

    class OutOfMemoryException {};

    inline ObjectStream (size_t chunk = 0)
      : _buf(0), _rb(0), _wb(0), _len(0), _bufChunk(chunk), _owner(true)
    {}

    inline ObjectStream (const ObjectStream & os)
      : _buf(0), _rb(0), _wb(0), _len(0), _bufChunk(os._bufChunk), _owner(true)
    {
      assign(os);
    }

    // reset write and read postitions
    inline void clear() { _wb = 0; _rb = 0; }

    // reset read position
    inline void resetReadPosition() { _rb = 0; }

    //! set position of write counter
    void seekp( const size_t pos )
    {
      _wb = pos;
      assert ( _wb <= _len );
    }

    // return's true if size > 0 and read position is zero
    // i.e. a read othe stream will result some valid data
    inline bool validToRead () const { return (_wb > 0) && (_rb == 0); }

    // return size of bytes allready written to stream
    inline int capacity() const { return _len; }

    // return size of bytes allready written to stream
    inline int size() const { return _wb; }

    // make sure that s bytes memory can be wrote without reallocation
    inline void reserve(size_t s)
    {
      const size_t newSize = _wb + s;
      if (newSize > _len) reallocateBuffer( newSize );
    }

    // delete stream
    inline ~ObjectStream () { removeObj(); }

    //! assign buffer from os the local buffer, os ownership is false afterwards
    //inline const ObjectStream & operator = (const ObjectStream & os)
    inline ObjectStream & operator = (const ObjectStream & os)
    {
      removeObj();
      assign(os);
      return *this;
    }

    // write value to stream
    template <class T>
    inline void write (const T & a)
    {
      writeT( a, true );
    }

    template <class T>
    inline void writeUnchecked( const T& a )
    {
      writeT( a, false );
    }

    ////////////////////////////////////
    // to behave like stringstream
    ////////////////////////////////////

    // put char
    inline void put (const signed char a)  { write(a); }

    // put char with checking buffer size (reserve size before usage)
    inline void putNoChk (const signed char a)  { writeUnchecked(a); }

    // get char
    inline signed char get ()
    {
      signed char a;
      read(a);
      return a;
    }

    // eof function
    bool eof () const { return (this->_rb >= this->_wb); }

    // good function
    bool good () const { return (this->_rb < this->_wb); }

    /////////////////////////////////////

  protected:
    template <class T>
    inline void writeT (const T & a, const bool checkLength )
    {
      assert ( _owner );
      const size_t ap = _wb;
      _wb += sizeof(T);

      // if buffer is to small, reallocate
      if (checkLength && _wb > _len)
      {
        reallocateBuffer(_wb);
      }
      assert ( _wb <= _len );

      // call assignment operator of type T
      Traits::copy( static_cast< void * >( getBuff( ap ) ), &a, 1 );
      return;
    }

    template<class T>
    inline void readT ( T& a, bool checkLength )
    {
      const size_t ap = _rb;
      _rb += sizeof(T);
      assert ( _rb <= _wb );

      // call assignment operator of type T
      Traits::copy( &a, static_cast< const void * >( getBuff( ap ) ), 1 );
      return;
    }

  public:
    // read value from stream
    template <class T>
    inline void read (T & a) { readT( a, true ); }

    template<class T>
    inline void readUnchecked ( T& a ) { readT( a, false ); }

    // read this stream and write to os
    inline void readStream (ObjectStream & os)
    {
      readStream(os,_wb);
    }

    // read length bytes from this stream and stores it to os
    inline void readStream (ObjectStream & os, const size_t length)
    {
      if( length == 0 ) return;
      // actual read position
      os.write( getBuff(_rb), length);
      removeObject(length);
    }

    // writes hole stream of os to this stream
    inline void writeStream (const ObjectStream & os)
    {
      write(os._buf, os._wb);
    }

    // increments the read position without actualy read data
    inline void removeObject(const size_t length)
    {
      _rb += length;
      assert ( _rb <= _wb );
    }

    //! free allocated memory
    inline void reset()
    {
      removeObj();
    }

    // static free for use with all buffers here
    inline static void freeBuffer(char * buffer)
    {
      // free buffer if not zero
      if( buffer ) free( buffer );
    }

    // compatibility with ostream
    inline void write(const char* buff, const size_t length )
    {
      assert ( _owner );
      if( length == 0 ) return;

      const size_t newWb = _wb + length;
      if (newWb > _len) reallocateBuffer(newWb);

      memcpy( getBuff(_wb), buff, length );
      _wb = newWb;
    }

    // compatibility with istream
    inline void read(char* buff, const size_t length )
    {
      if( length == 0 ) return;

      const size_t newRb = _rb + length;
      assert ( newRb <= _wb );

      memcpy( buff, getBuff(_rb), length );
      _rb = newRb;
    }

    // return pointer to buffer memory
    inline char* raw () { return getBuff( 0 ); }
    inline const char* raw () const { return getBuff( 0 ); }

    inline char * getBuff (const size_t ap) { return (_buf + ap); }
    inline const char * getBuff (const size_t ap) const { return (_buf + ap); }

  protected:
    // reallocated the buffer if necessary
    inline void reallocateBuffer(size_t newSize)
    {
      assert ( _owner );
      _len += _bufChunk;
      if(_len < newSize) _len = newSize;
      _buf = (char *) realloc (_buf, _len);
      if (!_buf)
      {
        perror ("**EXCEPTION in ObjectStream :: reallocateBuffer(size_t) ");
        throw OutOfMemoryException ();
      }
    }

    // delete buffer
    inline void removeObj()
    {
      if( _owner ) freeBuffer( _buf );
      _buf = 0; _len = 0; _wb = 0; _rb = 0; _owner = true;
      return;
    }

    // assign buffer
    inline void assign(const ObjectStream & os)
    {
      assert ( _buf == 0 );
      if( os._len > 0 )
      {
        _len = os._len;
        _wb  = os._wb;
        _rb  = os._rb;
        const_cast<size_t &> (_bufChunk) = os._bufChunk;

        // overtake buffer and set ownership of os to false
        _buf = os._buf;
        os._owner = false;
        // we are owner now
        _owner = true;
      }
      return;
    }

    inline void assign(char * buff, const size_t length)
    {
      if( length == 0 ) return;

      // if length > 0, buff should be valid
      assert ( buff );

      // set length
      _wb = _len = length;
      // set buffer
      _buf = buff;

      // read status is zero
      _rb = 0;

      // we are the owner
      _owner = true;
      return;
    }
  };

} // namespace MMeshImpl

} // namespace Dune

#endif // #ifndef DUNE_MMESH_MISC_OBJECTSTREAM_HH
