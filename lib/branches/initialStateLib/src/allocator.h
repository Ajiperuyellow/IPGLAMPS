//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/allocator.h $
//$LastChangedDate: 2015-02-26 17:06:51 +0100 (Do, 26. Feb 2015) $
//$LastChangedRevision: 2093 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


/** @file 
 * @brief Modifications of memory allocation
 *
 * Here we collect stuff connencted with memory allocations,
 * especially concerned with alignment neccessary for SSE (i.e. at 16
 * byte boundaries)
 *
 * Most of the following code is 'stolen' (and simplified) from 
 * the Vc library by Matthias Kretz. 
 */

#ifndef ALLOCATOR_H
#define ALLOCATOR_H

#include "configBAMPS.h"

#if defined(VECTOR_IMPL_SSE) || defined(VECTOR_IMPL_AVX)

#include <cstddef>
#include <new>

#define BAMPS_DECLARE_ALLOCATOR(Type)                                   \
  namespace std                                                         \
  {                                                                     \
    template<> class allocator<Type> : public ::Allocator<Type>         \
    {                                                                   \
    public:                                                             \
      template<typename U> struct rebind { typedef ::std::allocator<U> other; }; \
    };                                                                  \
  }
  
#define ALIGNOF(_TYPE_) __alignof(_TYPE_)

template<typename T> class Allocator
{
private:
  enum Constants {
    NaturalAlignment = sizeof(void *) > ALIGNOF(long double) ? sizeof(void *) :
    (ALIGNOF(long double) > ALIGNOF(long long) ? ALIGNOF(long double) : ALIGNOF(long long)),
    Alignment = ALIGNOF(T),
    ExtraBytes = Alignment > NaturalAlignment ? Alignment : 0,
    AlignmentMask = Alignment - 1
  };
    
public:
  typedef std::size_t    size_type;
  typedef std::ptrdiff_t difference_type;
  typedef T*             pointer;
  typedef const T*       const_pointer;
  typedef T&             reference;
  typedef const T&       const_reference;
  typedef T              value_type;
    
  template<typename U> struct rebind { typedef Allocator<U> other; };
    
  Allocator() throw() { }
  Allocator(const Allocator&) throw() { }
  template<typename U> Allocator(const Allocator<U>&) throw() { }
    
  pointer address(reference x) const { return &x; }
  const_pointer address(const_reference x) const { return &x; }
    
  pointer allocate(size_type n, const void* = 0)
  {
    if (n > this->max_size()) {
      throw std::bad_alloc();
    }
      
    char *p = static_cast<char *>(::operator new(n * sizeof(T) + ExtraBytes));
    if (ExtraBytes > 0) {
      char *const pp = p;
      p += ExtraBytes;
      const char *null = 0;
      p -= ((p - null) & AlignmentMask); // equivalent to p &= ~AlignmentMask;
      reinterpret_cast<char **>(p)[-1] = pp;
    }
    return reinterpret_cast<pointer>(p);
  }
    
  void deallocate(pointer p, size_type)
  {
    if (ExtraBytes > 0) {
      p = reinterpret_cast<pointer *>(p)[-1];
    }
    ::operator delete(p);
  }
    
  size_type max_size() const throw() { return size_t(-1) / sizeof(T); }
    
  void construct(pointer p, const T& __val) { ::new(p) T(__val); }
  void destroy(pointer p) { p->~T(); }
};

template<typename T> inline bool operator==(const Allocator<T>&, const Allocator<T>&) { return true;  }
template<typename T> inline bool operator!=(const Allocator<T>&, const Allocator<T>&) { return false; }

#if defined(VECTOR_IMPL_SSE) 
#define VectorAlignment 16u
#elif defined(VECTOR_IMPL_AVX)
#define VectorAlignment 32u
#endif



#define FREE_STORE_OPERATORS_ALIGNED                         \
  inline __attribute__((__always_inline__)) void *operator new(std::size_t size) { return _mm_malloc(size, VectorAlignment); } \
  inline __attribute__((__always_inline__)) void *operator new(std::size_t, void *p) { return p; } \
  inline __attribute__((__always_inline__)) void *operator new[](std::size_t size) { return _mm_malloc(size, VectorAlignment); } \
  inline __attribute__((__always_inline__)) void *operator new[](std::size_t , void *p) { return p; } \
  inline __attribute__((__always_inline__)) void operator delete(void *ptr, std::size_t) { _mm_free(ptr); } \
  inline __attribute__((__always_inline__)) void operator delete(void *, void *) {} \
  inline __attribute__((__always_inline__)) void operator delete[](void *ptr, std::size_t) { _mm_free(ptr); } \
  inline __attribute__((__always_inline__)) void operator delete[](void *, void *) {}


#else

// Do nothing:

#define VectorAlignment 1u
#define BAMPS_DECLARE_ALLOCATOR(Type)
#define FREE_STORE_OPERATORS_ALIGNED

#endif

#endif
