/*  radialprojection - tools to numerically compute the radial projection of point sets
 *  Copyright (C) 2012-2014 - Tobias Jakobi <tjakobi at math dot uni dash bielefeld dot de>
 *
 *  radialprojection is free software: you can redistribute it and/or modify it under the terms
 *  of the GNU General Public License as published by the Free Software Foundation, either
 *  version 3 of the License, or (at your option) any later version.
 *
 *  radialprojection is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 *  PURPOSE. See the GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along with radialprojection.
 *  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _POOLED_ALLOCATOR_H_
#define _POOLED_ALLOCATOR_H_

#include <cassert>
#include <cstddef>
#include <list>

class OneTimePool
{
private:
  size_t m_size;
  size_t m_used;

  char* m_storage;
  char* m_curpos;

  typedef void* voidptr;

  bool isFromPool(void const* instance) const {
    char const* block = reinterpret_cast<char const*>(instance);
    return m_storage <= block && block < m_curpos;
  }

public:
  OneTimePool(size_t size) : m_size(size), m_used(0),
                             m_storage(0), m_curpos(0) {
    if (m_size > 0) {
      m_storage = new char[m_size];
      m_curpos = m_storage;
    }
  }

  ~OneTimePool() {
    assert(m_used == 0 && "can't destroy a pool with outstanding allocations");

    delete [] m_storage;
  }

  size_t getByteSize() const {
    return m_size;
  }

  size_t getBytesUsed() const {
    return m_used;
  }

  bool isFree() const {
    return m_used == 0;
  }

  bool alloc(voidptr& block, size_t bytes) {
    if (m_curpos + bytes < m_storage + m_size) {
      block = m_curpos;
      m_curpos += bytes;
      m_used += bytes;
      return true;
    }

    return false;
  }

  bool dealloc(voidptr block, size_t bytes) {
    if (isFromPool(block)) {
      m_used -= bytes;
      return true;
    }

    return false;
  }

};

//! A standards-compliant pooled allocator.
template<typename T>
class PooledAllocator
{
public:
  typedef size_t      size_type;
  typedef ptrdiff_t   difference_type;

  typedef T           value_type;
  typedef T*          pointer;
  typedef T const*    const_pointer;
  typedef T&          reference;
  typedef T const&    const_reference;

  template<typename U> 
  struct rebind {
    typedef PooledAllocator<U> other;
  };

  PooledAllocator(OneTimePool* pool = 0) : m_pool(pool) {}

  template<typename U>
  PooledAllocator(PooledAllocator<U> const& arg) : m_pool(arg.m_pool) {}

  size_type max_size() const {
    return 0xffffffff;
  }

  pointer allocate(size_type count, std::allocator<void>::const_pointer = 0) const {
    void* memory;

    if (m_pool && m_pool->alloc(memory, count * sizeof(T))) {
      return reinterpret_cast<T*>(memory);
    } else {
      return reinterpret_cast<T*>(new char[count * sizeof(T)]);
    }
  }

  void deallocate(pointer block, size_type count) const throw() {
    assert(block && "null pointer argument");

    if(m_pool && m_pool->dealloc(block, count * sizeof(T))) {
      return;
    } else {
      delete [] reinterpret_cast<char*>(block);
    }
  }

  void construct(pointer element, T const& arg) {
    new(element) T(arg);
  }

  void destroy(pointer element) {
    element->~T();
  }

  pointer address(reference element) const {
    return &element;
  }

  const_pointer address(const_reference element) const {
    return &element;
  }

  OneTimePool* m_pool;
};

template<>
class PooledAllocator<void>
{
public:
  typedef size_t       size_type;
  typedef ptrdiff_t    difference_type;

  typedef void         value_type;
  typedef void*        pointer;
  typedef void const*  const_pointer;

  template<typename U> 
  struct rebind {
    typedef PooledAllocator<U> other;
  };

  PooledAllocator() : m_pool(0) {}

  PooledAllocator(OneTimePool* pool) : m_pool(pool) {}

  OneTimePool* m_pool;
};

template<typename T, typename U>
bool operator==(PooledAllocator<T> const& left, PooledAllocator<U> const& right) {
  return left.m_pool == right.m_pool;
}

template<typename T, typename U>
bool operator!=(PooledAllocator<T> const& left, PooledAllocator<U> const& right) {
  return left.m_pool != right.m_pool;
}

template<typename Value>
struct PooledList { 
  typedef std::list<Value, PooledAllocator<Value> > Type;
  static const size_t NodeByteSize = sizeof(std::_List_node<Value>); // WARNING: GCC specific!
};

#endif // _POOLED_ALLOCATOR_H_

