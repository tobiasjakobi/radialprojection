#ifndef _VISIBILITY_H_
#define _VISIBILITY_H_

#include "pooled_alloc.h"

#include "common.h"

namespace VisTest {

  template <typename T>
  class VisibleList {
  private:
    typedef typename T::invectype _invectype;

    typedef typename PooledList<_invectype>::Type _vlist;
    typedef typename _vlist::iterator _iter;
    typedef typename _vlist::const_iterator _citer;

    OneTimePool* mempool;

    _vlist* l; // internal point list (double-linked)
    _iter c; // last insert position

  public:

    VisibleList() : mempool(0), l(0) {}

    VisibleList(const VisibleList& vl) : mempool(0), l(0) {
      if (vl.mempool != 0) {
        mempool = new OneTimePool(vl.mempool->getByteSize());
      }

      if (vl.l != 0) {
        l = new PooledList<vec4s>::Type(mempool);
        c = l->end();
      }
    }

    ~VisibleList() {
      cerr << "info: before destructing the pool used " << mempool->getBytesUsed()
           << " out of " << mempool->getByteSize() << " bytes ("
           << (double(mempool->getBytesUsed()) / double(mempool->getByteSize()) * 100.0)
           << "%).\n";

      delete l;
      delete mempool;
    }

    void reserve(uint elements) {
      if (mempool != 0) return;

      const size_t bytes = elements * PooledList<_invectype>::NodeByteSize;
      cerr << "info: reserving memory pool of " << bytes << " bytes.\n";

      mempool = new OneTimePool(bytes);
    }

    void init() {
      if (l != 0) return;

      l = new typename PooledList<_invectype>::Type(mempool);
      c = l->end();
    }

    void insertSorted(const _invectype& v) {
      const double v_a = T::angle(v);
      const double c_a = T::angle(*c);

      if (l->empty()) {
        insert(v, l->end());
        return;
      }

      bool found = false;
      _iter p = c;

      if (v_a == c_a) {
        found = isearch(v);
      } else {
        if (v_a > c_a) {
          found = fsearch(v, p);
        } else {
          // v_a < c_a
          found = bsearch(v, p);
        }
      }

      if (!found) {
        insert(v, p);
      } else {
        c = p;
      }
    }

    // collect nodes that are "near" i
    void collectNearNodes(const _iter& i, vector<_iter>& list) {
      const double i_a = T::angle(*i);

      _iter j;

      // collect in forward direction
      j = i;
      ++j;
      while (j != l->end()) {
        if (T::angle(*j) > i_a + T::epsilon)
          break;

        list.push_back(j);
        ++j;
      }

      // collect in backward direction
      j = i;
      while (j != l->begin()) {
        --j;

        if (T::angle(*j) < i_a - T::epsilon)
          break;

        list.push_back(j);
      }
    }

    /* This DOESN'T remove invisible points but removes all but one vertex  *
     * on a ray (through zero, the origin) containing multiple vertices     *
     * from the patch. So this produces a vertex set where each vertex      *
     * is visible from the origin.                                          *
     * This suffices for the radial projection approach!                    */
    void removeInvisibleFast() {
      cerr << "info: computing (incorrect) visibility for "
           << l->size() << " vertices." << endl;

      vector<_iter> nodes;

      for (_iter i = l->begin(); i != l->end(); ++i) {
        collectNearNodes(i, nodes);

        for (typename vector<_iter>::iterator j = nodes.begin(); j != nodes.end(); ++j) {
          const _iter k = *j;
          if (T::rayTest(*i, *k)) l->erase(k);
        }

        nodes.clear();
      }

      cerr << "info: " << l->size() << " vertices are visible." << endl;
    }

    // This method really does what the name implies!
    // It only removes but the vertex being nearest to the origin.
    void removeInvisibleProper() {
      cerr << "info: computing (proper) visibility for "
           << l->size() << " vertices." << endl;

      vector<_iter> nodes;

      for (_iter i = l->begin(); i != l->end(); ++i) {
        collectNearNodes(i, nodes);

        for (typename vector<_iter>::iterator j = nodes.begin(); j != nodes.end(); ++j) {
          const _iter k = *j;
          if (!T::rayTest(*i, *k)) continue;

          if (T::toR2(*k).lengthSquared() >= T::toR2(*i).lengthSquared())
            l->erase(k);
        }

        nodes.clear();
      }

      cerr << "info: " << l->size() << " vertices are visible." << endl;
    }

    void toR2(Common::vec2dlist& output) const {
      for (_citer i = l->begin(); i != l->end(); ++i) {
        output.push_back(T::toR2(*i));
      }
    }

    void dump(vector<_invectype>& output) const {
      for (_citer i = l->begin(); i != l->end(); ++i) {
        output.push_back(*i);
      }
    }

    uint size() const {
      return l->size();
    }

    bool isInitialized() const {
      return l != 0;
    }

    bool isOrdered() const {
      if (l->size() == 0 || l->size() == 1) {
        return true;
      }

      _citer i = l->begin();

      bool ordered = true;
      double a = T::angle(*i);

      ++i;
      while (i != l->end()) {
        const double b = T::angle(*i);

        if (a > b) {
          ordered = false;
          break;
        }

        a = b;
        ++i;
      }

      return ordered;
    }

    void print() const {
      if (l->size() == 0) return;

      _citer i = l->begin();
      cerr << T::angle(*i);;
      ++i;

      while (i != l->end()) {
        cerr << ' ' << T::angle(*i);
        ++i;
      }
      cerr << endl;
    }

  private:

    bool isearch(const _invectype& v) {
      const double v_a = T::angle(v);
      _iter i;

      i = c;
      while (i != l->end()) {
        if (T::angle(*i) != v_a)
          break;

        if (v == *i)
          return true;

        ++i;
      }

      i = c;
      while (i != l->begin()) {
        --i;

        if (T::angle(*i) != v_a)
          break;

        if (v == *i)
          return true;
      }

      return false;
    }

    // search in forward direction for the entry v:
    // returns true if entry found
    // returns false otherwise together with an insert position
    bool fsearch(const _invectype& v, _iter& inspos) {
      const double v_a = T::angle(v);

      _iter i = c;

      while (i != l->end()) {
        if (v == *i) {
          return true;
        }

        if (T::angle(*i) > v_a) {
          inspos = i;
          return false;
        }

        ++i;
      }

      // reached the end of the list without finding anything
      // and v_a is greater than all entries in the list

      inspos = l->end();

      return false;
    }

    bool bsearch(const _invectype& v, _iter& inspos) {
      const double v_a = T::angle(v);

      _iter i = c;
      ++i;

      while (i != l->begin()) {
        --i;

        if (v == *i) {
          return true;
        }

        if (T::angle(*i) < v_a) {
          inspos = i;
          ++inspos;
          return false;
        }
      }

      // reached the start of the list without finding anything
      // and v_a is smaller than all entries in the list

      inspos = l->begin();

      return false;
    }

    void insert(const _invectype& v, const _iter& pos) {
      c = l->insert(pos, v);
    }

  };

};

#endif /* _VISIBILITY_H_ */

