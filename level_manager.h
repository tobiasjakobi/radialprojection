#ifndef _LEVEL_MANAGER_H_
#define _LEVEL_MANAGER_H_

#include "common.h"

/* TVLevel represents a "level" in a _T_iling _V_ertices list.         *
 * Used for constructing a tiling patch via the model set description. */
template <typename T>
class TVLevel {
private:
  vector<T>* vlist; // list of vertices

  uint lbegin;
  uint lend;

  uint range() const {
    return (lend + 1 - lbegin);
  }

  bool empty() const {
    return (lend < lbegin);
  }

public:
  TVLevel() : vlist(0) {
    lbegin = 1;
    lend = 0;
  }

  TVLevel(const TVLevel<T>& tvl) : vlist(tvl.vlist) {
    lbegin = tvl.lbegin;
    lend = tvl.lend;
  }

  ~TVLevel() {}

  void bind(vector<T>& list) {
    vlist = &list;
  }

  void init(uint a, uint b) {
    assert(b + 1 >= a);

    lbegin = a;
    lend = b;
  }

  void assign(const TVLevel& tvl) {
    lbegin = tvl.lbegin;
    lend = tvl.lend;
  }

  void insert(const T& item) {
    assert(vlist != 0);

    ++lend;
    vlist->push_back(item);
  }

  void shift(const TVLevel& tvl) {
    lbegin = tvl.lend + 1;
    lend = tvl.lend;
  }

  uint begin() const {
    return lbegin;
  }

  uint end() const {
    return lend;
  }

  bool locate(const T& target) const {
    if (this->empty()) return false;
    const uint r = this->range();

    for (uint k = 0; k < r; ++k) {
      if (vlist->at(lbegin + k - 1) == target) return true;
    }

    return false;
  }

};

/* Tiling Vertex (list) Level Manager:                                *
 * We assume that the list will always contain one (initial) element. */
template <typename T>
class TVLManager {
private:
  uint num_levels;
  TVLevel<T>* levels;

  // We use a different way to count levels here (insertion level is 0)
  void init(vector<T>& l) {
    assert(levels == 0);
    levels = new TVLevel<T>[num_levels];

    for (uint i = 0; i < num_levels; ++i) {
      levels[i].bind(l);
    }

    levels[0].init(2, 1);
    levels[1].init(1, 1);
  }

public:
  TVLManager(uint n, vector<T>& list) : num_levels(n), levels(0) {
    assert(num_levels >= 2);
    init(list);
  }

  // Forbid to use copy-constructor
  TVLManager(const TVLManager<T>& tvlm) {
    assert(false);
  }

  ~TVLManager() {
    delete [] levels;
  }

  bool insert(const T& item) {
    assert(levels != 0);

    for (uint i = 0; i < num_levels; ++i) {
      if (levels[i].locate(item)) return false;
    }

    levels[0].insert(item);
    return true;
  }

  void advance() {
    for (uint i = num_levels - 1; i > 0; --i) {
      levels[i].assign(levels[i - 1]);
    }

    levels[0].shift(levels[1]);
  }  

  uint begin() const {
    assert(levels != 0);
    // Compensate for the 1-shift
    return (levels[1].begin() - 1);
  }

  uint end() const {
    assert(levels != 0);
    return levels[1].end();
  }

};

#endif /* _LEVEL_MANAGER_H_ */

