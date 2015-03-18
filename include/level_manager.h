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

#ifndef _LEVEL_MANAGER_H_
#define _LEVEL_MANAGER_H_

#include "common.h"

/* TVLevel represents a "level" in a _T_iling _V_ertices list.         *
 * Used for constructing a tiling patch via the model set description. */
template <typename T>
class TVLevel {
private:
  vector<T>* vlist; // list of vertices
  //vector<T> vlist;

  //vector<T>::const_iterator begin;
  //vector<T>::const_iterator end;

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
template <typename T, unsigned int N>
class TVLManager {
private:
  TVLevel<T> levels[N];

  // We use a different way to count levels here (insertion level is 0)
  void init(vector<T>& l) {
    for (uint i = 0; i < N; ++i) {
      levels[i].bind(l);
    }

    levels[0].init(2, 1);
    levels[1].init(1, 1);
  }

public:
  TVLManager(vector<T>& list) {
    assert(N >= 2);
    init(list);
  }

  // Forbid to use copy-constructor
  TVLManager(const TVLManager<T, N>& tvlm) {
    assert(false);
  }

  bool insert(const T& item) {
    for (uint i = 0; i < N; ++i) {
      if (levels[i].locate(item)) return false;
    }

    levels[0].insert(item);
    return true;
  }

  void advance() {
    for (uint i = N - 1; i > 0; --i) {
      levels[i].assign(levels[i - 1]);
    }

    levels[0].shift(levels[1]);
  }  

  uint begin() const {
    // Compensate for the 1-shift
    return (levels[1].begin() - 1);
  }

  uint end() const {
    return levels[1].end();
  }

};

#endif /* _LEVEL_MANAGER_H_ */

