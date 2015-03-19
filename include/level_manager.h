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
  uint lbegin;
  uint lend;

  inline uint range() const {
    return (lend - lbegin);
  }

public:
  TVLevel() {
    lbegin = 0;
    lend = 0;
  }

  TVLevel(const TVLevel<T>& tvl) {
    lbegin = tvl.lbegin;
    lend = tvl.lend;
  }

  void init(uint a, uint b) {
    assert(b >= a);

    lbegin = a;
    lend = b;
  }

  void assign(const TVLevel& tvl) {
    lbegin = tvl.lbegin;
    lend = tvl.lend;
  }

  void insert(vector<T>& list, const T& item) {
    ++lend;
    list.push_back(item);
  }

  void shift(const TVLevel& tvl) {
    lbegin = tvl.lend;
    lend = tvl.lend;
  }

  inline uint begin() const {
    return lbegin;
  }

  inline uint end() const {
    return lend;
  }

  bool locate(const vector<T>& list, const T& target) const {
    const uint r = this->range();
    if (r == 0) return false;

    for (uint k = 0; k < r; ++k) {
      if (list.at(lbegin + k) == target) return true;
    }

    return false;
  }

};

/* Tiling Vertex (list) Level Manager:                                *
 * We assume that the list will always contain one (initial) element. */
template <typename T, unsigned int N>
class TVLManager {
private:
  vector<T>& list;
  TVLevel<T> levels[N];

public:
  TVLManager(vector<T>& l) : list(l) {
    assert(N >= 2);

    // We use a different way to count levels here (insertion level is 0)
    levels[0].init(1, 1);
    levels[1].init(0, 1);
  }

  // Forbid to use copy-constructor
  TVLManager(const TVLManager<T, N>& tvlm) {
    assert(false);
  }

  bool insert(const T& item) {
    for (uint i = 0; i < N; ++i) {
      if (levels[i].locate(list, item)) return false;
    }

    levels[0].insert(list, item);
    return true;
  }

  void advance() {
    for (uint i = N - 1; i > 0; --i) {
      levels[i].assign(levels[i - 1]);
    }

    levels[0].shift(levels[1]);
  }  

  uint begin() const {
    return levels[1].begin();
  }

  uint end() const {
    return levels[1].end();
  }

};

#endif /* _LEVEL_MANAGER_H_ */

