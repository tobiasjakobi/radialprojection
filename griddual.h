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

#ifndef _GRIDDUAL_H_
#define _GRIDDUAL_H_

#include "common.h"

#include <limits>
#include <fstream>

namespace GridDualizing {

  typedef vector<vec2d> vec2dlist;

  /* n = degree of rotational symmetry                                          *
   * This boils down to the Z-module over which all calculations are done.      *
   * In general the module is Z[xi_n], so without any simplifications           *
   * we need n integers to encode an element. But if n is even, then we         *
   * can reduce this amount to n/2 integers, since the latter half of the       *
   * xi^k is just the first half with a sign switch.                            *
   *                                                                            *
   * One can show that one only needs eulerphi(n) components to encode          *
   * an element, but for general n this reduction looks non-trivial, and        *
   * it also interfers with the main grid algorithm .                           */
  template <uint n>
  class GridVertex {
  public:

    // Store only first half of components when n is even.
    enum { size = ((n % 2 == 0) ? n/2 : n),
           realsize = n,
           reducedsize = Meta::eulerPhi<n>::value };

    GridVertex() {}
    GridVertex(const GridVertex<n>& gv) {
      for (uint i = 0; i < size; ++i) {coords[i] = gv.coords[i];}
    }

    GridVertex(int x) {
      for (uint i = 0; i < size; ++i) {coords[i] = x;}
    }

    GridVertex(int c[]) {
      assert(c != NULL);
      for (uint i = 0; i < size; ++i) {coords[i] = c[i];}
    }

    void clip(const GridVertex<n>& min, const GridVertex<n>& max) {
      for (uint i = 0; i < size; ++i) {
        coords[i] = ::max(::min(coords[i], max[i] + 1), min[i]);
      }
    }

    bool operator==(const GridVertex<n>& v) const {
      for (uint i = 0; i < size; ++i) {
        if (coords[i] != v.coords[i]) return false;
      }
      return true;
    }

    int operator[](uint i) const {
      assert(i < size);
      return coords[i];
    }

    int& operator[](uint i) {
      assert(i < size);
      return coords[i];
    }

    GridVertex<n>& operator=(const GridVertex<n>& gv) {
      for (uint i = 0; i < size; ++i) {coords[i] = gv.coords[i];}
      return *this;
    }

    GridVertex<n> operator+(const GridVertex<n>& gv) const {
      int newcoords[size];
      for (uint i = 0; i < size; ++i) {newcoords[i] = coords[i] + gv.coords[i];}
      return GridVertex<n>(newcoords);
    }

    vec2d to2D() const {
      assert(roots != NULL);
      vec2d res(roots[0] * double(coords[0]));
      for (uint i = 1; i < size; ++i) {res += roots[i] * double(coords[i]);}
      return res;
    }

    /* The internal representation (see also above) for a grid vertex is, even  *
     * after the n/2-"optimization", not unique. The components of the tuples   *
     * might differ, but the corresponding element is still the same. We can    *
     * achieve uniqueness by implementing the reduction to eulerphi(n)          *
     * components.                                                              *
     * Currently this is only done for n prime (where eulerphi(n) = n-1).       */
    GridVertex<reducedsize * 2> reduce() const {
      // TODO: Implement the other (non-prime) case!
      assert(size == reducedsize + 1);
      GridVertex<reducedsize * 2> res;

      for (uint i = 0; i < reducedsize; ++i) {
        res[i] = coords[i] - coords[reducedsize];
      }

      return res;
    }

    static void initRoots() {
      using namespace Common;

      if (roots != NULL) return;

      roots = new vec2d[size];

      for (uint i = 0; i < size; ++i) {
        roots[i].set(cos(2.0 * pi * double(i) / double(realsize)),
                     sin(2.0 * pi * double(i) / double(realsize)));
      }
    }

    static uint getSize() {
      return uint(size);
    } 

    static uint getReducedSize() {
      return uint(reducedsize);
    } 

  private:
    static vec2d* roots;
    int coords[size];
  };

  template <uint n>
  vec2d* GridVertex<n>::roots = NULL;

  template <typename T>
  class GridLine {
  public:
    typedef T vtype;

    GridLine() {}

    GridLine(const GridLine<T>& gl) {
      endpoints[0] = gl.endpoints[0];
      endpoints[1] = gl.endpoints[1];
    }

    GridLine(const vtype& a0, const vtype& a1) {
      endpoints[0] = a0;
      endpoints[1] = a1;
    }

    void set(const vtype& a0, const vtype& a1) {
      endpoints[0] = a0;
      endpoints[1] = a1;
    }

    bool operator==(const GridLine<T>& gl) const {
      return ((endpoints[0] == gl.endpoints[0] && endpoints[1] == gl.endpoints[1]) ||
              (endpoints[0] == gl.endpoints[1] && endpoints[1] == gl.endpoints[0]));
    }

    const vtype& getEndpoint(uint i) const {
      assert(i < 2);
      return endpoints[i];
    }

  private:
    vtype endpoints[2];
  };

  template <typename T>
  class GridTile {
  public:
    typedef T vtype;
    typedef GridLine<vtype> ltype; 

    GridTile() {}

    GridTile(const GridTile<T>& gt) {
      vertices[0] = gt.vertices[0];
      vertices[1] = gt.vertices[1];
      vertices[2] = gt.vertices[2];
      vertices[3] = gt.vertices[3];
    }

    GridTile(const vtype& t0, const vtype& t1,
             const vtype& t2, const vtype& t3) {
      vertices[0] = t0; vertices[1] = t1;
      vertices[2] = t2; vertices[3] = t3;
    }

    void clip(const vtype& min, const vtype& max) {
      vertices[0].clip(min, max);
      vertices[1].clip(min, max);
      vertices[2].clip(min, max);
      vertices[3].clip(min, max);
    }

    const vtype& getVertex(uint i) const {
      assert(i < 4);
      return vertices[i];
    }

    void getLines(ltype lines[4]) const {
      lines[0].set(vertices[0], vertices[1]);
      lines[1].set(vertices[1], vertices[2]);
      lines[2].set(vertices[2], vertices[3]);
      lines[3].set(vertices[3], vertices[0]);
    }

  private:
    vtype vertices[4];
  };

  template <typename T>
  struct GridTiling {
    typedef T vtype;
    typedef GridTile<T> ttype;
    typedef GridLine<T> ltype;

    vector<ttype> tiles;
    vector<vtype> vertices;
    vector<ltype> lines;
  };

  template <uint n>
  ostream& operator<<(ostream &os, const GridVertex<n>& gv) {
    os << '{' << gv[0];

    for (uint i = 1; i < GridVertex<n>::size; ++i) {
      os << ',' << gv[i];
    }

    os << '}';

    return os;
  }

  template <typename T>
  ostream& operator<<(ostream &os, const GridTile<T>& gt) {
    os << '{' << gt.getVertex(0) << ',' << gt.getVertex(1)
       << ',' << gt.getVertex(2) << ',' << gt.getVertex(3)
       << '}';

    return os;
  }

  template <typename T>
  ostream& operator<<(ostream &os, const GridLine<T>& gl) {
    os << '{' << gl.getEndpoint(0) << ','
       << gl.getEndpoint(1) << '}';

    return os;
  }

};

template <typename T>
ostream& operator<<(ostream &os, const GridDualizing::GridTiling<T>& tiling) {
  using namespace GridDualizing;

  os << '{' << GridTiling<T>::vtype::size << ','
     << tiling.vertices << ',' << tiling.lines << ','
     << tiling.tiles << '}';

  return os;
}

#endif

