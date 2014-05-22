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

#ifndef _HEXAGONAL_H_
#define _HEXAGONAL_H_

#include "common.h"

#include <algorithm>

namespace Triangular {

  const double radiusFactor = sqrt(3.0) * 0.5;

  void tiling(const vec2i& initpoint, uint maxstep,
              Common::vec2ilist& tilingpoints);

  void tilingVisLocal(const vec2i& initpoint, uint maxstep,
                      Common::vec2ilist& tilingpoints,
                      Common::vec2ilist& visiblepoints);

  void extractSector(const Common::vec2ilist& input,
                     Common::vec2ilist& output);

  void radialProj(const Common::vec2ilist& input,
                  Common::dlist& output, double& meandist);

};

/* Hexagonal tiling, also known as honeycomb structure. */
namespace Hexagonal {

  void tiling(const vec2i& initpoint, uint maxstep,
              Common::vec2ilist& tilingpoints);

};

// TODO: reimplement this in the style of cyclotomic_radial_xyz
class hexTiling {
public:
  typedef vec2i vertex;
  typedef vector<vertex> vlist;

private:

  struct level {
    uint a, b;
  };

  static bool locate(const vlist& a, const level& l, const vertex& v);

  void buildVertices(const vlist& a);

  vlist* vertices;
  bool onlyvis;
  bool onlysector;
  
public:
  hexTiling(uint N, const vertex& init, bool vis, bool sector);

  ~hexTiling() {
    delete vertices;
  }

  const vlist& getVertices() const {
    return *vertices;
  }
};

bool hexTiling::locate(const vlist& a, const level& l, const vertex& v) {
  if (l.a == l.b) return false;

  uint i = l.a;

  while (i != l.b) {
    if (a[i] == v) return true;

    ++i;
  }

  return false;
}

void hexTiling::buildVertices(const vlist& a) {
  assert(vertices == 0);

  cerr << "info: building vertex set from midpoints..." << endl;
  if (onlysector) {
    cerr << "info: restricting set to a 1/6-sector." << endl;
  }
  if (onlyvis) {
    cerr << "info: restricting to only visible vertices." << endl;
  }

  const uint numverts = 6;
  const vertex hexsteps[6] = {
    vertex(1, 0),  vertex(0, 1),  vertex(-1, 1),
    vertex(-1, 0), vertex(0, -1), vertex(1, -1)
  };

  vertices = new vlist;

  for (vlist::const_iterator i = a.begin(); i != a.end(); ++i) {
    for (uint j = 0; j < numverts; ++j) {
      const vertex t(*i + hexsteps[j]);

      // do sector test if enabled
      if (onlysector) {
        const vec2d x(t.transTriToR2());

        if (!x.inFirstQuadrant()) continue;
        if (!x.inSectorL3()) continue;
      }

      // do visibility test if enabled
      if (onlyvis && Coprime::gcdZ(abs(t.x), abs(t.y)) != 1)
        continue;

      // backward search
      // TODO: search can be optimized
      if (find(vertices->rbegin(), vertices->rend(), t) == vertices->rend()) {
        vertices->push_back(t);
      }
    }
  }

  cerr << "info: using " << vertices->size() << " out of " << vertices->capacity()
       << " reserved vertex elements ("
       << Common::vectorStats(*vertices) << "%)\n";
}

hexTiling::hexTiling(uint N, const vertex& init, bool vis,
  bool sector) : vertices(0), onlyvis(vis), onlysector(sector) {
  const uint numsteps = 6;
  const vertex steps[6] = {
    vertex(-1, 2), vertex(-2, 1), vertex(-1, -1),
    vertex(1, -2), vertex(2, -1), vertex(1, 1)
  };

  // list of midpoints for the hexagons
  vlist midlist;
  midlist.push_back(init);

  level lvls[3];

  // level 0 is empty
  lvls[0].a = 0;
  lvls[0].b = 0;

  // level 1 only contains the init midpoint
  lvls[1].a = 0;
  lvls[1].b = 1;

  // level 2 is empty as well (new midpoints are inserted here)
  lvls[2].a = 1;
  lvls[2].b = 1;

  for (uint n = 0; n < N; ++n) {
    for (uint i = lvls[1].a; i != lvls[1].b; ++i) {
      // we need to copy this since the push_back call in the
      // inner loop might trigger a reallocation of the vector
      const vertex p(midlist[i]);

      for (uint j = 0; j < numsteps; ++j) {
        const vertex pp(p + steps[j]); // apply step

        // Search in all three levels...
        if (locate(midlist, lvls[0], pp) ||
            locate(midlist, lvls[1], pp) ||
            locate(midlist, lvls[2], pp)) continue;

        // ...before adding the point to the list:
        midlist.push_back(pp);
        ++lvls[2].b;
      }
    }

    // level exhausted, shift to next:

    // level 1 becomes level 0
    // level 2 becomes level 1
    lvls[0] = lvls[1];
    lvls[1] = lvls[2];

    // new level 2 is empty again
    lvls[2].a = lvls[1].b;
    lvls[2].b = lvls[1].b;
  }

  cerr << "info: using " << midlist.size() << " out of " << midlist.capacity()
       << " reserved midpoint elements ("
       << Common::vectorStats(midlist) << "%)\n";

  buildVertices(midlist);
}

#endif

