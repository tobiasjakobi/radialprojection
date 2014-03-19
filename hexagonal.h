#ifndef _HEXAGONAL_H_
#define _HEXAGONAL_H_

#include "radial_math.h"

#include <vector>

typedef unsigned int uint;
typedef unsigned short ushort;

using namespace std;

typedef vector<vec2i> vec2ilist;

class triTiling {
public:
  typedef vec2i vertex;
  typedef vector<vertex> vlist;
  
private:

  struct level {
    uint a, b;
  };

  bool locate(const level& l, const vertex& v);

  vlist* vertices;
  
public:
  triTiling(uint N, const vertex& init);

  ~triTiling() {
    delete vertices;
  }

  const vlist& getVertices() const {
    return *vertices;
  }
};

bool triTiling::locate(const level& l, const vertex& v) {
  assert(vertices != NULL);
  if (l.a == l.b) return false;

  uint i = l.a;

  while (i != l.b) {
    if ((*vertices)[i] == v) return true;

    ++i;
  }

  return false;
}

triTiling::triTiling(uint N, const vertex& init) {
  const uint numsteps = 6;
  const vertex steps[6] = {
    vertex(1, 0),  vertex(0, 1),  vertex(-1, 1),
    vertex(-1, 0), vertex(0, -1), vertex(1, -1)
  };

  // list of triangle vertices
  vertices = new vlist;
  vertices->push_back(init);

  level lvls[3];

  // level 0 is empty
  lvls[0].a = 0;
  lvls[0].b = 0;

  // level 1 only contains the init vertex
  lvls[1].a = 0;
  lvls[1].b = 1;

  // level 2 is empty as well (new vertices are inserted here)
  lvls[2].a = 1;
  lvls[2].b = 1;

  for (uint n = 0; n < N; ++n) {
    for (uint i = lvls[1].a; i != lvls[1].b; ++i) {
      const vertex p((*vertices)[i]);

      for (uint j = 0; j < numsteps; ++j) {
        const vertex pp(p + steps[j]);

        // Search in all three levels...
        if (locate(lvls[0], pp) ||
            locate(lvls[1], pp) ||
            locate(lvls[2], pp)) continue;

        // ...before adding the point to the list:
        vertices->push_back(pp);
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

  cerr << "info: using " << vertices->size() << " out of " << vertices->capacity()
       << " reserved vertex elements ("
       << CommonRadial::vectorStats(*vertices) << "%)\n";
}


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
        const vec2d x(t.transHexTo2D());

        if (!x.inFirstQuadrant()) continue;
        if (!x.inSectorL3()) continue;
      }

      // do visibility test if enabled
      if (onlyvis && RadialCoprime::gcdZ(abs(t.x), abs(t.y)) != 1)
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
       << CommonRadial::vectorStats(*vertices) << "%)\n";
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
       << CommonRadial::vectorStats(midlist) << "%)\n";

  buildVertices(midlist);
}

#endif

