/*  radialprojection - tools to numerically compute the radial projection of point sets
 *  Copyright (C) 2012-2016 - Tobias Jakobi <tjakobi at math dot uni dash bielefeld dot de>
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

#include "chiral_radial.h"

#include <sstream>
#include <algorithm>

/*
 * Remove this define to enable large data structures for the coordinates
 * used to create patches for the planar chair tiling.
 *
 * The regular structures use signed short integers, which work up to a
 * inflation factor of 14, but overflow for higher factors.
 * Removing this define essentially replaces the shorts with regular 32-bit
 * signed integers.
 *
 * Be aware that an inflation factor 15 already consumes more than 6GB of
 * system memory during construction.
 */
#define CHAIR_2D_LARGE_DISABLE

namespace ChiralLB {

  // Inflation factor:
  const double lambda = sqrt(0.5 * (5.0 + sqrt(5.0)));

  typedef vector<vec4s> vec4slist;

  ostream& operator<<(ostream &os, const rhomb& r) {
    os << '{';

    os << ((r.type == 0) ? "True" : "False");
    os << ',' << r.ref << ',' << r.rot;

    os << '}';

    return os;
  }

  const double VisOp::epsilon = 2.0 * numeric_limits<double>::epsilon();

  // Compute number of rhombs (of type A and B)
  // after another inflation step:
  void numRhombs(uint& a, uint& b);
  void numRhombs(uint& a, uint& b, uint steps);

  void countRhombs(const rhomblist& patch, uint& a, uint& b);
  uint countRhombs(const rhomblist& patch, uint steps);

  void iterate(const rhomblist& patch, uint steps, rhomblist& output);

  // The "sun" (though looking more like a star) is the main object
  // used for inflation (since it's a simple patch with maximum symmetry).
  void constructSun(rhomblist& sunPatch);

  // This is mainly a function to check if a 16-bit (short) integer
  // is enough to store the data produced during inflation.
  void minmax(const rhomblist& patch, vec8s& min, vec8s& max);

  void createVertices(Common::vec2dlist& vertices,
                      const rhomblist& initial, uint steps);
  void createVerticesVis(Common::vec2dlist& vertices, const rhomblist& initial,
                      uint steps, bool cutAndReduce);

  /* An optimized (and vastly faster) version of 'createVerticesVis', which also *
   * works satisfactory for "large" inflation multipliers. Currenly only the     *
   * "incorrect" visibility computation is implemented.                          */
  template <typename T>
  void createVerticesVisFast(Common::vec2dlist& vertices, const rhomblist& initial,
                      uint steps, bool cutAndReduce);

  // Compute length (eA and eB) of a rhomb of type A and B
  void getLength(double& typeA, double& typeB);

  // Cut a sector from the patch generated from inflating the sun:
  void cutSector(const Common::vec2dlist& input,
                 Common::vec2dlist& output, uint steps);

};

namespace Chair2D {

#ifdef CHAIR_2D_LARGE_DISABLE
  typedef vec2s coord_type;
#else
  typedef vec2i coord_type;
#endif

  typedef chairL<coord_type> item_type;
  typedef item_type::list_type list_type;

  template <typename T>
  ostream& operator<<(ostream &os, const chairL<T>& l) {
    os << '{'
       << l.getRef() << ',' << l.getRot()
       << '}';

    return os;
  }

  uint numL(uint l, uint steps);
  uint countL(const list_type& patch, uint steps);
  double clipradius(uint steps);

  void iterate(const list_type& patch, uint steps, list_type& output);

  /*
   * Variant of iterate() that applies clipping with the function 'cfnc'
   * in each step.
   */
  void iterateClip(const list_type& patch, uint steps,
                   clipfunc cfnc, list_type& output);

  void constructCross(list_type& crossPatch, bool cut);
  void minmax(const list_type& patch, coord_type& min, coord_type& max);

  void createVertices(Common::vec2dlist& vertices, uint steps);
  void createVerticesVis(Common::vec2dlist& vertices, uint steps, bool cut);

  // Cut a sector from the patch generated from inflating the cross:
  void cutSector(const Common::vec2dlist& input,
                 Common::vec2dlist& output, uint steps);

};

void ChiralLB::rhomb::inflate(rhomblist& list) const {
  const vec8s mid(midpoint());

  if (type == 0) {
    // type A decomposes into 3xA + 1xB
    list.push_back(rhomb(0, rot + 1, mid));
    list.push_back(rhomb(0, rot + 9, mid));
    list.push_back(rhomb(0, rot + 5, mid));
    list.push_back(rhomb(1, rot + 0, mid));
  } else {
    // type A decomposes into 1xA + 2xB
    list.push_back(rhomb(0, rot + 10, mid));
    list.push_back(rhomb(1, rot + 9,  mid));
    list.push_back(rhomb(1, rot + 1,  mid));
  }
}

void ChiralLB::rhomb::getVertices(vec8s* list) const {
  assert(list != NULL);

  list[0] = ref;

  if (type == 0) {
    list[1] = ref + (vec8s(18)).shift(rot);
    list[2] = ref + (vec8s(18) + vec8s(2)).shift(rot);
    list[3] = ref + (vec8s(2)).shift(rot);
  } else {
    list[1] = ref + (vec8s(11)).shift(rot);
    list[2] = ref + (vec8s(11) + vec8s(19)).shift(rot);
    list[3] = ref + (vec8s(19)).shift(rot);
  }
}

vec8s ChiralLB::rhomb::midpoint() const {
  const vec8s scaled(ref.lambdaScale());

  if (type == 0)
    return scaled + (vec8s(1) + vec8s(17)).shift(rot);
  else
    return scaled + (vec8s(18) + vec8s(10)).shift(rot);
}

template <typename T>
void Chair2D::chairL<T>::inflate(list_type& list) const {
  const coord_type newref(ref * 2);

  list.push_back(item_type(rot + 0, newref));
  list.push_back(item_type(rot + 0, newref + coord_type(1, 1).shift(rot)));
  list.push_back(item_type(rot + 1, newref + coord_type(4, 0).shift(rot)));
  list.push_back(item_type(rot + 3, newref + coord_type(0, 4).shift(rot)));
}

template <typename T>
void Chair2D::chairL<T>::getVertices(T* list) const {
  assert(list != NULL);

  list[0] = ref;
  list[1] = ref + coord_type(2, 0).shift(rot);
  list[2] = ref + coord_type(2, 1).shift(rot);
  list[3] = ref + coord_type(1, 1).shift(rot);
  list[4] = ref + coord_type(1, 2).shift(rot);
  list[5] = ref + coord_type(0, 2).shift(rot);
}

template <typename T>
bool Chair2D::chairL<T>::clip(double r) const {
  coord_type verts[6];
  const double rsq = r*r;

  this->getVertices(verts);
  for (uint i = 0; i < 6; ++i) {
    if (verts[i].transZ2ToR2().lengthSquared() <= rsq)
      return true;
  }

  return false;
}

void ChiralLB::numRhombs(uint& a, uint& b) {
  const uint x = 3*a + b;
  const uint y = a + 2*b;
  a = x; b = y;
}

void ChiralLB::numRhombs(uint& a, uint& b, uint steps) {
  for (uint i = 0; i < steps; ++i)
    numRhombs(a, b);
}

void ChiralLB::countRhombs(const rhomblist& patch, uint& a, uint& b) {
  for (rhomblist::const_iterator i = patch.begin(); i != patch.end(); ++i) {
    if (i->type == 0)
      ++a;
    else
      ++b;
  }
}

uint ChiralLB::countRhombs(const rhomblist& patch, uint steps) {
  uint a = 0, b = 0;
  countRhombs(patch, a, b);
  numRhombs(a, b, steps);
  return a + b;
}

void ChiralLB::iterate(const rhomblist& patch, uint steps, rhomblist& output) {
  // compute number of rhombs of type A and B of the initial patch
  uint a = 0, b = 0;
  countRhombs(patch, a, b);

  cerr << "Starting with an initial patch of " << a
       << " rhombs of type A and " << b
       << " rhombs of type B.\n";

  numRhombs(a, b, steps);

  rhomblist temp;
  temp.reserve(a + b);

  // Clear output, reserve enough space and add the initial patch
  output.clear();
  output.reserve(a + b);
  output.insert(output.end(), patch.begin(), patch.end());

  for (uint i = 0; i < steps; ++i) {
    temp.clear();

    for (rhomblist::const_iterator j = output.begin(); j != output.end(); ++j)
      j->inflate(temp);

    output.swap(temp);
  }

  cerr << "After " << steps << " inflation steps the resulting patch has "
       << a << " rhombs of type A and " << b
       << " rhombs of type B.\n";
  cerr << "Amount of space used by data structures is "
       << sizeof(rhomb) * output.size() << " bytes.\n";
}

void ChiralLB::constructSun(rhomblist& sunPatch) {
  sunPatch.clear();
  sunPatch.reserve(5);

  sunPatch.push_back(rhomb(0, 15, vec8s(7)  + vec8s(3)));
  sunPatch.push_back(rhomb(0, 19, vec8s(7)  + vec8s(11)));
  sunPatch.push_back(rhomb(0, 3,  vec8s(11) + vec8s(15)));
  sunPatch.push_back(rhomb(0, 7,  vec8s(15) + vec8s(19)));
  sunPatch.push_back(rhomb(0, 11, vec8s(19) + vec8s(3)));
}

void ChiralLB::minmax(const rhomblist& patch, vec8s& min, vec8s& max) {
  min.set(0, 0, 0, 0, 0, 0, 0, 0);
  max.set(0, 0, 0, 0, 0, 0, 0, 0);

  for (rhomblist::const_iterator i = patch.begin(); i != patch.end(); ++i) {
    const vec8s& ref = i->ref;
    for (uint j = 0; j < 8; ++j) {
      if (ref[j] > max[j])
        max[j] = ref[j];

      if (ref[j] < min[j])
        min[j] = ref[j];
    }
  }
}

void ChiralLB::createVertices(Common::vec2dlist& vertices,
                 const rhomblist& initial, uint steps) {
  if (initial.empty())
    return;

  rhomblist* patch = new rhomblist;
  iterate(initial, steps, *patch);

  {
    vec8s min, max;
    minmax(*patch, min, max);
    cerr << "statistics: min = " << min << ", max = " << max << endl;
  }

  // the factor 1.3 derives from empirical evaluation
  vec4slist verts;
  verts.reserve(double(countRhombs(initial, steps)) * 1.3);

  bool reducemode;

  // Check if all vertices can be reduced into L10:
  {
    uint even = 0, odd = 0;
    for (rhomblist::const_iterator i = patch->begin(); i != patch->end(); ++i) {
      vec8s temp[4];
      i->getVertices(temp);

      for (uint j = 0; j < 4; ++j) {
        if (!temp[j].isInL10(even, odd)) {
          cerr << "error: vector found that isn't in L10" << endl;
          return;
        }
      }
    }

    if (even == patch->size() * 4) {
      reducemode = true;
    } else {
      if (odd == patch->size() * 4) {
        reducemode = false;
      } else {
        cerr << "error: vectors can't be reduced into L10" << endl;
        return;
      }
    }
  }

  for (rhomblist::const_iterator i = patch->begin(); i != patch->end(); ++i) {
    vec8s temp[4];
    i->getVertices(temp);

    for (uint j = 0; j < 4; ++j)
      verts.push_back(temp[j].reduceToL10(reducemode));
  }

  // Remove duplicate vertices
  sort(verts.begin(), verts.end());
  verts.erase(unique(verts.begin(), verts.end()), verts.end());

  cerr << "statistics: " << patch->size() << " rhombs reduced to "
       << verts.size() << " unique vertices" << endl;

  delete patch;
  patch = NULL;

  vertices.clear();
  vertices.reserve(verts.size());
  for (vec4slist::const_iterator i = verts.begin(); i != verts.end(); ++i) {
    vertices.push_back(i->transL10ToR2());
  }

  cerr << "info: used " << verts.size() << " out of " << verts.capacity()
       << " allocated (vector) elements ("
       <<  (double(verts.size()) / double(verts.capacity()) * 100.0)
       << "%).\n";
}

void ChiralLB::createVerticesVis(Common::vec2dlist& vertices,
                        const rhomblist& initial,
                        uint steps, bool cutAndReduce) {
  if (initial.empty())
    return;

  rhomblist* patch = new rhomblist;
  iterate(initial, steps, *patch);

  // the factors 0.15 and 1.4 were derived from empirical evaluations
  VisList* vlist = new VisList;
  vlist->reserve(double(countRhombs(initial, steps)) *
                 (cutAndReduce ? 0.15 : 1.4));

  vlist->init();

  bool reducemode;

  // Check if all vertices can be reduced into L10:
  {
    uint even = 0, odd = 0;
    for (rhomblist::const_iterator i = patch->begin(); i != patch->end(); ++i) {
      vec8s temp[4];
      i->getVertices(temp);

      for (uint j = 0; j < 4; ++j) {
        if (!temp[j].isInL10(even, odd)) {
          cerr << "error: vector found that isn't in L10" << endl;
          return;
        }
      }
    }

    if (even == patch->size() * 4) {
      reducemode = true;
    } else {
      if (odd == patch->size() * 4) {
        reducemode = false;
      } else {
        cerr << "error: vectors can't be reduced into L10" << endl;
        return;
      }
    }
  }
  
  if (cutAndReduce) {
    cerr << "info: trimming the tiling into a circular area\n";
    const double cutoff = Common::power(lambda, steps);
    for (rhomblist::const_iterator i = patch->begin(); i != patch->end(); ++i) {
      vec8s temp[4];
      i->getVertices(temp);

      for (uint j = 0; j < 4; ++j) {
        const vec4s reduced(temp[j].reduceToL10(reducemode));

        // extract a circular region from the patch
        const vec2d phys(reduced.transL10ToR2());
        if (phys.lengthSquared() > cutoff*cutoff) continue;

        // further reduce the region to remove redundant information (symmetry)
        if (!phys.inFirstQuadrant()) continue;
        if (!phys.inSectorL5()) continue;

        vlist->insertSorted(reduced);
      }
    }
  } else {
    for (rhomblist::const_iterator i = patch->begin(); i != patch->end(); ++i) {
      vec8s temp[4];
      i->getVertices(temp);

      for (uint j = 0; j < 4; ++j) {
        const vec4s reduced(temp[j].reduceToL10(reducemode));

        if (reduced.isZero()) continue;

        vlist->insertSorted(reduced);
      }
    }
  }

  cerr << "statistics: " << patch->size() << " rhombs reduced to "
       << vlist->size() << " unique vertices." << endl;

  delete patch;
  patch = NULL;

  vlist->removeInvisibleFast();

  vertices.clear();
  vertices.reserve(vlist->size());
  vlist->toR2(vertices);

  delete vlist;
  vlist = NULL;
}

template <typename T>
void ChiralLB::createVerticesVisFast(Common::vec2dlist& vertices,
                        const rhomblist& initial,
                        uint steps, bool cutAndReduce) {
  if (initial.empty())
    return;

  rhomblist* patch = new rhomblist;
  iterate(initial, steps, *patch);

  bool reducemode;

  // Check if all vertices can be reduced into L10:
  {
    uint even = 0, odd = 0;
    for (rhomblist::const_iterator i = patch->begin(); i != patch->end(); ++i) {
      vec8s temp[4];
      i->getVertices(temp);

      for (uint j = 0; j < 4; ++j) {
        if (!temp[j].isInL10(even, odd)) {
          cerr << "error: vector found that isn't in L10" << endl;
          return;
        }
      }
    }

    if (even == patch->size() * 4) {
      reducemode = true;
    } else {
      if (odd == patch->size() * 4) {
        reducemode = false;
      } else {
        cerr << "error: vectors can't be reduced into L10" << endl;
        return;
      }
    }
  }

  vector<T> vlist;
  vlist.reserve(patch->size() * 4);

  if (cutAndReduce) {
    cerr << "info: trimming the tiling into a circular area\n";
    const double cutoff = Common::power(lambda, steps);
    for (rhomblist::const_iterator i = patch->begin(); i != patch->end(); ++i) {
      vec8s temp[4];
      i->getVertices(temp);

      for (uint j = 0; j < 4; ++j) {
        const vec4s reduced(temp[j].reduceToL10(reducemode));

        // extract a circular region from the patch
        const vec2d phys(reduced.transL10ToR2());
        if (phys.lengthSquared() > cutoff*cutoff) continue;

        // further reduce the region to remove redundant information (symmetry)
        if (!phys.inFirstQuadrant()) continue;
        if (!phys.inSectorL5()) continue;

        vlist.push_back(reduced);
      }
    }
  } else {
    for (rhomblist::const_iterator i = patch->begin(); i != patch->end(); ++i) {
      vec8s temp[4];
      i->getVertices(temp);

      for (uint j = 0; j < 4; ++j) {
        const vec4s reduced(temp[j].reduceToL10(reducemode));

        if (reduced.isZero()) continue;

        vlist.push_back(reduced);
      }
    }
  }

  // First removal pass
  sort(vlist.begin(), vlist.end());
  vlist.erase(unique(vlist.begin(), vlist.end()), vlist.end());

  cerr << "statistics: " << patch->size() << " rhombs reduced to "
       << vlist.size() << " unique vertices." << endl;

  delete patch;
  patch = NULL;

  // TODO: implement correct visibility computation

  /*
   * This vector looks like it can be processed faster by using threads, but
   * performance evaluation has shown that the advantage is negligible at best.
   */
  for (typename vector<T>::iterator i = vlist.begin(); i != vlist.end(); ++i)
    *i = i->transL10ToDirect().directL10ToUnique();

  // Second removal pass
  sort(vlist.begin(), vlist.end());
  vlist.erase(unique(vlist.begin(), vlist.end()), vlist.end());

  cerr << "statistics: after (incorrect) visibility computation: "
       << vlist.size() << " vertices visible." << endl;

  vertices.clear();
  vertices.reserve(vlist.size());

  for (typename vector<T>::const_iterator i = vlist.begin(); i != vlist.end(); ++i)
    vertices.push_back(i->directL10ToR2());
}

void ChiralLB::getLength(double& typeA, double& typeB) {
  rhomb a(0, 0, vec8s());
  rhomb b(1, 0, vec8s());

  vec8s vertsA[4], vertsB[4];
  a.getVertices(vertsA);
  b.getVertices(vertsB);

  typeA = sqrt((vertsA[0].transL20ToR2() - vertsA[2].transL20ToR2()).lengthSquared());
  typeB = sqrt((vertsB[1].transL20ToR2() - vertsB[3].transL20ToR2()).lengthSquared());
}

void ChiralLB::cutSector(const Common::vec2dlist& input,
                   Common::vec2dlist& output, uint steps) {
  using namespace Common;

  const double cutoff = power(lambda, steps);

  output.clear();
  output.reserve(double(input.size()) * 0.75);

  for (vec2dlist::const_iterator i = input.begin(); i != input.end(); ++i) {
    if (i->lengthSquared() <= cutoff*cutoff)
      output.push_back(*i);
  }

  cerr << "After cutting off procedure " << output.size() << " vertices remain.\n";
}

uint Chair2D::numL(uint l, uint steps) {
  for (uint i = 0; i < steps; ++i) {
    l *= 4;
  }

  return l;
}

uint Chair2D::countL(const list_type& patch, uint steps) {
  const uint size = patch.size();

  return numL(size, steps);
}

double Chair2D::clipradius(uint steps) {
  return sqrt(2.0) * Common::power(2.0, steps);
}

void Chair2D::iterate(const list_type& patch, uint steps, list_type& output) {
  const uint l = countL(patch, steps);

  cerr << "Starting with an initial patch of " << patch.size()
       << " L-shaped tiles.\n";

  list_type temp;
  temp.reserve(l);

  // Clear output, reserve enough space and add the initial patch
  output.clear();
  output.reserve(l);
  output.insert(output.end(), patch.begin(), patch.end());

  for (uint i = 0; i < steps; ++i) {
    temp.clear();

    for (list_type::const_iterator j = output.begin(); j != output.end(); ++j)
      j->inflate(temp);

    output.swap(temp);
  }

  cerr << "After " << steps << " inflation steps the resulting patch has "
       << l << " L-shaped tiles.\n";
  cerr << "Amount of space used by data structures is "
       << sizeof(item_type) * output.size() << " bytes.\n";
}

void Chair2D::iterateClip(const list_type& patch, uint steps,
                  clipfunc cfnc, list_type& output) {
  const uint l = double(countL(patch, steps)) * 0.6;

  cerr << "Starting with an initial patch of " << patch.size()
       << " L-shaped tiles." << endl;
  cerr << "Clipping to circular shape is done in every iteration." << endl;
  cerr << "Allocated storage for " << l << " L-shaped tiles." << endl;

  list_type temp, to_clip;
  temp.reserve(l);

  // Clear output, reserve enough space and add the initial patch
  output.clear();
  output.reserve(l);
  output.insert(output.end(), patch.begin(), patch.end());

  for (uint i = 0; i < steps; ++i) {
    const double cr = cfnc(i + 1);
    temp.clear();

    for (list_type::const_iterator j = output.begin(); j != output.end(); ++j) {
      to_clip.clear();
      j->inflate(to_clip);

      for (list_type::const_iterator k = to_clip.begin(); k != to_clip.end(); ++k) {
        if (k->clip(cr))
          temp.push_back(*k);
      }
    }

    output.swap(temp);
  }

  cerr << "After " << steps << " inflation steps the resulting patch has "
       << output.size() << " L-shaped tiles ("
       << (100.0 * double(output.size()) / double(l))
       << " percent of allocated storage used)." << endl;
  cerr << "Amount of space used by data structures is "
       << sizeof(item_type) * output.size() << " bytes.\n";
}

void Chair2D::constructCross(list_type& crossPatch, bool cut) {
  typedef item_type::coord_type coord_type;

  crossPatch.clear();
  crossPatch.reserve(cut ? 1 : 4);

  crossPatch.push_back(item_type(0, coord_type(0, 0)));

  if (!cut) {
    crossPatch.push_back(item_type(1, coord_type(0, 0)));
    crossPatch.push_back(item_type(2, coord_type(0, 0)));
    crossPatch.push_back(item_type(3, coord_type(0, 0)));
  }
}

void Chair2D::minmax(const list_type& patch, coord_type& min, coord_type& max) {
  min.set(0, 0);
  max.set(0, 0);

  for (list_type::const_iterator i = patch.begin(); i != patch.end(); ++i) {
    const coord_type& ref = i->getRef();
    for (uint j = 0; j < 2; ++j) {
      if (ref[j] > max[j]) {
        max[j] = ref[j];
      }

      if (ref[j] < min[j]) {
        min[j] = ref[j];
      }
    }
  }
}

void Chair2D::createVertices(Common::vec2dlist& vertices, uint steps) {
  typedef vector<coord_type> vertexlist_type;
  typedef vertexlist_type::const_iterator iter_type;

  list_type seed;
  constructCross(seed, false);

  list_type* patch = new list_type;
  iterate(seed, steps, *patch);

  {
    coord_type min, max;
    minmax(*patch, min, max);
    cerr << "statistics: min = " << min << ", max = " << max << endl;
  }

  vertexlist_type verts;
  verts.reserve(double(countL(seed, steps)) * 3.0);

  for (list_type::const_iterator i = patch->begin(); i != patch->end(); ++i) {
    coord_type temp[6];
    i->getVertices(temp);

    for (uint j = 0; j < 6; ++j)
      verts.push_back(temp[j]);
  }

  // Remove duplicate vertices
  sort(verts.begin(), verts.end());
  verts.erase(unique(verts.begin(), verts.end()), verts.end());

  cerr << "statistics: " << patch->size() << " L-shaped tiles reduced to "
       << verts.size() << " unique vertices" << endl;

  delete patch;
  patch = NULL;

  vertices.clear();
  vertices.reserve(verts.size());
  for (iter_type i = verts.begin(); i != verts.end(); ++i)
    vertices.push_back(i->transZ2ToR2());
}

void Chair2D::createVerticesVis(Common::vec2dlist& vertices,
        uint steps, bool cut) {
  typedef vector<coord_type> vertexlist_type;
  typedef vertexlist_type::const_iterator iter_type;

  list_type seed;
  constructCross(seed, cut);

  list_type* patch = new list_type;

  if (cut)
    iterateClip(seed, steps, clipradius, *patch);
  else
    iterate(seed, steps, *patch);

  {
    coord_type min, max;
    minmax(*patch, min, max);
    cerr << "statistics: min = " << min << ", max = " << max << endl;
  }

  vertexlist_type verts;
  verts.reserve(double(countL(seed, steps)) * 3.0);

  // create a "occupation" map of the tiling vertices
  VisibilityMap<item_type>* occupied = new VisibilityMap<item_type>(*patch, steps);

  if (cut) {
    cerr << "info: trimming the tiling into a circular area\n";
    const double cutoff = sqrt(2.0) * Common::power(2.0, steps);
    for (list_type::const_iterator i = patch->begin(); i != patch->end(); ++i) {
      coord_type temp[6];
      i->getVertices(temp);

      for (uint j = 0; j < 6; ++j) {
        const coord_type current(temp[j]);

        // check for visibility
        if (!occupied->isVisible(current)) continue;

        // extract a circular region from the patch
        const vec2d phys(current.transZ2ToR2());
        if (phys.lengthSquared() > cutoff*cutoff) continue;

        // further reduce the region to remove redundant information (symmetry)
        if (phys.inFirstQuadrant())
          verts.push_back(current);
      }
    }
  } else {
    for (list_type::const_iterator i = patch->begin(); i != patch->end(); ++i) {
      coord_type temp[6];
      i->getVertices(temp);

      for (uint j = 0; j < 6; ++j) {
        const coord_type current(temp[j]);

        if (occupied->isVisible(current))
          verts.push_back(current);
      }
    }
  }

  // Remove duplicate vertices
  sort(verts.begin(), verts.end());
  verts.erase(unique(verts.begin(), verts.end()), verts.end());

  delete occupied;
  occupied = NULL;

  cerr << "statistics: " << patch->size() << " L-shaped tiles reduced to "
       << verts.size() << " unique vertices" << endl;

  delete patch;
  patch = NULL;

  vertices.clear();
  vertices.reserve(verts.size());
  for (iter_type i = verts.begin(); i != verts.end(); ++i)
    vertices.push_back(i->transZ2ToR2());
}

void Chair2D::cutSector(const Common::vec2dlist& input,
               Common::vec2dlist& output, uint steps) {
  using namespace Common;

  const double cutoff = sqrt(2.0) * power(2.0, steps);

  output.clear();
  output.reserve(double(input.size()) * 0.6);

  for (vec2dlist::const_iterator i = input.begin(); i != input.end(); ++i) {
    if (i->lengthSquared() <= cutoff*cutoff)
      output.push_back(*i);
  }

  cerr << "After cutting off procedure " << output.size() << " vertices remain.\n";
}

void print_usage() {
  cerr << "chiral_radial: usage:" << endl;

  cerr << "chiral_radial --chiral: selects chiral LB main mode" << endl;
  cerr << "chiral_radial --chair: selects 2D chair main mode" << endl;
  cerr << "(both modes share the same set of parameters)" << endl;

  cerr << "\tparameter 1: mode" << endl;
  cerr << "\t\t" "0 = tiles in Mathematica format" << endl;
  cerr << "\t\t" "1 = vertices in Mathematica format" << endl;
  cerr << "\t\t" "2 = spacings from radial projection" << endl;

  cerr << "\tparameter 2: steps (increase size of the patch)" << endl;
  cerr << "\tparameter 3: cut patch to circular shape" << endl;

  cerr << endl;

  cerr << "Passing --second-order as second argument switches from first to second"
       << endl << "order spacings (this only affects the radial projection mode)."
       << endl;
  cerr << "Passing --visible-vertex as second argument computes only visible"
       << endl << "vertices for mode 1." << endl;
}

int main_chiral(int argc, char* argv[]) {
  stringstream parser;

  uint steps = 2;
  uint mode = 0;
  bool cut = false;
  bool second_order = false;
  bool visible_vertex = false;

  using namespace ChiralLB;
  using namespace Common;

  while (argc >= 2) {
    const string arg(argv[1]);
    bool parsefail = false;

    if (arg == "--second-order")
      second_order = true;
    else if (arg == "--visible-vertex")
      visible_vertex = true;
    else
      parsefail = true;

    if (parsefail) break;

    argc--;
    argv++;
  }

  if (argc >= 2) {
    parser.str(argv[1]);
    parser.clear();
    parser >> mode;

    if (argc >= 3) {
      parser.str(argv[2]);
      parser.clear();
      parser >> steps;

      if (argc >= 4) {
        parser.str(argv[3]);
        parser.clear();
        parser >> cut;
      }
    }
  }

  cerr << "info: size of a single rhomb is " << sizeof(rhomb) << " bytes.\n";

  // We always use the "sun" (patch consisting of five rhombs
  // of type A) as the initial patch:
  rhomblist initialChiral;
  constructSun(initialChiral);

  rhomblist rhombs;
  vec2dlist verts;
  dlist spacings;
  double mean;
  vec8s min, max;

  switch (mode) {
  // Output rhomb data in Mathematica style (cut is ignored)
  case 0:
    iterate(initialChiral, steps, rhombs);
    minmax(rhombs, min, max);

    cerr << "statistics: min = " << min << ", max = " << max << endl;
    cout << rhombs;
    break;

  // Output vertex data in Mathematica style
  case 1:
    if (visible_vertex)
      createVerticesVis(verts, initialChiral, steps, cut);
    else
      createVertices(verts, initialChiral, steps);

    cout << verts;
    break;

  // Do radial projection and output data in raw mode
  case 2:
    // Switch to higher internal accuracy when using a least 13 inflation steps.
    if (steps > 12)
      createVerticesVisFast<vec4i>(verts, initialChiral, steps, cut);
    else
      createVerticesVisFast<vec4s>(verts, initialChiral, steps, cut);

    radialProj(verts, spacings, mean);

    meanDistanceMessage(verts.size(), mean);

    if (second_order) {
      cerr << "info: computing second-order spacings." << endl;

      vec2dlist spacings2;
      spacings2.reserve(spacings.size() - 1);

      secondOrderSpacings(spacings, spacings2);
      writeRawConsole(spacings2);
    } else {
      writeRawConsole(spacings);
    }
    break;

  default:
    cerr << "error: unsupported mode selected!\n";
    return 1;
  }

  return 0;
}

int main_chair(int argc, char* argv[]) {
  stringstream parser;

  uint steps = 2;
  uint mode = 0;
  bool cut = false;
  bool second_order = false;
  bool visible_vertex = false;

  using namespace Chair2D;
  using namespace Common;

  while (argc >= 2) {
    const string arg(argv[1]);
    bool parsefail = false;

    if (arg == "--second-order")
      second_order = true;
    else if (arg == "--visible-vertex")
      visible_vertex = true;
    else
      parsefail = true;

    if (parsefail) break;

    argc--;
    argv++;
  }

  if (argc >= 2) {
    parser.str(argv[1]);
    parser.clear();
    parser >> mode;

    if (argc >= 3) {
      parser.str(argv[2]);
      parser.clear();
      parser >> steps;

      if (argc >= 4) {
        parser.str(argv[3]);
        parser.clear();
        parser >> cut;
      }
    }
  }

#ifdef CHAIR_2D_LARGE_DISABLE
  /*
   * We use signed shorts in the data structure for the tiles of the
   * planar chair tiling to keep memory usage low.
   * With 14 inflation steps the coordinates reach a maximum value of
   * 23172. The signed short type can store values in the range
   * [-32768, 32767]. Applying another inflation step doubles the values,
   * leading to integer under- and overflow.
   */
  if (steps >= 15) {
    cerr << "error: maximum number of steps is 14 (at the moment)" << endl;
    return 1;
  }
#endif

  list_type tiles;
  vec2dlist verts;
  dlist spacings;
  double mean;
  coord_type min, max;

  switch (mode) {
  // Output rhomb data in Mathematica style (cut is ignored)
  case 0:
    {
      list_type seed;
      constructCross(seed, true);
      iterate(seed, steps, tiles);
      minmax(tiles, min, max);
    }

    cerr << "statistics: min = " << min << ", max = " << max << endl;
    cout << tiles;
    break;

  // Output vertex data in Mathematica style
  case 1:
    if (visible_vertex)
      createVerticesVis(verts, steps, cut);
    else
      createVertices(verts, steps);

    cout << verts;
    break;

  // Do radial projection and output data in raw mode
  case 2:
    createVerticesVis(verts, steps, cut);
    radialProj(verts, spacings, mean);

    meanDistanceMessage(verts.size(), mean);

    if (second_order) {
      cerr << "info: computing second-order spacings." << endl;

      vec2dlist spacings2;
      spacings2.reserve(spacings.size() - 1);

      secondOrderSpacings(spacings, spacings2);
      writeRawConsole(spacings2);
    } else {
      writeRawConsole(spacings);
    }
    break;

  default:
    cerr << "error: unsupported mode selected!\n";
    return 2;
  }

  return 0;
}

int main(int argc, char* argv[]) {
  stringstream parser;
  string tempstr;

  int main_mode = -1;

  if (argc >= 2) {
    parser.str(argv[1]);
    parser.clear();
    parser >> tempstr;
  }

  if (tempstr == "--chiral")
    main_mode = 0;
  else if (tempstr == "--chair")
    main_mode = 1;

  if (main_mode == -1) {
    print_usage();
    return 0;
  }

  switch (main_mode) {
  case 0:
    return main_chiral(argc - 1, argv + 1);

  case 1:
    return main_chair(argc - 1, argv + 1);

  default:
    assert(false);
    return 0;
  }
}
