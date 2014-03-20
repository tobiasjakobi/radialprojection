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

#include "tuebingen.h"

#include <sstream>

namespace TuebingenTri {

  ostream& operator<<(ostream &os, const tri& t) {
    os << '{';

    os << ((t.getType() == 1) ? "True" : "False") << ','
       << ((t.getOrient() == 1) ? "True" : "False");
    
    os << ',' << t.getRef() << ',' << t.getRot();

    os << '}';

    return os;
  }

  const double VisOp::epsilon = 2.0 * numeric_limits<double>::epsilon();

  // Compute number of triangles (of type A and B)
  // after another inflation step:
  void numTris(uint& a, uint& b);
  void numTris(uint& a, uint& b, uint steps);

  void countTris(const trilist& patch, uint& a, uint& b);
  uint countTris(const trilist& patch, uint steps);

  void iterate(const trilist& patch, uint steps, bool ttt, trilist& output);

  // The "sun" (consisting of 10 triangles of type A) is the main object
  // used for inflation (since it's a simple patch with maximum symmetry).
  void constructSun(trilist& sunPatch, bool ttt);

  // This is mainly a function to check if a 16-bit (short) integer
  // is enough to store the data produced during inflation.
  void minmax(const trilist& patch, vec4s& min, vec4s& max);

  /* Create set of visible vertices from an initial patch. Vertices are
   * output as double 2-tuples (cartesian format) and can therefore be
   * directly used for radial projection.
   * Also the visibility test uses some optimization for removing
   * invisible vertices, which works fine when using the output
   * for radial projection. Don't use this when you want the
   * correct set of visible vertices!
   *
   * steps = number of inflations steps applied to the initial patch
   * sector = only consider a 1/10-sector of the initial patch
   * ttt = apply TTT inflation rule (vs. PRT inflation rule)
   */
  void createVerticesVis(Common::vec2dlist& vertices, const trilist& initial,
                         uint steps, bool sector, bool ttt);

  /* Same as createVerticesVis but outputs the vertices as integer
   * 4-tuples (exact representation) and also correctly determines
   * invisible vertices. This of course slows down the computations,
   * so use with care!
   */
  void createVerticesVis2(Common::vec4slist& vertices, const trilist& initial,
                          uint steps, bool sector, bool ttt);

};

void TuebingenTri::tri::inflateTTT(trilist& list) const {
  if (type == 1) {
    // type A decomposes into 2xA + 1xB
    inflateTTTa(list);
  } else {
    // type B decomposes into 1xA + 1xB
    inflateTTTb(list);
  }
}

void TuebingenTri::tri::inflatePRT(trilist& list) const {
  if (type == 1) {
    // type A decomposes into 2xA + 1xB
    inflatePRTa(list);
  } else {
    // type B decomposes into 1xA + 1xB
    inflatePRTb(list);
  }
}

void TuebingenTri::tri::getVertices(vec4s* list) const {
  assert(list != NULL);

  if (type == 1)
    getVerticesA(list);
  else
    getVerticesB(list);
}

void TuebingenTri::tri::inflateTTTa(trilist& list) const {
  static const vec4s directions[6] = {
    vec4s(0, 0, 0, 0),
    vec4s(0, 0, -1, -1),
    vec4s(1, 1, 0, 0),
    vec4s(0, 0, 0, 0),
    vec4s(0, 0, 1, 1),
    vec4s(0, 1, 1, 1)
  };

  const vec4s newref(ref.scaleTauL5());

  if (orient == 1) {
    list.push_back(tri(0, 1, rot + 1, newref + directions[0].shiftL5ByL10(rot)));
    list.push_back(tri(1, 1, rot + 3, newref + directions[1].shiftL5ByL10(rot)));
    list.push_back(tri(1, 0, rot + 0, newref + directions[2].shiftL5ByL10(rot)));
  } else {
    list.push_back(tri(0, 0, rot + 9, newref + directions[3].shiftL5ByL10(rot)));
    list.push_back(tri(1, 0, rot + 7, newref + directions[4].shiftL5ByL10(rot)));
    list.push_back(tri(1, 1, rot + 0, newref + directions[5].shiftL5ByL10(rot)));
  }
}

void TuebingenTri::tri::inflateTTTb(trilist& list) const {
  static const vec4s directions[4] = {
    vec4s(1, 0, -1, -1),
    vec4s(1, 1, 0, 0),
    vec4s(-1, 0, 1, 1),
    vec4s(0, 1, 1, 1)
  };

  const vec4s newref(ref.scaleTauL5());

  if (orient == 1) {
    list.push_back(tri(0, 1, rot + 4, newref + directions[0].shiftL5ByL10(rot)));
    list.push_back(tri(1, 0, rot + 3, newref + directions[1].shiftL5ByL10(rot)));
  } else {
    list.push_back(tri(0, 0, rot + 6, newref + directions[2].shiftL5ByL10(rot)));
    list.push_back(tri(1, 1, rot + 7, newref + directions[3].shiftL5ByL10(rot)));
  }
}

void TuebingenTri::tri::inflatePRTa(trilist& list) const {
  static const vec4s directions[6] = {
    vec4s(0, 0, -1, -1),         
    vec4s(1, 1, 1, 0),
    vec4s(1, 2, 1, 0),
    vec4s(0, 0, 1, 1),
    vec4s(0, 1, 1, 0),
    vec4s(1, 2, 2, 1)
  };

  const vec4s newref(ref.scaleTauL5());

  if (orient == 1) {
    list.push_back(tri(1, 1, rot + 3, newref + directions[0].shiftL5ByL10(rot)));
    list.push_back(tri(1, 0, rot + 4, newref + directions[1].shiftL5ByL10(rot)));
    list.push_back(tri(0, 0, rot + 3, newref + directions[2].shiftL5ByL10(rot)));
  } else {
    list.push_back(tri(1, 0, rot + 7, newref + directions[3].shiftL5ByL10(rot)));
    list.push_back(tri(1, 1, rot + 6, newref + directions[4].shiftL5ByL10(rot)));
    list.push_back(tri(0, 1, rot + 7, newref + directions[5].shiftL5ByL10(rot)));
  }
};

void TuebingenTri::tri::inflatePRTb(trilist& list) const {
  static const vec4s directions[4] = {
    vec4s(1, 0, -1, -1),
    vec4s(0, 0, -1, -1),
    vec4s(-1, 0, 1, 1),
    vec4s(0, 0, 1, 1)
  };

  const vec4s newref(ref.scaleTauL5());

  if (orient == 1) {
    list.push_back(tri(0, 1, rot + 4, newref + directions[0].shiftL5ByL10(rot)));
    list.push_back(tri(1, 1, rot + 3, newref + directions[1].shiftL5ByL10(rot)));
  } else {
    list.push_back(tri(0, 0, rot + 6, newref + directions[2].shiftL5ByL10(rot)));
    list.push_back(tri(1, 0, rot + 7, newref + directions[3].shiftL5ByL10(rot)));
  }
};

void TuebingenTri::tri::getVerticesA(vec4s* list) const {
  static const vec4s directions[4] = {
    vec4s(1, 0, 0, 0),
    vec4s(1, 1, 1, 0),
    vec4s(-1, 0, 0, 0),
    vec4s(0, 1, 1, 0)
  };

  list[0] = ref;

  if (orient == 1) {
    list[1] = ref + directions[0].shiftL5ByL10(rot);
    list[2] = ref + directions[1].shiftL5ByL10(rot);
  } else {
    list[1] = ref + directions[2].shiftL5ByL10(rot);
    list[2] = ref + directions[3].shiftL5ByL10(rot);
  }
}

void TuebingenTri::tri::getVerticesB(vec4s* list) const {
  static const vec4s directions[4] = {
    vec4s(0, 0, -1, -1),
    vec4s(0, 0, 0, -1),
    vec4s(0, 0, 1, 1),
    vec4s(0, 0, 1, 0)
  };

  list[0] = ref;

  if (orient == 1) {
    list[1] = ref + directions[0].shiftL5ByL10(rot);
    list[2] = ref + directions[1].shiftL5ByL10(rot);
  } else {
    list[1] = ref + directions[2].shiftL5ByL10(rot);
    list[2] = ref + directions[3].shiftL5ByL10(rot);
  }
}

void TuebingenTri::numTris(uint& a, uint& b) {
  const uint x = 2*a + b;
  const uint y = a + b;
  a = x; b = y;
}

void TuebingenTri::numTris(uint& a, uint& b, uint steps) {
  for (uint i = 0; i < steps; ++i) {
    numTris(a, b);
  }
}

void TuebingenTri::countTris(const trilist& patch, uint& a, uint& b) {
  for (trilist::const_iterator i = patch.begin(); i != patch.end(); ++i) {
    if (i->getType() == 1)
      ++a;
    else
      ++b;
  }
}

uint TuebingenTri::countTris(const trilist& patch, uint steps) {
  uint a = 0, b = 0;
  countTris(patch, a, b);
  numTris(a, b, steps);
  return a + b;
}

void TuebingenTri::iterate(const trilist& patch, uint steps, bool ttt,
                           trilist& output) {

  // compute number of triangles of type A and B of the initial patch
  uint a = 0, b = 0;
  countTris(patch, a, b);

  cerr << "Starting with an initial patch of " << a
       << " triangles of type A and " << b
       << " triangles of type B.\n";

  numTris(a, b, steps);

  trilist temp;
  temp.reserve(a + b);

  // Clear output, reserve enough space and add the initial patch
  output.clear();
  output.reserve(a + b);
  output.insert(output.end(), patch.begin(), patch.end());

  if (ttt) {
    cerr << "Applying TTT / Tuebingen triangle tiling inflation rules!\n";
    for (uint i = 0; i < steps; ++i) {
      temp.clear();

      for (trilist::const_iterator j = output.begin(); j != output.end(); ++j)
        j->inflateTTT(temp);

      output.swap(temp);
    }
  } else {
    cerr << "Applying PRT / Penrose-Robinson tiling inflation rules!\n";
    for (uint i = 0; i < steps; ++i) {
      temp.clear();

      for (trilist::const_iterator j = output.begin(); j != output.end(); ++j)
        j->inflatePRT(temp);

      output.swap(temp);
    }
  }

  cerr << "After " << steps << " inflation steps the resulting patch has "
       << a << " triangles of type A and " << b
       << " triangles of type B.\n";
  cerr << "Amount of space used by data structures is "
       << sizeof(tri) * output.size() << " bytes.\n";
}

void TuebingenTri::constructSun(trilist& sunPatch, bool ttt) {
  static const vec4s directionsTTT[7] = {
    vec4s(-1, -1, -1, 0),
    vec4s(0, -1, -1, -1),
    vec4s(0, 0, -1, -1),
    vec4s(1, 1, 0, 0),
    vec4s(0, 1, 1, 0),
    vec4s(0, 1, 1, 1),
    vec4s(-1, -1, 0, 0)
  };

  static const vec4s directionsPRT[5] = {
    vec4s(-1, -1, -1, 0),
    vec4s(0, -1, -1, -1),
    vec4s(1, 1, 0, 0),
    vec4s(0, 1, 1, 0),
    vec4s(0, 0, 1, 1)
  };

  sunPatch.clear();
  sunPatch.reserve(10);

  if (ttt) {
    // Initial TTT patch:
    sunPatch.push_back(tri(1, 1, 0, directionsTTT[0]));
    sunPatch.push_back(tri(1, 0, 1, directionsTTT[1]));
    sunPatch.push_back(tri(1, 0, 2, directionsTTT[2]));
    sunPatch.push_back(tri(1, 1, 3, directionsTTT[2]));
    sunPatch.push_back(tri(1, 1, 4, directionsTTT[3]));
    sunPatch.push_back(tri(1, 0, 5, directionsTTT[4]));
    sunPatch.push_back(tri(1, 1, 6, directionsTTT[4]));
    sunPatch.push_back(tri(1, 1, 7, directionsTTT[5]));
    sunPatch.push_back(tri(1, 0, 8, directionsTTT[6]));
    sunPatch.push_back(tri(1, 0, 9, directionsTTT[0]));
  } else {
    // Initial PRT patch:
    sunPatch.push_back(tri(1, 1, 0, directionsPRT[0]));
    sunPatch.push_back(tri(1, 0, 9, directionsPRT[0]));
    sunPatch.push_back(tri(1, 0, 1, directionsPRT[1]));
    sunPatch.push_back(tri(1, 1, 2, directionsPRT[1]));
    sunPatch.push_back(tri(1, 0, 3, directionsPRT[2]));
    sunPatch.push_back(tri(1, 1, 4, directionsPRT[2]));
    sunPatch.push_back(tri(1, 0, 5, directionsPRT[3]));
    sunPatch.push_back(tri(1, 1, 6, directionsPRT[3]));
    sunPatch.push_back(tri(1, 0, 7, directionsPRT[4]));
    sunPatch.push_back(tri(1, 1, 8, directionsPRT[4]));
  }
}

void TuebingenTri::minmax(const trilist& patch, vec4s& min, vec4s& max) {
  min.set(0, 0, 0, 0);
  max.set(0, 0, 0, 0);

  for (trilist::const_iterator i = patch.begin(); i != patch.end(); ++i) {
    const vec4s& ref = i->getRef();
    for (uint j = 0; j < 4; ++j) {
      if (ref[j] > max[j]) {
        max[j] = ref[j];
      }

      if (ref[j] < min[j]) {
        min[j] = ref[j];
      }
    }
  }
}

void TuebingenTri::createVerticesVis(Common::vec2dlist& vertices,
             const trilist& initial, uint steps, bool sector, bool ttt) {
  if (initial.empty())
    return;

  trilist* patch = new trilist;
  uint numtris;
  
  if (sector) {
    trilist initialsector;
    initialsector.push_back(initial[0]);
    numtris = countTris(initialsector, steps);

    iterate(initialsector, steps, ttt, *patch);
  } else {
    numtris = countTris(initial, steps);

    iterate(initial, steps, ttt, *patch);
  }

  VisList* vlist = new VisList;
  vlist->reserve(double(numtris) * 0.4);

  vlist->init();

  for (trilist::const_iterator i = patch->begin(); i != patch->end(); ++i) {
    vec4s temp[3];
    i->getVertices(temp);

    for (uint j = 0; j < 3; ++j) {
      // never include zero, the reference point for the radial projection
      if (temp[j].isZero()) continue;

      vlist->insertSorted(temp[j]);
    }
  }

  cerr << "statistics: " << patch->size() << " triangles reduced to "
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

void TuebingenTri::createVerticesVis2(Common::vec4slist& vertices, const trilist& initial,
             uint steps, bool sector, bool ttt) {
  if (initial.empty())
    return;

  trilist* patch = new trilist;
  uint numtris;
  
  if (sector) {
    trilist initialsector;
    initialsector.push_back(initial[0]);
    numtris = countTris(initialsector, steps);

    iterate(initialsector, steps, ttt, *patch);
  } else {
    numtris = countTris(initial, steps);

    iterate(initial, steps, ttt, *patch);
  }

  VisList* vlist = new VisList;
  vlist->reserve(double(numtris) * 0.4);

  vlist->init();

  for (trilist::const_iterator i = patch->begin(); i != patch->end(); ++i) {
    vec4s temp[3];
    i->getVertices(temp);

    for (uint j = 0; j < 3; ++j) {
      // never include zero, the reference point for the radial projection
      if (temp[j].isZero()) continue;

      vlist->insertSorted(temp[j]);
    }
  }

  cerr << "statistics: " << patch->size() << " triangles reduced to "
       << vlist->size() << " unique vertices." << endl;

  delete patch;
  patch = NULL;

  vlist->removeInvisibleProper();

  vertices.clear();
  vertices.reserve(vlist->size());
  vlist->dump(vertices);

  delete vlist;
  vlist = NULL;
}

int main(int argc, char* argv[]) {
  uint steps = 40;
  uint mode = 0;
  bool sector = false;

  using namespace Common;
  
  if (argc >= 2) {
    stringstream ss(argv[1]);
    ss >> steps;
  }

  if (argc >= 3) {
    stringstream ss(argv[2]);
    ss >> mode;
  }

  if (argc >= 4) {
    stringstream ss(argv[3]);
    ss >> sector;
  }

  const bool ttt = (mode >= 0 && mode <= 2);

  // Construct initial patch, this one is needed for all modes
  TuebingenTri::trilist initialTueb;
  TuebingenTri::constructSun(initialTueb, ttt);

  // First 3 modes are for the TTT. The other (last) 3 modes are for the PRT.
  switch (mode) {

    // Output triangle data in Mathematica style (0 = TTT / 3 = PRT)
    case 0: // fallthrough!
    case 3: {
              TuebingenTri::trilist output;

              vec4s min, max;
              iterate(initialTueb, steps, ttt, output);
              minmax(output, min, max);

              cerr << "statistics: min = " << min
                   << ", max = " << max << endl;

              cout << output;
            }
            break;

    // Output visible vertices of patch in Mathematica style (1 = TTT / 4 = PRT)
    // Vertices are output as integer 4-tuples (_not_ as double 2-tuples / cartesian)
    case 1: // fallthrough!
    case 4: {
              vec4slist output;
              TuebingenTri::createVerticesVis2(output, initialTueb, steps, sector, ttt);

              cout << output;
            }
            break;

    // Apply radial projection method and output histogram data (2 = TTT / 5 = PRT)
    case 2: // fallthrough!
    case 5: {
              vec2dlist verts;
              TuebingenTri::createVerticesVis(verts, initialTueb, steps, sector, ttt);

              dlist output;
              double mean;
              radialProj(verts, output, mean);

              cerr << "mean distance " << mean
                   << " during radial projection of "
                   << verts.size() << " vertices.\n";
              writeRawConsole(output);
            }
            break;

    default: cerr << "error: unsupported mode selected!\n";
             return 0;
  }

  return 0;
}

