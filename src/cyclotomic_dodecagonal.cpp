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

#include "cyclotomic_dodecagonal.h"

#include "level_manager.h"

#include <algorithm>

const double Dodecagonal::VisOp::epsilon = 2.0 * numeric_limits<double>::epsilon();

bool Dodecagonal::checkProjInSector(const vec2d& orthpoint, bool useAlt) {
  using namespace Common;

  const vec2d v(orthpoint.reduceIntoSectorL12());
  double test;

  const vec2d* const verts = (useAlt ? verticesAlt : vertices);

  for (uint i = 0; i < 2; ++i) {
    test = checkPosition(verts[i], verts[i+1], v);
    if (test < -Constants::eps) return false;
#ifndef NDEBUG
    if (test <= Constants::eps) {
      cerr << "Warning: Insufficient accuracy in function checkProjInSector.\n";
    }
#endif
  }

  return true;
}

bool Dodecagonal::checkProjInWindow(const vec4i& point, bool useCircle) {
  using namespace Common;

  const vec2d pt(point.orthProjShiftL12());
  const double pt1 = pt.lengthSquared();

  if (useCircle) {
    return circularCheck(refCircleRadiusSquared, pt1);
  }

  if (innerRadSquared - pt1 > Constants::eps) {
    return true;
  } else {
    if (outerRadSquared - pt1 < -Constants::eps) {
      return false;
    } else {
      return checkProjInSector(pt, windowBookOrientation);
    }
  }
}

bool Dodecagonal::checkScaledProjInWindow(const vec4i& point,
                    bool gcdNormTwo, bool useCircle) {
  using namespace Common;

  /* Empirical tests with large patches indicate that the test  *
   * against a rescaled window depends on the gcd-norm of the   *
   * vertex. The scaling factors look a bit weird, but become   *
   * nicer if we don't consider the dodecagon with outer radius *
   * sqrt(unitZ3) but the normalized dodecagon (outer rad = 1). */
  const double scaler = gcdNormTwo ? sqrt(Constants::unitZ3 * 0.5) :
                                     sqrt(Constants::unitZ3 * 2.0);

  /* Interesting observation: Multiply the two scaling factors and *
   * we recover unitZ3 again. Hmm, nice!                           */

  const vec2d pt(point.orthProjShiftL12(scaler, gcdNormTwo));
  const double pt1 = pt.lengthSquared();

  if (useCircle) {
    return circularCheck(refCircleRadiusSquared, pt1);
  }

  if (innerRadSquared - pt1 > Constants::eps) {
    return true;
  } else {
    if (outerRadSquared - pt1 < -Constants::eps) {
      return false;
    } else {
      return checkProjInSector(pt, windowBookOrientation);
    }
  }
}

void Dodecagonal::projTiling(const vec4i& initpoint, uint maxstep,
                     Common::vec4ilist& tilingpoints) {
  using namespace Common;

  vec4i p, pp;
  const uint numsteps = 12;
  const vec4i hyperstep[12] = {vec4i(1,0,0,0),  vec4i(0,1,0,0),
                               vec4i(0,0,1,0),  vec4i(0,0,0,1),
                               vec4i(-1,0,1,0), vec4i(0,-1,0,1),
                               vec4i(-1,0,0,0), vec4i(0,-1,0,0),
                               vec4i(0,0,-1,0), vec4i(0,0,0,-1),
                               vec4i(1,0,-1,0), vec4i(0,1,0,-1)};

  tilingpoints.clear();
  tilingpoints.push_back(initpoint);

  if (!checkProjInWindow(initpoint, circularWindow)) {
    cerr << "Initial point not in projection window.\n";
    return;
  }

  // We need 2 + 1 levels to avoid going "back" (into the wrong direction) when creating the patch.
  TVLManager<vec4ilist, 2 + 1> lvlman(tilingpoints);

  for (uint n = 0; n < maxstep; ++n) {
    for (uint i = lvlman.begin(); i < lvlman.end(); ++i) {
      p = tilingpoints[i];

      for (uint j = 0; j < numsteps; ++j) {
        pp = p + hyperstep[j];

        if (checkProjInWindow(pp, circularWindow)) lvlman.insert(pp);
      }
    }

    lvlman.advance();
  }

  cerr << "Constructed patch of dodecagonal tiling with "
       << tilingpoints.size() << " vertices.\n";
}

void Dodecagonal::projTilingVisLocal(const vec4i& initpoint, uint maxstep,
                     Common::vec4ilist& tilingpoints,
                     Common::vec4ilist& visiblepoints) {
  using namespace Common;
  using namespace Coprime;

  vec4i p, pp;
  const uint numsteps = 12;
  const vec4i hyperstep[12] = {vec4i(1,0,0,0),  vec4i(0,1,0,0),
                               vec4i(0,0,1,0),  vec4i(0,0,0,1),
                               vec4i(-1,0,1,0), vec4i(0,-1,0,1),
                               vec4i(-1,0,0,0), vec4i(0,-1,0,0),
                               vec4i(0,0,-1,0), vec4i(0,0,0,-1),
                               vec4i(1,0,-1,0), vec4i(0,1,0,-1)};

  tilingpoints.clear();
  visiblepoints.clear();

  tilingpoints.push_back(initpoint);

  if (!checkProjInWindow(initpoint, circularWindow)) {
    cerr << "Initial point not in projection window.\n";
    return;
  }

  // We need 2 + 1 levels to avoid going "back" (into the wrong direction) when creating the patch.
  TVLManager<vec4ilist, 2 + 1> lvlman(tilingpoints);

  for (uint n = 0; n < maxstep; ++n) {
    for (uint i = lvlman.begin(); i < lvlman.end(); ++i) {
      p = tilingpoints[i];

      for (uint j = 0; j < numsteps; ++j) {
        pp = p + hyperstep[j];

        if (!checkProjInWindow(pp, circularWindow)) continue;
        if (!lvlman.insert(pp)) continue;

        /* By empirical tests the only type of vertices that can be visible *
         * are those with a gcd-norm of 1 (coprime elements) or gcd-norm 2. *
         * The norm 2 vertices are the exceptional ones which don't appear  *
         * in the octagonal (AB) and decagonal (TTT) case.                  */
        const int gcdnorm = Coprime::gcdZ3(pp.transL12ToDirect()).normZ3();
        switch (gcdnorm) {
          case 1: if (checkScaledProjInWindow(pp, false, circularWindow)) continue;
                  break;
          case 2: if (checkScaledProjInWindow(pp, true, circularWindow)) continue;
                  break;
          default: continue;
        }

        visiblepoints.push_back(pp);
      }
    }

    lvlman.advance();
  }

  cerr << "Constructed patch of dodecagonal tiling with "
       << tilingpoints.size() << " vertices and "
       << visiblepoints.size() << " visible ones.\n";
}

void Dodecagonal::projTilingVis(const vec4i& initpoint,
                     const vec4i& origin,
                     uint maxstep, bool radialproj,
                     Common::vec4ilist& tilingpoints,
                     Common::vec4ilist& visiblepoints) {
  using namespace Common;

  vec4i p, pp;
  const uint numsteps = 12;
  const vec4i hyperstep[12] = {vec4i(1,0,0,0),  vec4i(0,1,0,0),
                               vec4i(0,0,1,0),  vec4i(0,0,0,1),
                               vec4i(-1,0,1,0), vec4i(0,-1,0,1),
                               vec4i(-1,0,0,0), vec4i(0,-1,0,0),
                               vec4i(0,0,-1,0), vec4i(0,0,0,-1),
                               vec4i(1,0,-1,0), vec4i(0,1,0,-1)};

  tilingpoints.clear();
  tilingpoints.push_back(initpoint);

  if (!checkProjInWindow(initpoint, circularWindow)) {
    cerr << "Initial point not in projection window.\n";
    return;
  }

  // We need 2 + 1 levels to avoid going "back" (into the wrong direction) when creating the patch.
  TVLManager<vec4ilist, 2 + 1> lvlman(tilingpoints);

  for (uint n = 0; n < maxstep; ++n) {
    for (uint i = lvlman.begin(); i < lvlman.end(); ++i) {
      p = tilingpoints[i];

      for (uint j = 0; j < numsteps; ++j) {
        pp = p + hyperstep[j];

        if (checkProjInWindow(pp, circularWindow)) lvlman.insert(pp);
      }
    }

    lvlman.advance();
  }

  cerr << "Constructed patch of dodecagonal tiling with "
       << tilingpoints.size() << " vertices.\n";

  extractVisible(origin, radialproj, tilingpoints, visiblepoints);
}

vec4i Dodecagonal::sqDist(const vec4i& v, const vec4i& w) {
  const vec4i z(v - w);

  return z.multL12(z.conjL12());
}

void Dodecagonal::projTilingPatch(const vec4i& initpoint, uint maxstep,
                     Common::vec4ilist& tilingpoints,
                     Common::edgelist& edges) {
  using namespace Common;

  vec4i p, pp;
  const uint numsteps = 12;
  const vec4i hyperstep[12] = {vec4i(1,0,0,0),  vec4i(0,1,0,0),
                               vec4i(0,0,1,0),  vec4i(0,0,0,1),
                               vec4i(-1,0,1,0), vec4i(0,-1,0,1),
                               vec4i(-1,0,0,0), vec4i(0,-1,0,0),
                               vec4i(0,0,-1,0), vec4i(0,0,0,-1),
                               vec4i(1,0,-1,0), vec4i(0,1,0,-1)};

  tilingpoints.clear();
  tilingpoints.push_back(initpoint);

  if (!checkProjInWindow(initpoint, circularWindow)) {
    cerr << "Initial point not in projection window.\n";
    return;
  }

  TVLManager<vec4ilist, 2 + 1> lvlman(tilingpoints);

  for (uint n = 0; n < maxstep; ++n) {
    for (uint i = lvlman.begin(); i < lvlman.end(); ++i) {
      p = tilingpoints[i];

      for (uint j = 0; j < numsteps; ++j) {
        pp = p + hyperstep[j];

        if (checkProjInWindow(pp, circularWindow)) lvlman.insert(pp);
      }
    }

    lvlman.advance();
  }

  cerr << "Constructed patch of dodecagonal tiling with "
       << tilingpoints.size() << " vertices.\n";

  // Connect the vertices with edgelength Sqrt[2 - Sqrt[3]]
  edges.clear();
  const vec4i edgedist(vec4i(2, -1, 0, 0).transDirectToL12());

  for (uint i = 0; i < tilingpoints.size(); ++i) {
    for (uint j = i + 1; j < tilingpoints.size(); ++j) {
      if (sqDist(tilingpoints[i], tilingpoints[j]) != edgedist) continue;

      edges.push_back(tilingEdge(i, j));
    }
  }

  cerr << "Introduced " << edges.size() << " edges to the patch.\n"; 
}

void Dodecagonal::extractSector(const Common::vec4ilist& input,
                     Common::vec4ilist& output) {

  using namespace Common;

  // Warning: Due to the shift into a non-singular position, there
  // is no 12-star at the origin of the tiling.
  output.clear();
  output.reserve(input.size() / 6);

  for (vec4ilist::const_iterator i = input.begin(); i != input.end(); ++i) {
    const vec2d physProj(i->paraProjL12());

    if (physProj.inFirstQuadrant() && physProj.inSectorL12()) {
      output.push_back(*i);
    }
  }
}

void Dodecagonal::extractVisible(const vec4i& origin, bool radialproj,
                      const Common::vec4ilist& input,
                      Common::vec4ilist& output) {
  using namespace Common;

  VisList* vlist = new VisList;
  vlist->reserve(input.size() - 1);

  vlist->init();

  for (vec4ilist::const_iterator i = input.begin(); i != input.end(); ++i) {
    const vec4i shifted(*i - origin);

    if (shifted.isZero()) continue; // zero = reference point (not visible)
    vlist->insertSorted(shifted);
  }

  if (radialproj)
    vlist->removeInvisibleFast();
  else
    vlist->removeInvisibleProper();

  output.clear();
  output.reserve(vlist->size());
  vlist->dump(output);

  delete vlist;
}


void Dodecagonal::extractVisibleFast(const vec4i& origin,
                const Common::vec4ilist& input, Common::vec4ilist& output) {
  using namespace Common;

  if (input.empty())
    return;

  vec4ilist temp;
  temp.reserve(input.size());

  for (vec4ilist::const_iterator i = input.begin(); i != input.end(); ++i) {
    const vec4i shifted(*i - origin);

    if (shifted.isZero()) continue;
    temp.push_back(i->transL12ToDirect().directL12ToUnique());
  }

  sort(temp.begin(), temp.end());
  temp.erase(unique(temp.begin(), temp.end()), temp.end());

  cerr << "statistics: after (incorrect) visibility computation: "
       << temp.size() << " vertices visible." << endl;

  output.clear();
  output.reserve(temp.size());

  for (vec4ilist::const_iterator i = temp.begin(); i != temp.end(); ++i)
    output.push_back(i->transDirectToL12());
}

uint Dodecagonal::estimateGrowth(uint input, bool steps) {
  const double x = double(input);

  if (steps) {
    static const double params[4] = {0.051852, -10.6339, 113.08, 38.5713};

    return uint(params[0] * (params[1] + sqrt(params[2] + params[3] * x)));
  } else {
    // Parameters were computed using Mathematica's FindFit
    static const double params[2] = {10.6339, 9.64282};

    return uint(params[0] * x + params[1] * x * x);
  }
}

void Dodecagonal::radialProj(const Common::vec4ilist& input,
                  Common::dlist& output,
                  double& meandist) {
  using namespace Common;

  output.clear();
  output.reserve(input.size());

  dlist angles;
  angles.reserve(input.size());

  for (vec4ilist::const_iterator i = input.begin(); i != input.end(); ++i) {
    const vec2d physProj(i->paraProjL12());
    angles.push_back(physProj.angle());
  }

  sort(angles.begin(), angles.end());
  neighbourDiff(angles, output, meandist);
  normalizeAngDists(output, meandist);
}

void Dodecagonal::testWindow(Common::vec2ilist& output, uint resolution) {
  using namespace Common;

  // Scan the area [-3,3]^2 (outer dodecagon radius is approx. 1.93)
  const double step = 6.0 / double(resolution);

  output.clear();
  output.reserve(resolution * resolution);

  for (uint i = 0; i < resolution; ++i) {
    for (uint j = 0; j < resolution; ++j) {
      const vec2d pos(-3.0 + i * step, -3.0 + j * step);

      if (checkProjInSector(pos, windowBookOrientation)) {
        output.push_back(vec2i(i, j));
      }
    }
  }
}

bool Dodecagonal::VisOp::rayTest(const vec4i& a, const vec4i& b) {
  // let Z<3> be Z[sqrt[3]]
  // transform into the Z<3>*1 + Z<3>*xi
  // representation (this is a direct sum)
  const vec4i pa(a.transL12ToDirect());
  const vec4i pb(b.transL12ToDirect());

  // first filter the trivial cases
  if (pa.isFirstZero()) {
    return pb.isFirstZero();
  }

  if (pb.isFirstZero()) {
    return pa.isFirstZero();
  }

  if (pa.isSecondZero()) {
    return pb.isSecondZero();
  }

  if (pb.isSecondZero()) {
    return pa.isSecondZero();
  }

  // pa = z_a + w_a * xi
  // pb = z_b + w_b * xi
  // with z_a, z_b, w_a, w_b elements in Z<3>
  vec2i c, d;

  // now compute:
  // c = z_a * w_b
  // d = z_b * w_a
  Coprime::multZ3(vec2i(pa[0], pa[1]),
                  vec2i(pb[2], pb[3]), c);
  Coprime::multZ3(vec2i(pb[0], pb[1]),
                  vec2i(pa[2], pa[3]), d);

  return (c == d);
}

