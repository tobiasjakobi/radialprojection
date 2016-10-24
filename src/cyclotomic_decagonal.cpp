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

#include "cyclotomic_decagonal.h"

#include "level_manager.h"

#include <algorithm>

const double Decagonal::VisOp::epsilon = 2.0 * numeric_limits<double>::epsilon();

namespace {

  bool checkPhyInSector(const vec2d& phy){
    return (phy.inFirstQuadrant() && phy.inSectorL5());
  }

  bool checkPhyInSectorEps(const vec2d& phy){
    return (phy.inFirstQuadrant() && phy.inSectorL5Eps());
  }

};

bool Decagonal::checkProjInSector(const vec2d& orthpoint, bool useAlt) {
  using namespace Common;

  const vec2d v(orthpoint.abs());
  double test;

  const vec2d* const verts = (useAlt ? verticesAlt : vertices);

  for (uint i = 0; i < 3; ++i) {
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

bool Decagonal::checkProjInWindow(const vec4i& point, bool useCircle) {
  using namespace Common;

  const vec2d pt(point.orthProjShiftL5());
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

bool Decagonal::checkScaledProjInWindow(const vec4i& point, bool useCircle) {
  using namespace Common;

  const vec2d pt(point.orthProjShiftL5(Constants::unitGM, true)); /* invert the shift here */
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

void Decagonal::projTiling(const vec4i& initpoint, uint maxstep,
                     Common::vec4ilist& tilingpoints) {
  using namespace Common;

  /* see 'projTilingVisLocal' for comments */
  const uint numsteps = 10;

  vec4i p, pp;
  const vec4i hyperstep[10] = {vec4i(1,0,0,0),  vec4i(0,1,0,0),
                               vec4i(0,0,1,0),  vec4i(0,0,0,1),
                               vec4i(1,1,1,1),  vec4i(-1,0,0,0),
                               vec4i(0,-1,0,0), vec4i(0,0,-1,0),
                               vec4i(0,0,0,-1), vec4i(-1,-1,-1,-1)};

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

  cerr << "Constructed patch of decagonal tiling with "
       << tilingpoints.size() << " vertices.\n";
}

void Decagonal::projTilingVis(const vec4i& initpoint,
                     const vec4i& origin,
                     uint maxstep, bool radialproj,
                     Common::vec4ilist& tilingpoints,
                     Common::vec4ilist& visiblepoints) {
  using namespace Common;

  vec4i p, pp;
  const uint numsteps = 10;
  const vec4i hyperstep[10] = {vec4i(1,0,0,0),  vec4i(0,1,0,0),
                               vec4i(0,0,1,0),  vec4i(0,0,0,1),
                               vec4i(1,1,1,1),  vec4i(-1,0,0,0),
                               vec4i(0,-1,0,0), vec4i(0,0,-1,0),
                               vec4i(0,0,0,-1), vec4i(-1,-1,-1,-1)};

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

        if (checkProjInWindow(pp, circularWindow)) lvlman.insert(pp);
      }
    }

    lvlman.advance();
  }

  cerr << "Constructed patch of decagonal tiling with "
       << tilingpoints.size() << " vertices.\n";

  extractVisible(origin, radialproj, tilingpoints, visiblepoints);
}

void Decagonal::projTilingVisLocal(const vec4i& initpoint, uint maxstep,
                     bool sector,
                     Common::vec4ilist& tilingpoints,
                     Common::vec4ilist& visiblepoints) {
  using namespace Common;
  using namespace Coprime;

  // The trick based on reducing the used hypersteps to the ones in positive
  // direction apparantly doesn't work anymore for the decagonal case. At
  // least it doesn't produce a single sector anymore but just thins out
  // the resulting tiling.
  const uint numsteps = 10;

  vec4i p, pp;
  const vec4i hyperstep[10] = {vec4i(1,0,0,0),  vec4i(0,1,0,0),
                               vec4i(0,0,1,0),  vec4i(0,0,0,1),
                               vec4i(1,1,1,1),  vec4i(-1,0,0,0),
                               vec4i(0,-1,0,0), vec4i(0,0,-1,0),
                               vec4i(0,0,0,-1), vec4i(-1,-1,-1,-1)};

  tilingpoints.clear();
  visiblepoints.clear();

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

        if (!checkProjInWindow(pp, circularWindow)) continue;
        if (!lvlman.insert(pp)) continue;

        /* Apply the local test for visibility:                                              *
         * The test is similar to the AB one. It always again a test for coprimality of the  *
         * coordinates (in the respective module), and a (scaled) window test. However here  *
         * we are working with a shifted window. For the local visibility test to agree with *
         * the generic / ray-test one, we have to invert the direction of the shift.         */
        if (!checkScaledProjInWindow(pp, circularWindow) &&
            coprimeZTau(pp.transL5ToDirect())) {
           visiblepoints.push_back(pp);
        }
      }
    }

    lvlman.advance();
  }

  cerr << "Constructed patch of decagonal tiling with "
       << tilingpoints.size() << " vertices and "
       << visiblepoints.size() << " visible ones.\n";
}

void Decagonal::extractSector(const Common::vec4ilist& input,
                 Common::vec4ilist& output) {
  using namespace Common;

  output.clear();
  output.reserve(input.size() / 4);

  for (vec4ilist::const_iterator i = input.begin(); i != input.end(); ++i) {
    const vec2d physProj(i->paraProjL5());

    if (physProj.inFirstQuadrant() && physProj.inSectorL5()) {
      output.push_back(*i);
    }
  }
}

void Decagonal::extractVisible(const vec4i& origin, bool radialproj,
                          const Common::vec4ilist& input,
                          Common::vec4ilist& output) {
  using namespace Common;

  // We're not removing vertices in this case, so allocate the full amount.
  VisList* vlist = new VisList;
  vlist->reserve(input.size() - 1);

  vlist->init();

  for (vec4ilist::const_iterator i = input.begin(); i != input.end(); ++i) {
    const vec4i shifted(*i - origin);

    if (shifted.isZero()) continue;
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

void Decagonal::extractVisibleFast(const vec4i& origin,
                const Common::vec4ilist& input, Common::vec4ilist& output) {
  using namespace Common;

  if (input.empty())
    return;

  vec4ilist temp;
  temp.reserve(input.size());

  for (vec4ilist::const_iterator i = input.begin(); i != input.end(); ++i) {
    const vec4i shifted(*i - origin);

    if (shifted.isZero()) continue;
    temp.push_back(i->transL5ToDirect().directL5ToUnique());
  }

  sort(temp.begin(), temp.end());
  temp.erase(unique(temp.begin(), temp.end()), temp.end());

  cerr << "statistics: after (incorrect) visibility computation: "
       << temp.size() << " vertices visible." << endl;

  output.clear();
  output.reserve(temp.size());

  for (vec4ilist::const_iterator i = temp.begin(); i != temp.end(); ++i)
    output.push_back(i->transDirectToL5());
}

uint Decagonal::estimateGrowth(uint input, bool steps) {
  const double x = double(input);

  if (steps) {
    static const double params[4] = {0.117292, -5.94039, 35.2883, 17.0515};

    return uint(params[0] * (params[1] + sqrt(params[2] + params[3] * x)));
  } else {
    // Parameters were computed using Mathematica's FindFit
    static const double params[2] = {5.94039, 4.26288};

    return uint(params[0] * x + params[1] * x * x);
  }
}

void Decagonal::radialProj(const Common::vec4ilist& input,
                          Common::dlist& output,
                          double& meandist) {
  using namespace Common;

  output.clear();
  output.reserve(input.size());

  dlist angles;
  angles.reserve(input.size());

  for (vec4ilist::const_iterator i = input.begin(); i != input.end(); ++i) {
    const vec2d physProj(i->paraProjL5());
    angles.push_back(physProj.angle());
  }

  sort(angles.begin(), angles.end());
  neighbourDiff(angles, output, meandist);
  normalizeAngDists(output, meandist);
}

void Decagonal::radialProj(const Common::vec4ilist& input,
                   const vec4i& origin, double radius,
                   Common::dlist& output, double& meandist) {
  using namespace Common;

  const double rSq = radius * radius;
  vec4ilist visiblepoints;
  dlist angles;

  VisList* vlist = new VisList;
  vlist->reserve(input.size() - 1);

  vlist->init();

  for (vec4ilist::const_iterator i = input.begin(); i != input.end(); ++i) {
    const vec4i shifted(*i - origin);

    if (shifted.isZero() || shifted.paraProjL5().lengthSquared() > rSq)
      continue;

    vlist->insertSorted(shifted);
  }

  vlist->removeInvisibleFast();

  visiblepoints.reserve(vlist->size());
  vlist->dump(visiblepoints);

  delete vlist;
  vlist = NULL;

  angles.reserve(visiblepoints.size());

  for (vec4ilist::const_iterator i = visiblepoints.begin(); i != visiblepoints.end(); ++i) {
    angles.push_back(i->paraProjL5().angle());
  }

  output.clear();
  output.reserve(visiblepoints.size());

  sort(angles.begin(), angles.end());
  neighbourDiff(angles, output, meandist);
  normalizeAngDists(output, meandist);
}

void Decagonal::projTilingVisLocal2(const vec4i& initpoint, uint maxstep,
          Common::vec4ilist& tilingpoints, Common::vec4ilist& visiblepoints) {
  projTilingVisLocal(initpoint, maxstep, false, tilingpoints, visiblepoints);
}

void Decagonal::testWindow(Common::vec2ilist& output, uint resolution) {
  using namespace Common;

  // Scan the area [-2,2]^2 (outer decagon radius is approx. 1.37)
  const double step = 4.0 / double(resolution);

  output.clear();
  output.reserve(resolution * resolution);

  for (uint i = 0; i < resolution; ++i) {
    for (uint j = 0; j < resolution; ++j) {
      const vec2d pos(-2.0 + i * step, -2.0 + j * step);

      if (checkProjInSector(pos, windowBookOrientation)) {
        output.push_back(vec2i(i, j));
      }
    }
  }
}

void Decagonal::innerOuterRadius(const Common::vec4ilist& tilingpoints,
                        double& inner, double& outer) {
  using namespace Common;

  double out = 0.0;

  for (vec4ilist::const_iterator i = tilingpoints.begin(); i != tilingpoints.end(); ++i) {
    const double lSq = i->paraProjL5().lengthSquared();

    if (lSq > out) out = lSq;
  }

  outer = sqrt(out);
  inner = cos(Constants::pi / 10.0) * outer;
}

bool Decagonal::VisOp::rayTest(const vec4i& a, const vec4i& b) {
  // let Z<t> be Z[tau]
  // transform into the Z<t>*1 + Z<t>*xi
  // representation (this is a direct sum)
  const vec4i pa(a.transL5ToDirect());
  const vec4i pb(b.transL5ToDirect());

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
  // with z_a, z_b, w_a, w_b elements in Z<t>
  vec2i c, d;

  // now compute:
  // c = z_a * w_b
  // d = z_b * w_a
  Coprime::multZTau(vec2i(pa[0], pa[1]),
                    vec2i(pb[2], pb[3]), c);
  Coprime::multZTau(vec2i(pb[0], pb[1]),
                    vec2i(pa[2], pa[3]), d);

  return (c == d);
}

