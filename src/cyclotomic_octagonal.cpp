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

#include "cyclotomic_octagonal.h"

#include "level_manager.h"

#include <algorithm>

const double Octagonal::VisOp::epsilon = 2.0 * numeric_limits<double>::epsilon();

bool Octagonal::checkProjInSector(const vec2d& orthpoint) {
  using namespace Common;

  const double v[2] = {abs(orthpoint.x), abs(orthpoint.y)};

  const double tt = min(min(sqrt((2.0 * silverMean + 1.0) / 8.0) - v[0],
                            sqrt((2.0 * silverMean + 1.0) / 8.0) - v[1]),
                            silverMean / 2.0 - (v[0] + v[1]));

  if (tt > eps) {
    return true;
  } else {
    if (tt < -eps) {
      return false;
    }
  }

  cerr << "Warning: Insufficient accuracy in function checkProjInSector.\n";
  return false;
}

bool Octagonal::checkProjInWindow(const vec4i& point, bool useCircle) {
  using namespace Common;

  const vec2d pt(point.orthProjL8());
  const double pt1 = pt.lengthSquared();

  if (useCircle) {
    return circularCheck(refCircleRadiusSquared, pt1);
  }

  if (innerRadiusSquared - pt1 > eps) {
    return true;
  } else {
    if (outerRadiusSquared - pt1 < -eps) {
      return false;
    } else {
      return checkProjInSector(pt);
    }
  }
}

bool Octagonal::checkScaledProjInWindow(const vec4i& point, bool useCircle) {
  using namespace Common;

  const vec2d pt(point.orthProjL8() * silverMean);
  const double pt1 = pt.lengthSquared();

  if (useCircle) {
    return circularCheck(refCircleRadiusSquared, pt1);
  }

  if (innerRadiusSquared - pt1 > eps) {
    return true;
  } else {
    if (outerRadiusSquared - pt1 < -eps) {
      return false;
    } else {
      return checkProjInSector(pt);
    }
  }
}

void Octagonal::projTiling(const vec4i& initpoint, uint maxstep,
                 Common::vec4ilist& tilingpoints) {
  using namespace Common;

  vec4i p, pp;
  const uint numsteps = 8;
  const vec4i hyperstep[8] = {vec4i(1,0,0,0),  vec4i(0,1,0,0),
                              vec4i(0,0,1,0),  vec4i(0,0,0,1),
                              vec4i(-1,0,0,0), vec4i(0,-1,0,0),
                              vec4i(0,0,-1,0), vec4i(0,0,0,-1)};

  tilingpoints.clear();
  tilingpoints.push_back(initpoint);

  if (!checkProjInWindow(initpoint, circularWindow)) {
    cerr << "Initial point not in projection window.\n";
    return;
  }

  TVLManager<vec4i> lvlman(2 + 1, tilingpoints);

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

  cerr << "Constructed patch of octagonal tiling with "
       << tilingpoints.size() << " vertices.\n";
}

void Octagonal::projTilingVis(const vec4i& initpoint,
                 const vec4i& origin,
                 uint maxstep, bool radialproj,
                 Common::vec4ilist& tilingpoints,
                 Common::vec4ilist& visiblepoints) {
  using namespace Common;

  vec4i p, pp;
  const uint numsteps = 8;
  const vec4i hyperstep[8] = {vec4i(1,0,0,0),  vec4i(0,1,0,0),
                              vec4i(0,0,1,0),  vec4i(0,0,0,1),
                              vec4i(-1,0,0,0), vec4i(0,-1,0,0),
                              vec4i(0,0,-1,0), vec4i(0,0,0,-1)};

  tilingpoints.clear();
  visiblepoints.clear();

  tilingpoints.push_back(initpoint);

  if (!checkProjInWindow(initpoint, circularWindow)) {
    cerr << "Initial point not in projection window.\n";
    return;
  }

  // We need 2 + 1 levels to avoid going "back" (into the wrong direction) when creating the patch.
  TVLManager<vec4i> lvlman(2 + 1, tilingpoints);

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

  cerr << "Constructed patch of octagonal tiling with "
       << tilingpoints.size() << " vertices.\n";

  extractVisible(origin, radialproj, tilingpoints, visiblepoints);
}

// See projTilingVis for more comments on the code.
void Octagonal::projTilingVisLocal(const vec4i& initpoint,
                 uint maxstep, bool onlySector,
                 Common::vec4ilist& tilingpoints,
                 Common::vec4ilist& visiblepoints) {
  using namespace Common;
  using namespace Coprime;

  // Only use the hypersteps in positive directions.
  // This results in an (eighth) sector after projection into physical space.
  const uint numsteps = onlySector ? 4 : 8;

  vec4i p, pp;
  const vec4i hyperstep[8] = {vec4i(1,0,0,0),  vec4i(0,1,0,0),
                              vec4i(0,0,1,0),  vec4i(0,0,0,1),
                              vec4i(-1,0,0,0), vec4i(0,-1,0,0),
                              vec4i(0,0,-1,0), vec4i(0,0,0,-1)};

  tilingpoints.clear();
  visiblepoints.clear();

  tilingpoints.push_back(initpoint);

  if (!checkProjInWindow(initpoint, circularWindow)) {
    cerr << "Initial point not in projection window.\n";
    return;
  }

  TVLManager<vec4i> lvlman(2 + 1, tilingpoints);

  for (uint n = 0; n < maxstep; ++n) {
    for (uint i = lvlman.begin(); i < lvlman.end(); ++i) {
      p = tilingpoints[i];

      for (uint j = 0; j < numsteps; ++j) {
        pp = p + hyperstep[j];

        if (!checkProjInWindow(pp, circularWindow)) continue;
        if (!lvlman.insert(pp)) continue;

        // Apply the local test for visibility.
        if (!checkScaledProjInWindow(pp, circularWindow) &&
            coprimeZ2(pp.transL8ToDirect())) {
           visiblepoints.push_back(pp);
        }
      }
    }

    lvlman.advance();
  }

  cerr << "Constructed patch of octagonal tiling with "
       << tilingpoints.size() << " vertices and "
       << visiblepoints.size() << " visible ones.\n";
}

void Octagonal::extractVisible(const vec4i& origin, bool radialproj,
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

uint Octagonal::estimateSize(uint maxstep) {
  // Parameters were computed using Mathematica's FindFit
  static const double params[2] = {2.23594, 2.34416};

  const double x = double(maxstep);
  return uint(params[0] * x + params[1] * x * x);
}

uint Octagonal::estimateSteps(uint size) {
  static const double params[4] = {0.213296, -2.23594, 4.99944, 9.37663};

  const double x = double(size);
  return uint(params[0] * (params[1] + sqrt(params[2] + params[3] * x)));
}

void Octagonal::radialProj(const Common::vec4ilist& input,
                          Common::dlist& output,
                          double& meandist, bool onlySector) {
  using namespace Common;

  output.clear();
  output.reserve(input.size());

  dlist angles;
  angles.reserve(input.size());

  // Do radial projection and (optionally) check if the
  // projection falls into the sector:
  for (vec4ilist::const_iterator i = input.begin(); i != input.end(); ++i) {
    const vec2d physProj(i->paraProjL8());

    if (onlySector) {
      if (!physProj.inSectorL8()) continue;
    }

    angles.push_back(physProj.angle());
  }

  sort(angles.begin(), angles.end());
  neighbourDiff(angles, output, meandist);
  normalizeAngDists(output, meandist);
}

void Octagonal::radialProj(const Common::vec4ilist& input,
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

    if (shifted.isZero() || shifted.paraProjL8().lengthSquared() > rSq)
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
    angles.push_back(i->paraProjL8().angle());
  }

  output.clear();
  output.reserve(visiblepoints.size());

  sort(angles.begin(), angles.end());
  neighbourDiff(angles, output, meandist);
  normalizeAngDists(output, meandist);
}

void Octagonal::testWindow(Common::vec2ilist& output,
                          uint resolution) {
  using namespace Common;

  // Scan the area [-2,2]^2 (outer octagon radius is approx. 0.92)
  const double step = 4.0 / double(resolution);

  output.clear();
  output.reserve(resolution * resolution);

  for (uint i = 0; i < resolution; ++i) {
    for (uint j = 0; j < resolution; ++j) {
      const vec2d pos(-2.0 + i * step, -2.0 + j * step);

      if (checkProjInSector(pos)) {
        output.push_back(vec2i(i, j));
      }
    }
  }
}

void Octagonal::innerOuterRadius(const Common::vec4ilist& tilingpoints,
                          double& inner, double& outer) {
  using namespace Common;

  double out = 0.0;

  for (vec4ilist::const_iterator i = tilingpoints.begin(); i != tilingpoints.end(); ++i) {
    const double lSq = i->paraProjL8().lengthSquared();

    if (lSq > out) out = lSq;
  }

  outer = sqrt(out);
  inner = cos(pi / 8.0) * outer;
}

bool Octagonal::VisOp::rayTest(const invectype& a, const invectype& b) {
  // let Z<2> be Z[sqrt[2]]
  // transform into the Z<2>*1 + Z<2>*xi
  // representation (this is a direct sum)
  const invectype pa(a.transL8ToDirect());
  const invectype pb(b.transL8ToDirect());

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
  // with z_a, z_b, w_a, w_b elements in Z<2>
  vec2i c, d;

  // now compute:
  // c = z_a * w_b
  // d = z_b * w_a
  Coprime::multZ2(vec2i(pa[0], pa[1]),
                  vec2i(pb[2], pb[3]), c);
  Coprime::multZ2(vec2i(pb[0], pb[1]),
                  vec2i(pa[2], pa[3]), d);

  return (c == d);
}

