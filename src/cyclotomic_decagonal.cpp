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

#include "cyclotomic_decagonal.h"

#include "level_manager.h"

#include <algorithm>

bool Decagonal::checkProjInSector(const vec2d& orthpoint, bool useAlt) {
  using namespace Common;

  const vec2d v(abs(orthpoint.x), abs(orthpoint.y));
  double test;

  const vec2d* const verts = (useAlt ? verticesAlt : vertices);

  for (uint i = 0; i < 3; ++i) {
    test = checkPosition(verts[i], verts[i+1], v);
    if (test < -eps) return false;
    if (test <= eps) {
      cerr << "Warning: Insufficient accuracy in function checkProjInSector.\n";
    }
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

  if (innerRadSquared - pt1 > eps) {
    return true;
  } else {
    if (outerRadSquared - pt1 < -eps) {
      return false;
    } else {
      return checkProjInSector(pt, windowBookOrientation);
    }
  }
}

bool Decagonal::checkScaledProjInWindow(const vec4i& point, bool useCircle) {
  using namespace Common;

  const vec2d pt(point.orthProjShiftL5(tau));
  const double pt1 = pt.lengthSquared();

  if (useCircle) {
    return circularCheck(refCircleRadiusSquared, pt1);
  }

  if (innerRadSquared - pt1 > eps) {
    return true;
  } else {
    if (outerRadSquared - pt1 < -eps) {
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

  cerr << "Constructed patch of decagonal tiling with "
       << tilingpoints.size() << " vertices.\n";
}

void Decagonal::projTilingVisLocal(const vec4i& initpoint, uint maxstep,
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

  TVLManager<vec4i> lvlman(2 + 1, tilingpoints);

  for (uint n = 0; n < maxstep; ++n) {
    for (uint i = lvlman.begin(); i < lvlman.end(); ++i) {
      p = tilingpoints[i];

      for (uint j = 0; j < numsteps; ++j) {
        pp = p + hyperstep[j];

        if (!checkProjInWindow(pp, circularWindow)) continue;
        if (!lvlman.insert(pp)) continue;

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

