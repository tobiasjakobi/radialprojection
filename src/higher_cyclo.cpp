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

#include "higher_cyclo.h"

#include <sstream>
#include "level_manager.h"

#include <algorithm>

namespace Common {

  typedef vector<vec6s> vec6slist;

  bool locate(const vec6slist& list, const vec6s& target,
              uint i, uint j, uint range);

};

namespace Heptagonal {

  // Minimal polynomial is p(x) = x^3 + x^2 - 2*x - 1
  const double lambda = 2.0 * cos(2.0 * Common::pi / 7.0);

  // Use a ball of radius R in 4-space
  const double refBallRadiusSquared = 7.25; // TODO: adjust

  const double VisOp::epsilon = 2.0 * numeric_limits<double>::epsilon();

  bool checkProjInWindow(const vec6s& point);
  //bool checkScaledProjInWindow(const vec6s& point);

  void projTiling(const vec6s& initpoint, uint maxstep,
                  Common::vec6slist& tilingpoints);

  void projTilingVis(const vec6s& initpoint,
                     const vec6s& origin,
                     uint maxstep, bool radialproj,
                     Common::vec6slist& tilingpoints,
                     Common::vec6slist& visiblepoints);

  void radialProj(const Common::vec6slist& input,
                  Common::dlist& output, double& meandist);

};

bool Common::locate(const vec6slist& list, const vec6s& target,
                          uint i, uint j, uint range) {
  if (j < i) return false;
  assert((j - i + 1) == range);

  for (int k = 0; k < int(j - i + 1); ++k) {
    if (list[i - 1 + k] == target) return true;
  }

  return false;
}

bool Heptagonal::checkProjInWindow(const vec6s& point) {
  using namespace Common;

  const vec4d pt(point.toInternalL7());
  const double pt1 = pt.lengthSquared();

  if (refBallRadiusSquared - pt1 > eps) {
    return true;
  } else {
    if (refBallRadiusSquared - pt1 < -eps) return false;
  }

  cerr << "Warning: Insufficient accuracy in function checkProjInWindow.\n";
  return false;
}

/*bool HeptagonalRadial::checkScaledProjInWindow(const vec6s& point) {
}*/

void Heptagonal::projTiling(const vec6s& initpoint, uint maxstep, 
                   Common::vec6slist& tilingpoints) {
  using namespace Common;

  vec6s p, pp;
  const uint numsteps = 7;
  const vec6s hyperstep[7] = {vec6s(1,0,0,0,0,0), vec6s(0,1,0,0,0,0),
                              vec6s(0,0,1,0,0,0), vec6s(0,0,0,1,0,0),
                              vec6s(0,0,0,0,1,0), vec6s(0,0,0,0,0,1),
                              vec6s(-1,-1,-1,-1,-1,-1)};

  tilingpoints.clear();
  tilingpoints.push_back(initpoint);

  if (!checkProjInWindow(initpoint)) {
    cerr << "Initial point not in projection window.\n";
    return;
  }

  TVLManager<vec6s> lvlman(7 + 1, tilingpoints);

  for (uint n = 0; n < maxstep; ++n) {
    for (uint i = lvlman.begin(); i < lvlman.end(); ++i) {
      p = tilingpoints[i];

      for (uint j = 0; j < numsteps; ++j) {
        pp = p + hyperstep[j];

        if (checkProjInWindow(pp)) lvlman.insert(pp);
      }
    }

    lvlman.advance();
  }

  cerr << "Constructed patch of heptagonal tiling with "
       << tilingpoints.size() << " vertices.\n";
}

void Heptagonal::projTilingVis(const vec6s& initpoint,
                   const vec6s& origin, uint maxstep, bool radialproj,
                   Common::vec6slist& tilingpoints,
                   Common::vec6slist& visiblepoints) {
  using namespace Common;

  vec6s p, pp;
  const uint numsteps = 7;
  const vec6s hyperstep[7] = {vec6s(1,0,0,0,0,0), vec6s(0,1,0,0,0,0),
                              vec6s(0,0,1,0,0,0), vec6s(0,0,0,1,0,0),
                              vec6s(0,0,0,0,1,0), vec6s(0,0,0,0,0,1),
                              vec6s(-1,-1,-1,-1,-1,-1)};

  tilingpoints.clear();
  tilingpoints.push_back(initpoint);

  if (!checkProjInWindow(initpoint)) {
    cerr << "Initial point not in projection window.\n";
    return;
  }

  TVLManager<vec6s> lvlman(7 + 1, tilingpoints);

  for (uint n = 0; n < maxstep; ++n) {
    for (uint i = lvlman.begin(); i < lvlman.end(); ++i) {
      p = tilingpoints[i];

      for (uint j = 0; j < numsteps; ++j) {
        pp = p + hyperstep[j];

        if (checkProjInWindow(pp)) lvlman.insert(pp);
      }
    }

    lvlman.advance();
  }

  cerr << "Constructed patch of heptagonal tiling with "
       << tilingpoints.size() << " vertices.\n";

  VisList* vlist = new VisList;

  // For RP-mode with the default origin, we can apply the usual
  // sector-reduction optimization.
  if (radialproj && origin.isZero()) {
    assert((tilingpoints.size() - 1) % 7 == 0);
    vlist->reserve((tilingpoints.size() - 1) / 7);
  } else {
    vlist->reserve(tilingpoints.size() - 1);
  }

  vlist->init();

  if (radialproj) {
    cerr << "info: processing heptagonal tiling in RP mode.\n";

    // If we use the default origin, we can reduce the patch to the usual sector.
    if (origin.isZero()) {
      for (vec6slist::const_iterator i = tilingpoints.begin(); i != tilingpoints.end(); ++i) {
        if (i->isZero()) continue;

        const vec2d phys(i->toPhysicalL7());

        if (phys.inFirstQuadOpen() && phys.inSectorL7()) {
          vlist->insertSorted(*i);
        }
      }
    } else {
      for (vec6slist::const_iterator i = tilingpoints.begin(); i != tilingpoints.end(); ++i) {
        const vec6s shifted(*i - origin);
        if (!shifted.isZero()) vlist->insertSorted(shifted);
      }
    }

    vlist->removeInvisibleFast();
  } else {
    for (vec6slist::const_iterator i = tilingpoints.begin(); i != tilingpoints.end(); ++i) {
      const vec6s shifted(*i - origin);

      if (!shifted.isZero()) vlist->insertSorted(shifted);
    }
    vlist->removeInvisibleProper();
  }

  visiblepoints.clear();
  visiblepoints.reserve(vlist->size());
  vlist->dump(visiblepoints);

  delete vlist;
  vlist = NULL;
}

void Heptagonal::radialProj(const Common::vec6slist& input,
                  Common::dlist& output, double& meandist) {
  using namespace Common;

  output.clear();
  output.reserve(input.size());

  dlist angles;
  angles.reserve(input.size());

  for (vec6slist::const_iterator i = input.begin(); i != input.end(); ++i) {
    const vec2d phys(i->toPhysicalL7());
    angles.push_back(phys.angle());
  }

  sort(angles.begin(), angles.end());
  neighbourDiff(angles, output, meandist);
  normalizeAngDists(output, meandist);
}

int main(int argc, char* argv[]) {
  const vec6s init(0, 0, 0, 0, 0, 0);
  uint steps = 10;
  uint mode = 0;

  if (argc >= 2) {
    stringstream ss(argv[1]);
    ss >> steps;
  }

  if (argc >= 3) {
    stringstream ss(argv[2]);
    ss >> mode;
  }

  Common::vec6slist tiling, visible;

  if (mode >= 3) {
    cerr << "Unknown mode selected:\n"
         << "0 = generate full tiling\n"
         << "1 = only visible points\n"
         << "2 = radial projection\n";
    return 0;
  }

  if (mode == 0) {
    Heptagonal::projTiling(init, steps, tiling);
    cout << tiling;
  } else
  if (mode == 1) {
    Heptagonal::projTilingVis(init, init, steps, false, tiling, visible);
    cout << visible;
  } else {
    Common::dlist out;
    double mean;

    Heptagonal::projTilingVis(init, init, steps, true, tiling, visible);
    Heptagonal::radialProj(visible, out, mean);

    Common::meanDistanceMessage(out.size() + 1, mean);
    Common::writeRawConsole(out);
  }

  return 0;
}

