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

#include "higher_cyclo.h"

#include <sstream>
#include "level_manager.h"

#include <algorithm>

namespace Common {

  bool locate(const vec6slist& list, const vec6s& target,
              uint i, uint j, uint range);

};

namespace Heptagonal {

  // Minimal polynomial is p(x) = x^3 + x^2 - 2*x - 1
  const double lambda = 2.0 * cos(2.0 * Constants::pi / 7.0);

  // Use a ball of radius R in 4-space.
  const double refBallRadiusSquared = 6.5;

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

namespace Elevenfold {

  // Minimal polynomial is p(x) = x^5 + x^4 - 4*x^3 - 3*x^2 + 3*x + 1
  const double lambda = 2.0 * cos(2.0 * Constants::pi / 11.0);

  // Use a ball of radius R in 8-space
  const double refBallRadiusSquared = 10.5;

  const double VisOp::epsilon = 2.0 * numeric_limits<double>::epsilon();

  bool checkProjInWindow(const vec10s& point);

  void projTiling(const vec10s& initpoint, uint maxstep,
                  Common::vec10slist& tilingpoints);

  void projTilingVis(const vec10s& initpoint,
                     const vec10s& origin,
                     uint maxstep, bool radialproj,
                     Common::vec10slist& tilingpoints,
                     Common::vec10slist& visiblepoints);

  void radialProj(const Common::vec10slist& input,
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

  if (refBallRadiusSquared - pt1 > Constants::eps) {
    return true;
  } else {
    if (refBallRadiusSquared - pt1 < -Constants::eps) return false;
  }

  cerr << "Warning: Insufficient accuracy in function checkProjInWindow.\n";
  return false;
}

void Heptagonal::projTiling(const vec6s& initpoint, uint maxstep, 
                   Common::vec6slist& tilingpoints) {
  using namespace Common;

  vec6s p, pp[2];
  const uint numsteps = 7;

  tilingpoints.clear();
  tilingpoints.push_back(initpoint);

  if (!checkProjInWindow(initpoint)) {
    cerr << "Initial point not in projection window.\n";
    return;
  }

  TVLManager<vec6slist, 2 + 1> lvlman(tilingpoints);

  for (uint n = 0; n < maxstep; ++n) {
    for (uint i = lvlman.begin(); i < lvlman.end(); ++i) {
      p = tilingpoints[i];

      for (uint j = 0; j < numsteps; ++j) {
        pp[0] = p.step(j, false);
        pp[1] = p.step(j, true);

        if (checkProjInWindow(pp[0])) lvlman.insert(pp[0]);
        if (checkProjInWindow(pp[1])) lvlman.insert(pp[1]);
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

  vec6s p, pp[2];
  const uint numsteps = 7;

  tilingpoints.clear();
  tilingpoints.push_back(initpoint);

  if (!checkProjInWindow(initpoint)) {
    cerr << "Initial point not in projection window.\n";
    return;
  }

  TVLManager<vec6slist, 2 + 1> lvlman(tilingpoints);

  for (uint n = 0; n < maxstep; ++n) {
    for (uint i = lvlman.begin(); i < lvlman.end(); ++i) {
      p = tilingpoints[i];

      for (uint j = 0; j < numsteps; ++j) {
        pp[0] = p.step(j, false);
        pp[1] = p.step(j, true);

        if (checkProjInWindow(pp[0])) lvlman.insert(pp[0]);
        if (checkProjInWindow(pp[1])) lvlman.insert(pp[1]);
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
    cerr << "info: processing heptagonal tiling in radial projection mode.\n";

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

bool Elevenfold::checkProjInWindow(const vec10s& point) {
  using namespace Common;

  const vec8d pt(point.toInternalL11());
  const double pt1 = pt.lengthSquared();

  if (refBallRadiusSquared - pt1 > Constants::eps) {
    return true;
  } else {
    if (refBallRadiusSquared - pt1 < -Constants::eps) return false;
  }

  cerr << "Warning: Insufficient accuracy in function checkProjInWindow.\n";
  return false;
}

void Elevenfold::projTiling(const vec10s& initpoint, uint maxstep,
                  Common::vec10slist& tilingpoints) {
  using namespace Common;

  vec10s p, pp[2];
  const uint numsteps = 11;

  tilingpoints.clear();
  tilingpoints.push_back(initpoint);

  if (!checkProjInWindow(initpoint)) {
    cerr << "Initial point not in projection window.\n";
    return;
  }

  TVLManager<vec10slist, 2 + 1> lvlman(tilingpoints);

  for (uint n = 0; n < maxstep; ++n) {
    for (uint i = lvlman.begin(); i < lvlman.end(); ++i) {
      p = tilingpoints[i];

      for (uint j = 0; j < numsteps; ++j) {
        pp[0] = p.step(j, false);
        pp[1] = p.step(j, true);

        if (checkProjInWindow(pp[0])) lvlman.insert(pp[0]);
        if (checkProjInWindow(pp[1])) lvlman.insert(pp[1]);
      }
    }

    lvlman.advance();
  }

  cerr << "Constructed patch of elevenfold tiling with "
       << tilingpoints.size() << " vertices.\n";
}

void Elevenfold::projTilingVis(const vec10s& initpoint,
                     const vec10s& origin,
                     uint maxstep, bool radialproj,
                     Common::vec10slist& tilingpoints,
                     Common::vec10slist& visiblepoints) {
  using namespace Common;

  vec10s p, pp[2];
  const uint numsteps = 11;

  tilingpoints.clear();
  tilingpoints.push_back(initpoint);

  if (!checkProjInWindow(initpoint)) {
    cerr << "Initial point not in projection window.\n";
    return;
  }

  TVLManager<vec10slist, 2 + 1> lvlman(tilingpoints);

  for (uint n = 0; n < maxstep; ++n) {
    for (uint i = lvlman.begin(); i < lvlman.end(); ++i) {
      p = tilingpoints[i];

      for (uint j = 0; j < numsteps; ++j) {
        pp[0] = p.step(j, false);
        pp[1] = p.step(j, true);

        if (checkProjInWindow(pp[0])) lvlman.insert(pp[0]);
        if (checkProjInWindow(pp[1])) lvlman.insert(pp[1]);
      }
    }

    lvlman.advance();
  }

  cerr << "Constructed patch of elevenfold tiling with "
       << tilingpoints.size() << " vertices.\n";

  VisList* vlist = new VisList;

  // For RP-mode with the default origin, we can apply the usual
  // sector-reduction optimization.
  if (radialproj && origin.isZero()) {
    assert((tilingpoints.size() - 1) % 11 == 0);
    vlist->reserve((tilingpoints.size() - 1) / 11);
  } else {
    vlist->reserve(tilingpoints.size() - 1);
  }

  vlist->init();

  if (radialproj) {
    cerr << "info: processing elevenfold tiling in radial projection mode.\n";

    // If we use the default origin, we can reduce the patch to the usual sector.
    if (origin.isZero()) {
      for (vec10slist::const_iterator i = tilingpoints.begin(); i != tilingpoints.end(); ++i) {
        if (i->isZero()) continue;

        const vec2d phys(i->toPhysicalL11());

        if (phys.inFirstQuadOpen() && phys.inSectorL11()) {
          vlist->insertSorted(*i);
        }
      }
    } else {
      for (vec10slist::const_iterator i = tilingpoints.begin(); i != tilingpoints.end(); ++i) {
        const vec10s shifted(*i - origin);
        if (!shifted.isZero()) vlist->insertSorted(shifted);
      }
    }

    vlist->removeInvisibleFast();
  } else {
    for (vec10slist::const_iterator i = tilingpoints.begin(); i != tilingpoints.end(); ++i) {
      const vec10s shifted(*i - origin);

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

void Elevenfold::radialProj(const Common::vec10slist& input,
                  Common::dlist& output, double& meandist) {
  // TODO: implement
}

ostream& operator<<(ostream &os, const vec4d& v) {
  os << '{' << fixed << setprecision(3) << v[0];
  for (uint i = 1; i < 4; ++i)
    os << ',' << v[i];
  os << '}';

  return os;
}

ostream& operator<<(ostream &os, const vec8d& v) {
  os << '{' << fixed << setprecision(3) << v[0];
  for (uint i = 1; i < 8; ++i)
    os << ',' << v[i];
  os << '}';

  return os;
}

void print_usage() {
  cerr << "higher_cyclo: usage:" << endl;

  cerr << "higher_cyclo --hepta: selects heptagonal (n=7) tiling" << endl;
  cerr << "higher_cyclo --eleven: selects elevenfold tiling (n=11)" << endl;

  cerr << "\tparameter 1: mode" << endl;
    cerr << "\t\t" << "0 = tiling vertices" << endl;
    cerr << "\t\t" << "1 = visible vertices" << endl;
    cerr << "\t\t" << "2 = radial projection" << endl;
  cerr << "\tparameter 2: steps" << endl;
  cerr << "\tparameter 3: sector (if possible reduce point set to sector)" << endl;

  cerr << endl;
}

int main_hepta(int argc, char* argv[]) {
  const vec6s init(0, 0, 0, 0, 0, 0);

  stringstream parser;

  uint steps = 20;
  uint mode = 0;

  using namespace Heptagonal;

  if (argc >= 2) {
    parser.str(argv[1]);
    parser.clear();
    parser >> steps;
  }

  if (argc >= 3) {
    parser.str(argv[2]);
    parser.clear();
    parser >> mode;
  }

  Common::vec6slist tiling, visible;
  Common::dlist spacings;
  double mean;

  if (mode >= 3) {
    cerr << "error: unknown mode " << mode << " selected for heptagonal tiling.\n";
    return 1;
  }

  switch (mode) {
  case 0:
    projTiling(init, steps, tiling);
    cout << tiling;
    break;

  case 1:
    projTilingVis(init, init, steps, false, tiling, visible);
    cout << visible;
    break;

  case 2:
    projTilingVis(init, init, steps, true, tiling, visible);
    radialProj(visible, spacings, mean);
    Common::meanDistanceMessage(spacings.size() + 1, mean);
    Common::writeRawConsole(spacings);
    break;

  default:
    assert(false);
    return 1;
  }

  return 0;
}

int main_eleven(int argc, char* argv[]) {
  const vec10s init(0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

  stringstream parser;

  uint steps = 10;
  uint mode = 0;

  using namespace Elevenfold;

  if (argc >= 2) {
    parser.str(argv[1]);
    parser.clear();
    parser >> steps;
  }

  if (argc >= 3) {
    parser.str(argv[2]);
    parser.clear();
    parser >> mode;
  }

  Common::vec10slist tiling, visible;
  Common::dlist spacings;
  double mean;

  if (mode >= 3) {
    cerr << "error: unknown mode " << mode << " selected for elevenfold tiling.\n";
    return 1;
  }

  switch (mode) {
  case 0:
    projTiling(init, steps, tiling);
    cout << tiling;
    break;

  case 1:
    projTilingVis(init, init, steps, false, tiling, visible);
    cout << visible;
    break;

  case 2:
    projTilingVis(init, init, steps, true, tiling, visible);
    radialProj(visible, spacings, mean);
    Common::meanDistanceMessage(spacings.size() + 1, mean);
    Common::writeRawConsole(spacings);
  break;

  default:
    assert(false);
    return 1;
  }

  return 0;
}

int main(int argc, char* argv[]) {
  stringstream parser;
  string main_mode;

  if (argc >= 2) {
    parser.str(argv[1]);
    parser.clear();
    parser >> main_mode;
  }

  if (main_mode == "--hepta")
    return main_hepta(argc - 1, argv + 1);
  else if (main_mode == "--eleven")
    return main_eleven(argc - 1, argv + 1);
  else
    print_usage();

  return 0;
}
