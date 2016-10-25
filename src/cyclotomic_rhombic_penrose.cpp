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

#include "cyclotomic_rhombic_penrose.h"

#include "level_manager.h"
#include "cyclotomic_decagonal.h"

#include <algorithm>

/*
 * Beginning of anonymous namespace.
 */
namespace {

/*
 * The rhombic penrose case uses four windows W1, W2, W3, W4:
 * Let P be the convex hull of {1, xi, xi^2, xi^3, xi^4} with
 * xi = exp(2*pi*i/5), then let
 * W1 = P, W4 = -P, W3 = tau*P, W2 = -tau*P.
 */

// Pentagon radii:
const double innerRadius[2] = {
  0.5 * Constants::unitGM,
  0.25 * (Constants::unitGM + 1.0)
};
const double outerRadius[2] = {1.0, Constants::unitGM}; 

const double innerRadSquared[2] = {
  0.25 * (Constants::unitGM + 1.0),
  (7.0 + 3.0*sqrt(5.0)) / 8.0
};
const double outerRadSquared[2] = {1.0, Constants::unitGM*Constants::unitGM};

const vec2d vertices[4][5] = {
  {
    vec2d(1.0, 0.0),
    vec2d(0.5 * (Constants::unitGM - 1.0), 0.5 * sqrt(Constants::unitGM + 2.0)),
    vec2d(-0.5 * Constants::unitGM, 0.5 * sqrt(3.0 - Constants::unitGM)),
    vec2d(-0.5 * Constants::unitGM, -0.5 * sqrt(3.0 - Constants::unitGM)),
    vec2d(0.5 * (Constants::unitGM - 1.0), -0.5 * sqrt(Constants::unitGM + 2.0))
  }, {
    vec2d(0.5 * (Constants::unitGM + 1.0), -0.5 * sqrt(Constants::unitGM + 2.0)),
    vec2d(0.5 * (Constants::unitGM + 1.0), 0.5 * sqrt(Constants::unitGM + 2.0)),
    vec2d(-0.5, 0.5 * sqrt(4.0*Constants::unitGM + 3.0)),
    vec2d(-Constants::unitGM, 0.0),
    vec2d(-0.5, -0.5 * sqrt(4.0*Constants::unitGM + 3.0))
  }, {
    vec2d(Constants::unitGM, 0.0),
    vec2d(0.5, 0.5 * sqrt(4.0*Constants::unitGM + 3.0)),
    vec2d(-0.5 * (Constants::unitGM + 1.0), 0.5 * sqrt(Constants::unitGM + 2.0)),
    vec2d(-0.5 * (Constants::unitGM + 1.0), -0.5 * sqrt(Constants::unitGM + 2.0)),
    vec2d(0.5, -0.5 * sqrt(4.0*Constants::unitGM + 3.0))
  }, {
    vec2d(0.5 * Constants::unitGM, -0.5 * sqrt(3.0 - Constants::unitGM)),
    vec2d(0.5 * Constants::unitGM, 0.5 * sqrt(3.0 - Constants::unitGM)),
    vec2d(-0.5 * (Constants::unitGM - 1.0), 0.5 * sqrt(Constants::unitGM + 2.0)),
    vec2d(-1.0, 0.0),
    vec2d(-0.5 * (Constants::unitGM - 1.0), -0.5 * sqrt(Constants::unitGM + 2.0))
  }
};

void getInnerOuterSquared(double& inner, double& outer, uint window) {
  switch (window) {
    case 0:
    case 3:
      inner = innerRadSquared[0];
      outer = outerRadSquared[0];
      break;

    case 1:
    case 2:
      inner = innerRadSquared[1];
      outer = outerRadSquared[1];
      break;

    default:
      assert(false);
  }
}

bool checkProjInSector(const vec2d& orthpoint, uint window) {
  using namespace Common;

  const vec2d v(orthpoint.abs());
  double test;

  const vec2d* const verts = vertices[window];

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

bool checkProjInWindow(const vec4i& point, uint window) {
  using namespace Common;

  const vec2d pt(point.orthProjShiftL5());
  const double pt1 = pt.lengthSquared();

  double innerSquared, outerSquared;
  getInnerOuterSquared(innerSquared, outerSquared, window);

  if (innerSquared - pt1 > Constants::eps) {
    return true;
  } else {
    if (outerSquared - pt1 < -Constants::eps) {
      return false;
    } else {
      return checkProjInSector(pt, window);
    }
  }
}

bool checkScaledProjInWindow(const vec4i& point, uint window) {
  using namespace Common;

  // Note that scaling with -tau (MINUS!) is correct!
  // Usually this doesn't matter since the window is symmetric under
  // transformation with -1, but not here.
  const vec2d pt(point.orthProjShiftL5() * -Constants::unitGM);
  const double pt1 = pt.lengthSquared();

  double innerSquared, outerSquared;
  getInnerOuterSquared(innerSquared, outerSquared, window);

  if (innerSquared - pt1 > Constants::eps) {
    return true;
  } else {
    if (outerSquared - pt1 < -Constants::eps) {
      return false;
    } else {
      return checkProjInSector(pt, window);
    }
  }
}

struct VisOp {
  typedef Common::vec4ilist list_type;
  static const double epsilon;

  static inline double angle(const vec4i& a) {
    return a.paraProjL5().angle();
  }

  static inline vec2d toR2(const vec4i& a) {
    return a.paraProjL5();
  }

  static bool rayTest(const vec4i& a, const vec4i& b);
};

bool VisOp::rayTest(const vec4i& a, const vec4i& b) {
  // transform into the Z[tau]*1 + Z[tau]*xi
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
  // with z_a, z_b, w_a, w_b elements in Z[tau]
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

const double VisOp::epsilon = 2.0 * numeric_limits<double>::epsilon();

typedef VisTest::VisibleList<VisOp> VisList;

};

/*
 * End of anonymous namespace.
 */


void RhombicPenrose::projTiling(const vec4i& initpoint, uint maxstep,
             Common::vec4ilist& tilingpoints) {
  using namespace Common;

  vec4i p, pp;
  const uint numsteps = 10;
  const vec4i hyperstep[10] = {vec4i(1,0,0,0),  vec4i(0,1,0,0),
                               vec4i(0,0,1,0),  vec4i(0,0,0,1),
                               vec4i(1,1,1,1),  vec4i(-1,0,0,0),
                               vec4i(0,-1,0,0), vec4i(0,0,-1,0),
                               vec4i(0,0,0,-1), vec4i(-1,-1,-1,-1)};

  if (initpoint.kappaL5() != 0) {
    cerr << "Initial point not of zero-parity.\n";
    return;
  }

  tilingpoints.clear();
  tilingpoints.push_back(initpoint);

  TVLManager<vec4ilist, 2 + 1> lvlman(tilingpoints);

  for (uint n = 0; n < maxstep; ++n) {
    for (uint i = lvlman.begin(); i < lvlman.end(); ++i) {
      p = tilingpoints[i];

      for (uint j = 0; j < numsteps; ++j) {
        pp = p + hyperstep[j];

        const uint parity = pp.kappaL5();
        if (parity == 0) continue;

        if (checkProjInWindow(pp, parity - 1)) lvlman.insert(pp);
      }
    }

    lvlman.advance();
  }

  cerr << "Constructed patch of rhombic Penrose tiling with "
       << tilingpoints.size() << " vertices.\n";
}

void RhombicPenrose::projTilingVis(const vec4i& initpoint,
             const vec4i& origin,
             uint maxstep, enum Common::proj_tiling_hint hint,
             Common::vec4ilist& tilingpoints,
             Common::vec4ilist& visiblepoints) {
  using namespace Common;

  vec4ilist tiling_sector;

  vec4i p, pp;
  const uint numsteps = 10;
  const vec4i hyperstep[10] = {vec4i(1,0,0,0),  vec4i(0,1,0,0),
                               vec4i(0,0,1,0),  vec4i(0,0,0,1),
                               vec4i(1,1,1,1),  vec4i(-1,0,0,0),
                               vec4i(0,-1,0,0), vec4i(0,0,-1,0),
                               vec4i(0,0,0,-1), vec4i(-1,-1,-1,-1)};

  if (initpoint.kappaL5() != 0) {
    cerr << "Initial point not of zero-parity.\n";
    return;
  }

  tilingpoints.clear();
  visiblepoints.clear();

  tilingpoints.push_back(initpoint);

  // We need 2 + 1 levels to avoid going "back" (into the wrong direction) when creating the patch.
  TVLManager<vec4ilist, 2 + 1> lvlman(tilingpoints);

  for (uint n = 0; n < maxstep; ++n) {
    for (uint i = lvlman.begin(); i < lvlman.end(); ++i) {
      p = tilingpoints[i];

      for (uint j = 0; j < numsteps; ++j) {
        pp = p + hyperstep[j];

        const uint parity = pp.kappaL5();
        if (parity == 0) continue;

        if (checkProjInWindow(pp, parity - 1)) lvlman.insert(pp);
      }
    }

    lvlman.advance();
  }

  cerr << "Constructed patch of rhombic penrose tiling with "
       << tilingpoints.size() << " vertices.\n";

  switch (hint) {
  case proj_tiling_none:
    extractVisible(origin, false, tilingpoints, visiblepoints);

  case proj_tiling_radialprojection:
    /* Visibility computation is expensive here, so reduce the tiling
     * to a sector before proceeding. */
    extractSector(tilingpoints, tiling_sector);
    tilingpoints.swap(tiling_sector);
    cerr << "Reduced full tiling to a sector containing "
         << tilingpoints.size() << " vertices.\n";

    /* If tiling vertices are going to be used for radial projection, we can
     * proceed using the fastvis routines. */
    extractVisibleFast(origin, tilingpoints, visiblepoints);
    break;

  case proj_tiling_onlysector:
    assert(origin.isZero());
    extractVisible(origin, true, tilingpoints, visiblepoints);
    break;

  default:
    assert(false);
    break;
  }
}

void RhombicPenrose::extractSector(const Common::vec4ilist& input,
             Common::vec4ilist& output) {
  // Identical to the decagonal case:
  Decagonal::extractSector(input, output);
}

void RhombicPenrose::extractVisible(const vec4i& origin, bool sector,
                      const Common::vec4ilist& input, Common::vec4ilist& output) {
  using namespace Common;

  VisList* vlist = new VisList;

  // Assert that we're called with valid input parameters.
  if (sector)
    assert(origin.isZero());

  if (sector)
    vlist->reserve((input.size() - 1) / 5);
  else
    vlist->reserve(input.size() - 1);

  vlist->init();

  if (sector) {
    for (vec4ilist::const_iterator i = input.begin(); i != input.end(); ++i) {
      if (i->isZero()) continue;

      const vec2d physProj(i->paraProjL5());

      if (physProj.inFirstQuadrant() && physProj.inSectorL5()) {
        vlist->insertSorted(*i);
      }
    }
  } else {
    for (vec4ilist::const_iterator i = input.begin(); i != input.end(); ++i) {
      const vec4i shifted(*i - origin);

      if (shifted.isZero()) continue;
      vlist->insertSorted(shifted);
    }
  }

  /* This function is never used for radial projection, so always
   * apply proper/correct visibility computation. */
  vlist->removeInvisibleProper();

  output.clear();
  output.reserve(vlist->size());
  vlist->dump(output);

  delete vlist;
}

void RhombicPenrose::extractVisibleFast(const vec4i& origin,
                const Common::vec4ilist& input, Common::vec4ilist& output) {
  // Identical to the decagonal case:
  Decagonal::extractVisibleFast(origin, input, output);
}

uint RhombicPenrose::estimateGrowth(uint input, bool steps) {
  const double x = double(input);

  if (steps) {
    static const double params[4] = {0.209245, -2.00585, 4.02342, 9.55819};

    return uint(params[0] * (params[1] + sqrt(params[2] + params[3] * x)));
  } else {
    // Parameters were computed using Mathematica's FindFit
    static const double params[2] = {2.00585, 2.38955};

    return uint(params[0] * x + params[1] * x * x);
  }
}

void RhombicPenrose::radialProj(const Common::vec4ilist& input,
             Common::dlist& output, double& meandist) {
  // Identical to the decagonal case:
  Decagonal::radialProj(input, output, meandist);
}

void RhombicPenrose::radialProj(const Common::vec4ilist& input,
             Common::dlist& output, double& meandist, uint window) {
  using namespace Common;

  output.clear();
  output.reserve(input.size());

  dlist angles;
  angles.reserve(input.size());

  const uint parity = window + 1;
  for (vec4ilist::const_iterator i = input.begin(); i != input.end(); ++i) {
    if (i->kappaL5() == parity) {
      const vec2d physProj(i->paraProjL5());
      angles.push_back(physProj.angle());
    }
  }

  sort(angles.begin(), angles.end());
  neighbourDiff(angles, output, meandist);
  normalizeAngDists(output, meandist);
}

// The 'sector' argument is ignored.
void RhombicPenrose::projTilingVis2(const vec4i& initpoint,
          uint maxstep, bool sector,
          Common::vec4ilist& tilingpoints, Common::vec4ilist& visiblepoints) {
  static const vec4i def_origin(0, 0, 0, 0);

  projTilingVis(initpoint, def_origin, maxstep, Common::proj_tiling_radialprojection,
                tilingpoints, visiblepoints);
}

void RhombicPenrose::testWindow(Common::vec2ilist& output,
              uint resolution, uint window) {
  using namespace Common;

  assert(window < 4);

  // Scan the area [-2,2]^2 (outer radius of largest window is approx. 1.62)
  const double step = 4.0 / double(resolution);

  output.clear();
  output.reserve(resolution * resolution);

  for (uint i = 0; i < resolution; ++i) {
    for (uint j = 0; j < resolution; ++j) {
      const vec2d pos(-2.0 + i * step, -2.0 + j * step);

      if (checkProjInSector(pos, window)) {
        output.push_back(vec2i(i, j));
      }
    }
  }
}

