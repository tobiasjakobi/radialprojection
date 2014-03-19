#include "cyclotomic_rhombic_penrose.h"

#include "level_manager.h"
#include "cyclotomic_decagonal.h"

#include <algorithm>

const double RhombicPenrose::VisOp::epsilon = 2.0 * numeric_limits<double>::epsilon();

void RhombicPenrose::getInnerOuterSquared(double& inner, double& outer, uint window) {
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

bool RhombicPenrose::checkProjInSector(const vec2d& orthpoint, uint window) {
  using namespace Common;

  const vec2d v(orthpoint.x, abs(orthpoint.y));
  double test;

  const vec2d* const verts = vertices[window];

  for (uint i = 0; i < 3; ++i) {
    test = checkPosition(verts[i], verts[i+1], v);
    if (test < -eps) return false;
    if (test <= eps) {
      cerr << "Warning: Insufficient accuracy in function checkProjInSector.\n";
    }
  }

  return true;
}

bool RhombicPenrose::checkProjInWindow(const vec4i& point, uint window) {
  using namespace Common;

  const vec2d pt(point.orthProjShiftL5());
  const double pt1 = pt.lengthSquared();

  double innerSquared, outerSquared;
  getInnerOuterSquared(innerSquared, outerSquared, window);

  if (innerSquared - pt1 > eps) {
    return true;
  } else {
    if (outerSquared - pt1 < -eps) {
      return false;
    } else {
      return checkProjInSector(pt, window);
    }
  }
}

bool RhombicPenrose::checkScaledProjInWindow(const vec4i& point, uint window) {
  using namespace Common;

  // Note that scaling with -tau (MINUS!) is correct!
  // Usually this doesn't matter since the window is symmetric under
  // transformation with -1, but not here.
  const vec2d pt(point.orthProjShiftL5() * -tau);
  const double pt1 = pt.lengthSquared();

  double innerSquared, outerSquared;
  getInnerOuterSquared(innerSquared, outerSquared, window);

  if (innerSquared - pt1 > eps) {
    return true;
  } else {
    if (outerSquared - pt1 < -eps) {
      return false;
    } else {
      return checkProjInSector(pt, window);
    }
  }
}

void RhombicPenrose::projTilingAll(const vec4i& initpoint, uint maxstep,
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

  TVLManager<vec4i> lvlman(2 + 1, tilingpoints);

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
}

void RhombicPenrose::selectVisible(const Common::vec4ilist& patch,
             Common::vec4ilist& visiblepoints, bool radialproj) {
  using namespace Common;

  VisList* vlist = new VisList;
  vlist->reserve(patch.size());
  vlist->init();

  for (vec4ilist::const_iterator i = patch.begin(); i != patch.end(); ++i) {
    vlist->insertSorted(*i);
  }

  if (radialproj)
    vlist->removeInvisibleFast();
  else
    vlist->removeInvisibleProper();

  visiblepoints.reserve(vlist->size());
  vlist->dump(visiblepoints);

  delete vlist;
  vlist = NULL;
}

void RhombicPenrose::extractSector(const Common::vec4ilist& input,
             Common::vec4ilist& output) {
  // Identical to the decagonal case:
  Decagonal::extractSector(input, output);
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

bool RhombicPenrose::VisOp::rayTest(const invectype& a, const invectype& b) {
  // transform into the Z[tau]*1 + Z[tau]*xi
  // representation (this is a direct sum)
  const invectype pa(a.transL5ToDirect());
  const invectype pb(b.transL5ToDirect());

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

