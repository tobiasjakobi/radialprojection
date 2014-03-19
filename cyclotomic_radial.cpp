#include "cyclotomic_radial.h"

#include <vector>
#include <limits>
#include <sstream>
#include <fstream>
#include <algorithm>

#include "level_manager.h"

vec2d vec4i::shift(0.0, 0.0);

namespace CommonRadial {

  const double eps = numeric_limits<double>::epsilon();

  /* Use the decagon/dodecagon orientation which results from connecting the *
   * ten/twelve roots of unity (false), or the orientation found in the      *
   * book (true), where the right-most edge is aligned with the y-axis.      */
  const bool windowBookOrientation = true;

  // Replace the window by a circle of the same area.
  const bool circularWindow = false;

  double checkPosition(const vec2d& a, const vec2d& b, const vec2d& c);

  bool circularCheck(double radSquared, double xSquared);

  uint modulo(int input, uint mod) {
    const int temp = input % int(mod);
    return (temp < 0) ? (int(mod) + temp) : temp;
  }

  /* Select vertices based on condition specified in 'S'. */
  template <typename T, typename S>
  void selectVertices(const vector<T>& in, vector<T>& out);

  /* Select the vertices which are contained in a ball of radius R. *
   * The radius is set in radiusSq (=R^2)                           */
  struct RadiusSelector {
    static double radiusSq; /* squared radius */

    static bool eval(const vec4i& x) {
       return (x.paraProjL8().lengthSquared() <= radiusSq);
    }
  };

  /* Select a number of origins from the tiling, based on the following conditions:   *
   * Let R = radius * percentage, then only origins are selected, where the ball of   *
   * radius R is still contained in the tiling (radius = radius of the tiling itself. */
  template <typename T, typename S>
  void selectOrigins(const vector<T>& tiling, vector<T>& origins,
                     uint samples, double radius, float percentage);
};

namespace Octogonal {

  const double silverMean = 1.0 + sqrt(2.0);

  /* octagon radii */
  const double innerRadiusSquared = (2.0 * silverMean + 1.0) / 8.0;
  const double outerRadiusSquared = (silverMean + 1.0) / 4.0;

  const double refCircleRadiusSquared = silverMean / CommonRadial::pi;

  const double VisOp::epsilon = 2.0 * numeric_limits<double>::epsilon();

  bool checkProjInSector(const vec2d& orthpoint);
  bool checkProjInWindow(const vec4i& point, bool useCircle);
  bool checkScaledProjInWindow(const vec4i& point, bool useCircle);

  /* projTiling only constructs the tiling, starting from the initpoint */
  void projTiling(const vec4i& initpoint, uint maxstep,
                  CommonRadial::vec4ilist& tilingpoints);

  void projTilingVis(const vec4i& initpoint,
                     const vec4i& origin,
                     uint maxstep, bool radialproj,
                     CommonRadial::vec4ilist& tilingpoints,
                     CommonRadial::vec4ilist& visiblepoints);

  void projTilingVisLocal(const vec4i& initpoint,
                          uint maxstep, bool onlySector,
                          CommonRadial::vec4ilist& tilingpoints,
                          CommonRadial::vec4ilist& visiblepoints);

  void radialProj(const CommonRadial::vec4ilist& input,
                  CommonRadial::dlist& output,
                  double& meandist, bool onlySector);

  /* radial projection of a given input tiling, but computed from a non-default origin/radius */
  void radialProj(const CommonRadial::vec4ilist& input,
                  const vec4i& origin, double radius,
                  CommonRadial::dlist& output, double& meandist);

  void testWindow(CommonRadial::vec2ilist& output, uint resolution);

  /* Compute inner and outer radius of tiling. */
  void innerOuterRadius(const CommonRadial::vec4ilist& tilingpoints,
                        double& inner, double& outer);

};

namespace Decagonal {

  const double tau = 0.5 * (1.0 + sqrt(5.0));

  /* decagon radii */
  const double innerRadius = 0.5 * (tau + 1.0);
  const double outerRadius = tau*tau / sqrt(tau*tau + 1.0); /* = sqrt(1 + 2/sqrt(5)) */

  const double innerRadSquared = 0.25 * (3.0 * tau + 2.0);
  const double outerRadSquared = 1.0 + 2.0/sqrt(5.0);

  const double refCircleRadiusSquared = sqrt((5.0 / 4.0) * (11.0 * tau + 7.0)) / CommonRadial::pi;

  /* Regular decagon with edge length sqrt((tau + 2)/5). *
   * The upper-most edge is aligned with the x-axis.     */
  const vec2d vertices[10] = {
    vec2d(outerRadius, 0.0),                          vec2d(sqrt((11.0*tau + 7.0)/5.0)*0.5, tau*0.5),
    vec2d(sqrt((tau + 2.0)/5.0)*0.5, innerRadius),    vec2d(-sqrt((tau + 2.0)/5.0)*0.5, innerRadius),
    vec2d(-sqrt((11.0*tau + 7.0)/5.0)*0.5, tau*0.5),  vec2d(-outerRadius, 0.0),
    vec2d(-sqrt((11.0*tau + 7.0)/5.0)*0.5, -tau*0.5), vec2d(-sqrt((tau + 2.0)/5.0)*0.5, -innerRadius),
    vec2d(sqrt((tau + 2.0)/5.0)*0.5, -innerRadius),   vec2d(sqrt((11.0*tau + 7.0)/5.0)*0.5, -tau*0.5)};

  /* This is essentially the other decagon rotated by 18 degrees. *
   * Results in alignment of the right-most edge with the y-axis. */
  const vec2d verticesAlt[10] = {
    vec2d(innerRadius, -0.5 * sqrt((tau + 1.0)/(tau + 2.0))),    vec2d(innerRadius, 0.5 * sqrt((tau + 1.0)/(tau + 2.0))),
    vec2d(0.5 * tau, 0.5 * sqrt((8.0 * tau + 5.0)/(tau + 2.0))), vec2d(0.0, outerRadius), vec2d(-0.5 * tau,
    0.5 * sqrt((8.0 * tau + 5.0)/(tau + 2.0))),                  vec2d(-innerRadius, 0.5 * sqrt((tau + 1.0)/(tau + 2.0))),
    vec2d(-innerRadius, -0.5 * sqrt((tau + 1.0)/(tau + 2.0))),   vec2d(-0.5 * tau, -0.5 * sqrt((8.0 * tau + 5.0)/(tau + 2.0))),
    vec2d(0.0, -outerRadius),                                    vec2d(0.5 * tau, -0.5 * sqrt((8.0 * tau + 5.0)/(tau + 2.0)))
  };

  /* Important: When using decagonVerticesAlt as window, it needs to be shifted *
   * by a small epsilon. Otherwise FP precision issues appear since certain     *
   * lattice vector are projected onto the window boundary.                     */

  bool checkProjInSector(const vec2d& orthpoint, bool useAlt);
  bool checkProjInWindow(const vec4i& point, bool useCircle);
  bool checkScaledProjInWindow(const vec4i& point, bool useCircle);

  void projTilingVisLocal(const vec4i& initpoint, uint maxstep,
                          CommonRadial::vec4ilist& tilingpoints,
                          CommonRadial::vec4ilist& visiblepoints);

  // TODO: add ray visibility test

  void extractSector(const CommonRadial::vec4ilist& input,
                     CommonRadial::vec4ilist& output);

  void radialProj(const CommonRadial::vec4ilist& input,
                  CommonRadial::dlist& output,
                  double& meandist);

  void testWindow(CommonRadial::vec2ilist& output, uint resolution);

};

namespace Dodecagonal {

  const double tau = 2.0 + sqrt(3.0);

  /* dodecagon radii */
  const double innerRadius = tau * 0.5;
  const double outerRadius = sqrt(tau);

  const double innerRadSquared = tau - 0.25;
  const double outerRadSquared = tau;

  const double refCircleRadiusSquared = 3.0 * tau / CommonRadial::pi;

  const double VisOp::epsilon = 2.0 * numeric_limits<double>::epsilon();

  /* Regular dodecagon (12 sides) with edge length one. *
   * Orientation is the one resulting from connecting   *
   * the twelve roots of unity (plus scaling).          */
  const vec2d vertices[12] = {
    vec2d(outerRadius, 0.0), /* v1 */
    vec2d(0.5 * sqrt(3.0) * outerRadius, 0.5 * outerRadius), /* v2 */
    vec2d(0.5 * outerRadius, 0.5 * sqrt(3.0) * outerRadius), /* v3 */
    vec2d(0.0, outerRadius), /* v4 */
    vec2d(-0.5 * outerRadius, 0.5 * sqrt(3.0) * outerRadius), /* v5 */
    vec2d(-0.5 * sqrt(3.0) * outerRadius, 0.5 * outerRadius), /* v6 */
    vec2d(-outerRadius, 0.0), /* v7 = -v1 */
    vec2d(-0.5 * sqrt(3.0) * outerRadius, -0.5 * outerRadius), /* v8 = -v2 */
    vec2d(-0.5 * outerRadius, -0.5 * sqrt(3.0) * outerRadius), /* v9 = -v3 */
    vec2d(0.0, -outerRadius), /* v10 = -v4 */
    vec2d(0.5 * outerRadius, -0.5 * sqrt(3.0) * outerRadius), /* v11 = -v5 */
    vec2d(0.5 * sqrt(3.0) * outerRadius, -0.5 * outerRadius)  /* v12 = -v6 */
  };

  const vec2d verticesAlt[12] = {
    vec2d(0.5 * tau, -0.5), vec2d(0.5 * tau, 0.5), /* v12, v1 */
    vec2d(0.5 * (tau - 1.0), 0.5 * (tau - 1.0)), vec2d(0.5, 0.5 * tau), /* v2, v3 */
    vec2d(-0.5, 0.5 * tau), vec2d(-0.5 * (tau - 1.0), 0.5 * (tau - 1.0)), /* v4, v5 */
    vec2d(-0.5 * tau, 0.5), vec2d(-0.5 * tau, -0.5), /* v6, v7 */
    vec2d(-0.5 * (tau - 1.0), -0.5 * (tau - 1.0)), vec2d(-0.5, -0.5 * tau), /* v8, v9 */
    vec2d(0.5, -0.5 * tau), vec2d(0.5 * (tau - 1.0), -0.5 * (tau - 1.0)) /* v10, v11 */
  };

  bool checkProjInSector(const vec2d& orthpoint, bool useAlt);
  bool checkProjInWindow(const vec4i& point, bool useCircle);
  bool checkScaledProjInWindow(const vec4i& point, bool gcdNormTwo, bool useCircle);

  /* projTilingVisLocal uses a modified local test to determine         *
   * visibility. The test was derived partly from numerical             *
   * simulations and also from direct computation with the              *
   * geometrical visibility condition plus the model set                *
   * description.                                                       */
  void projTilingVisLocal(const vec4i& initpoint, uint maxstep,
                          CommonRadial::vec4ilist& tilingpoints,
                          CommonRadial::vec4ilist& visiblepoints);

  /* projTilingVis using a regular ray test to determine visibility:     *
   * This method is the most exact one and should be used for            *
   * reference. It is of course also vastly slower than the              *
   * local test.                                                         */
  void projTilingVis(const vec4i& initpoint, uint maxstep,
                     bool radialproj,
                     CommonRadial::vec4ilist& tilingpoints,
                     CommonRadial::vec4ilist& visiblepoints);

  /* Compute the squared distance of v and w.                                     *
   * Both v and w have to be in standard coordinatization (_not_ direct-sum one). */
  vec4i sqDist(const vec4i& v, const vec4i& w);

  /* Creates a patch of the tiling via the projection method.            *
   * Edges are computed by connecting vertices with a specific distance. */
  void projTilingPatch(const vec4i& initpoint, uint maxstep,
                       CommonRadial::vec4ilist& tilingpoints,
                       CommonRadial::edgelist& edges);

  void extractSector(const CommonRadial::vec4ilist& input,
                     CommonRadial::vec4ilist& output);

  void radialProj(const CommonRadial::vec4ilist& input,
                  CommonRadial::dlist& output,
                  double& meandist);

  void testWindow(CommonRadial::vec2ilist& output, uint resolution);

  /* Select the vertices where the gcd of the direct sum representation *
   * has norm two. These vertices are problematic when computing the    *
   * visibility of a given vertex (based on local information).         */
  struct NormTwoSelector {
    static bool eval(const vec4i& x) {
       const vec2i gcd(RadialCoprime::gcdZ3(x.transL12ToDirect()));
       return (gcd.normZ3() == 2);
    }
  };

};

namespace RhombicPenrose {

  const double VisOp::epsilon = 2.0 * numeric_limits<double>::epsilon();

  const double tau = 0.5 * (1.0 + sqrt(5.0));

  /* The rhombic penrose case uses four windows W1, W2, W3, W4: *
   * Let P be the convex hull of {1, xi, xi^2, xi^3, xi^4} with *
   * xi = exp(2*pi*i/5), then let                               *
   * W1 = P, W4 = -P, W3 = tau*P, W2 = -tau*P.                  */

  /* pentagon radii */
  const double innerRadius[2] = {0.5 * tau, 0.25 * (tau + 1.0)};
  const double outerRadius[2] = {1.0, tau}; 

  const double innerRadSquared[2] = {0.25 * (tau + 1.0), (7.0 + 3.0*sqrt(5.0)) / 8.0};
  const double outerRadSquared[2] = {1.0, tau*tau};

  const vec2d vertices[4][5] = {
    {vec2d(1.0, 0.0), vec2d(0.5 * (tau - 1.0), 0.5 * sqrt(tau + 2.0)),
     vec2d(-0.5 * tau, 0.5 * sqrt(3.0 - tau)), vec2d(-0.5 * tau, -0.5 * sqrt(3.0 - tau)),
     vec2d(0.5 * (tau - 1.0), -0.5 * sqrt(tau + 2.0))},
    {vec2d(0.5 * (tau + 1.0), -0.5 * sqrt(tau + 2.0)),
     vec2d(0.5 * (tau + 1.0), 0.5 * sqrt(tau + 2.0)), vec2d(-0.5, 0.5 * sqrt(4.0*tau + 3.0)), vec2d(-tau, 0.0), vec2d(-0.5, -0.5 * sqrt(4.0*tau + 3.0))},
    {vec2d(tau, 0.0), vec2d(0.5, 0.5 * sqrt(4.0*tau + 3.0)),
     vec2d(-0.5 * (tau + 1.0), 0.5 * sqrt(tau + 2.0)),
     vec2d(-0.5 * (tau + 1.0), -0.5 * sqrt(tau + 2.0)), vec2d(0.5, -0.5 * sqrt(4.0*tau + 3.0))},
    {vec2d(0.5 * tau, -0.5 * sqrt(3.0 - tau)), vec2d(0.5 * tau, 0.5 * sqrt(3.0 - tau)),
     vec2d(-0.5 * (tau - 1.0), 0.5 * sqrt(tau + 2.0)), vec2d(-1.0, 0.0),
     vec2d(-0.5 * (tau - 1.0), -0.5 * sqrt(tau + 2.0))}
  };

  void getInnerOuterSquared(double& inner, double& outer, uint window);

  bool checkProjInSector(const vec2d& orthpoint, uint window);
  bool checkProjInWindow(const vec4i& point, uint window);
  bool checkScaledProjInWindow(const vec4i& point, uint window);

  void projTilingAll(const vec4i& initpoint, uint maxstep,
                     CommonRadial::vec4ilist& tilingpoints);
  void selectVisible(const CommonRadial::vec4ilist& patch,
                     CommonRadial::vec4ilist& visiblepoints, bool radialproj);

  void extractSector(const CommonRadial::vec4ilist& input,
                     CommonRadial::vec4ilist& output);

  /* The alt(ernate) version of radialProj only projects *
   * the vertices which come from a specific window (and *
   * therefore have a specific "parity".                 */
  void radialProj(const CommonRadial::vec4ilist& input,
                  CommonRadial::dlist& output,
                  double& meandist);
  void radialProj(const CommonRadial::vec4ilist& input,
                  CommonRadial::dlist& output,
                  double& meandist, uint window);

  void testWindow(CommonRadial::vec2ilist& output, uint resolution, uint window);

};

ostream& operator<<(ostream &os, const CommonRadial::vec4ilist& list)
{
  using namespace CommonRadial;

  if (!list.empty()) {
    vec4ilist::const_iterator i = list.begin();

    os << '{' << *i;
    ++i;

    while (i != list.end()) {
      os << ',' << *i;
      ++i;
    }

    os << '}';
  } else {
    os << "{}";
  }

  return os;
}

ostream& operator<<(ostream &os, const CommonRadial::vec2ilist& list)
{
  using namespace CommonRadial;

  if (!list.empty()) {
    vec2ilist::const_iterator i = list.begin();

    os << '{' << *i;
    ++i;

    while (i != list.end()) {
      os << ',' << *i;
      ++i;
    }

    os << '}';
  } else {
    os << "{}";
  }

  return os;
}

double CommonRadial::RadiusSelector::radiusSq = 0.0;

double CommonRadial::checkPosition(const vec2d& a, const vec2d& b, const vec2d& c) {
  const double pos = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
  return pos;
}

bool CommonRadial::circularCheck(double radSquared, double xSquared) {
  if (radSquared - xSquared > eps) {
    return true;
  } else {
    if (radSquared - xSquared < -eps) return false;
  }

  cerr << "Warning: Insufficient accuracy in function circularCheck.\n";
  return false;
}

template <typename T, typename S>
void CommonRadial::selectVertices(const vector<T>& in, vector<T>& out) {
  for (typename vector<T>::const_iterator i = in.begin(); i != in.end(); ++i) {
    if (S::eval(*i)) out.push_back(*i);
  }
}

template <typename T, typename S>
void CommonRadial::selectOrigins(const vector<T>& tiling, vector<T>& origins,
                     uint samples, double radius, float percentage) {
  // TODO: implement

}

bool Octogonal::checkProjInSector(const vec2d& orthpoint) {
  using namespace CommonRadial;

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

bool Octogonal::checkProjInWindow(const vec4i& point, bool useCircle) {
  using namespace CommonRadial;

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

bool Octogonal::checkScaledProjInWindow(const vec4i& point, bool useCircle) {
  using namespace CommonRadial;

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

void Octogonal::projTiling(const vec4i& initpoint, uint maxstep,
                 CommonRadial::vec4ilist& tilingpoints) {
  using namespace CommonRadial;

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

void Octogonal::projTilingVis(const vec4i& initpoint,
                 const vec4i& origin,
                 uint maxstep, bool radialproj,
                 CommonRadial::vec4ilist& tilingpoints,
                 CommonRadial::vec4ilist& visiblepoints) {
  using namespace CommonRadial;

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

  // We're not removing vertices in this case, so allocate the full amount.
  VisList* vlist = new VisList;
  vlist->reserve(tilingpoints.size() - 1);

  vlist->init();

  for (vec4ilist::const_iterator i = tilingpoints.begin(); i != tilingpoints.end(); ++i) {
    const vec4i shifted(*i - origin);

    if (shifted.isZero()) continue;
    vlist->insertSorted(shifted);
  }

  if (radialproj)
    vlist->removeInvisibleFast();
  else
    vlist->removeInvisibleProper();

  visiblepoints.clear();
  visiblepoints.reserve(vlist->size());
  vlist->dump(visiblepoints);

  delete vlist;
  vlist = NULL;
}

// See projTilingVis for more comments on the code.
void Octogonal::projTilingVisLocal(const vec4i& initpoint,
                 uint maxstep, bool onlySector,
                 CommonRadial::vec4ilist& tilingpoints,
                 CommonRadial::vec4ilist& visiblepoints) {
  using namespace CommonRadial;
  using namespace RadialCoprime;

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

void Octogonal::radialProj(const CommonRadial::vec4ilist& input,
                          CommonRadial::dlist& output,
                          double& meandist, bool onlySector) {
  using namespace CommonRadial;

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

void Octogonal::radialProj(const CommonRadial::vec4ilist& input,
                          const vec4i& origin, double radius,
                          CommonRadial::dlist& output, double& meandist) {
  using namespace CommonRadial;

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

void Octogonal::testWindow(CommonRadial::vec2ilist& output,
                          uint resolution) {
  using namespace CommonRadial;

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

void Octogonal::innerOuterRadius(const CommonRadial::vec4ilist& tilingpoints,
                          double& inner, double& outer) {
  using namespace CommonRadial;

  double out = 0.0;

  for (vec4ilist::const_iterator i = tilingpoints.begin(); i != tilingpoints.end(); ++i) {
    const double lSq = i->paraProjL8().lengthSquared();

    if (lSq > out) out = lSq;
  }

  outer = sqrt(out);
  inner = cos(pi / 8.0) * outer;

  cerr << "debug: outer = " << outer << endl;
  cerr << "debug: inner = " << inner << endl;
}

bool Decagonal::checkProjInSector(const vec2d& orthpoint, bool useAlt) {
  using namespace CommonRadial;

  const vec2d v(abs(orthpoint.x), abs(orthpoint.y));
  double test;

  const vec2d* const vertices = (useAlt ? verticesAlt : vertices);

  for (uint i = 0; i < 3; ++i) {
    test = checkPosition(vertices[i], vertices[i+1], v);
    if (test < -eps) return false;
    if (test <= eps) {
      cerr << "Warning: Insufficient accuracy in function checkProjInSector.\n";
    }
  }

  return true;
}

bool Decagonal::checkProjInWindow(const vec4i& point, bool useCircle) {
  using namespace CommonRadial;

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
  using namespace CommonRadial;

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

void Decagonal::projTilingVisLocal(const vec4i& initpoint, uint maxstep,
                     CommonRadial::vec4ilist& tilingpoints,
                     CommonRadial::vec4ilist& visiblepoints) {
  using namespace CommonRadial;
  using namespace RadialCoprime;

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

void Decagonal::extractSector(const CommonRadial::vec4ilist& input,
                 CommonRadial::vec4ilist& output) {
  using namespace CommonRadial;

  output.clear();
  output.reserve(input.size() / 4);

  for (vec4ilist::const_iterator i = input.begin(); i != input.end(); ++i) {
    const vec2d physProj(i->paraProjL5());

    if (physProj.inFirstQuadrant() && physProj.inSectorL5()) {
      output.push_back(*i);
    }
  }
}

void Decagonal::radialProj(const CommonRadial::vec4ilist& input,
                          CommonRadial::dlist& output,
                          double& meandist) {
  using namespace CommonRadial;

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

void Decagonal::testWindow(CommonRadial::vec2ilist& output, uint resolution) {
  using namespace CommonRadial;

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

bool Dodecagonal::checkProjInSector(const vec2d& orthpoint, bool useAlt) {
  using namespace CommonRadial;

  const vec2d v(orthpoint.reduceIntoSectorL12());
  double test;

  const vec2d* const vertices = (useAlt ? verticesAlt : vertices);

  for (uint i = 0; i < 2; ++i) {
    test = checkPosition(vertices[i], vertices[i+1], v);
    if (test < -eps) return false;
    if (test <= eps) {
      cerr << "Warning: Insufficient accuracy in function checkProjInSector.\n";
    }
  }

  return true;
}

bool Dodecagonal::checkProjInWindow(const vec4i& point, bool useCircle) {
  using namespace CommonRadial;

  const vec2d pt(point.orthProjShiftL12());
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

bool Dodecagonal::checkScaledProjInWindow(const vec4i& point,
                    bool gcdNormTwo, bool useCircle) {
  using namespace CommonRadial;

  /* Empirical tests with large patches indicate that the test  *
   * against a rescaled window depends on the gcd-norm of the   *
   * vertex. The scaling factors look a bit weird, but become   *
   * nicer if we don't consider the dodecagon with outer radius *
   * sqrt(tau) but the normalized dodecagon (outer rad = 1).    */
  const double scaler = gcdNormTwo ? sqrt(tau * 0.5) :
                                     sqrt(tau * 2.0);

  /* Interesting observation: Multiply the two scaling factors and *
   * we recover tauDode again. Hmm, nice!                          */

  const vec2d pt(point.orthProjShiftL12(scaler, gcdNormTwo));
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

void Dodecagonal::projTilingVisLocal(const vec4i& initpoint, uint maxstep,
                     CommonRadial::vec4ilist& tilingpoints,
                     CommonRadial::vec4ilist& visiblepoints) {
  using namespace CommonRadial;
  using namespace RadialCoprime;

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

  TVLManager<vec4i> lvlman(2 + 1, tilingpoints);

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
         * in the octogonal (AB) and decagonal (TTT) case.                  */
        const int gcdnorm = RadialCoprime::gcdZ3(pp.transL12ToDirect()).normZ3();
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

void Dodecagonal::projTilingVis(const vec4i& initpoint, uint maxstep,
                     bool radialproj,
                     CommonRadial::vec4ilist& tilingpoints,
                     CommonRadial::vec4ilist& visiblepoints) {
  using namespace CommonRadial;

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

  cerr << "Constructed patch of dodecagonal tiling with "
       << tilingpoints.size() << " vertices.\n";

  // TODO: Due to the shift into a non-singular position, there is no 12-star at the
  //       origin of the tiling. Can one still reduce the tiling to a, say 1/6 sector?
  VisList* vlist = new VisList;
  vlist->reserve(double(tilingpoints.size() - 1) * (radialproj ? 0.17 : 1.0));

  vlist->init();

  if (radialproj) {
    cerr << "info: processing dodecagonal tiling in RP mode (reduction to sector).\n";

    for (vec4ilist::const_iterator i = tilingpoints.begin(); i != tilingpoints.end(); ++i) {
      if (i->isZero()) continue; // zero = reference point (not visible)

      const vec2d physProj(i->paraProjL12());

      if (physProj.inFirstQuadrant() && physProj.inSectorL12()) {
        vlist->insertSorted(*i);
      }
    }

    vlist->removeInvisibleFast();
  } else {
    for (vec4ilist::const_iterator i = tilingpoints.begin(); i != tilingpoints.end(); ++i) {
      if (i->isZero()) continue;
      vlist->insertSorted(*i);
    }

    vlist->removeInvisibleProper();
  }

  visiblepoints.clear();
  visiblepoints.reserve(vlist->size());
  vlist->dump(visiblepoints);

  delete vlist;
  vlist = NULL;
}

vec4i Dodecagonal::sqDist(const vec4i& v, const vec4i& w) {
  const vec4i z(v - w);

  return z.multL12(z.conjL12());
}

void Dodecagonal::projTilingPatch(const vec4i& initpoint, uint maxstep,
                     CommonRadial::vec4ilist& tilingpoints,
                     CommonRadial::edgelist& edges) {
  using namespace CommonRadial;

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

void Dodecagonal::extractSector(const CommonRadial::vec4ilist& input,
                     CommonRadial::vec4ilist& output) {

  using namespace CommonRadial;

  output.clear();
  output.reserve(input.size() / 6);

  for (vec4ilist::const_iterator i = input.begin(); i != input.end(); ++i) {
    const vec2d physProj(i->paraProjL12());

    if (physProj.inFirstQuadrant() && physProj.inSectorL12()) {
      output.push_back(*i);
    }
  }
}

void Dodecagonal::radialProj(const CommonRadial::vec4ilist& input,
                  CommonRadial::dlist& output,
                  double& meandist) {
  using namespace CommonRadial;

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

void Dodecagonal::testWindow(CommonRadial::vec2ilist& output, uint resolution) {
  using namespace CommonRadial;

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
  using namespace CommonRadial;

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
  using namespace CommonRadial;

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
  using namespace CommonRadial;

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
             CommonRadial::vec4ilist& tilingpoints) {
  using namespace CommonRadial;

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

void RhombicPenrose::selectVisible(const CommonRadial::vec4ilist& patch,
             CommonRadial::vec4ilist& visiblepoints, bool radialproj) {
  using namespace CommonRadial;

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

void RhombicPenrose::extractSector(const CommonRadial::vec4ilist& input,
             CommonRadial::vec4ilist& output) {
  // Identical to the decagonal case:
  Decagonal::extractSector(input, output);
}

void RhombicPenrose::radialProj(const CommonRadial::vec4ilist& input,
             CommonRadial::dlist& output, double& meandist) {
  // Identical to the decagonal case:
  Decagonal::radialProj(input, output, meandist);
}

void RhombicPenrose::radialProj(const CommonRadial::vec4ilist& input,
             CommonRadial::dlist& output, double& meandist, uint window) {
  using namespace CommonRadial;

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

void RhombicPenrose::testWindow(CommonRadial::vec2ilist& output,
              uint resolution, uint window) {
  using namespace CommonRadial;

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

/* multimachine routine that is executed on master machine */
int main_mmachine_master(int argc, char* argv[]) {
  using namespace CommonRadial;

  stringstream parser;
  vec4ilist fulltiling, tiling;
  double inner, outer;
  uint steps, samples;
  uint percentage;

  const vec4i init(0, 0, 0, 0);

  if (argc != 4) return 1;

  parser.str(argv[1]);
  parser.clear();
  parser >> steps;

  parser.str(argv[2]);
  parser.clear();
  parser >> samples;

  parser.str(argv[3]);
  parser.clear();
  parser >> percentage;



  Octogonal::projTiling(init, steps, fulltiling);
  Octogonal::innerOuterRadius(fulltiling, inner, outer);

  RadiusSelector::radiusSq = inner * inner;

  selectVertices<vec4i, RadiusSelector>(fulltiling, tiling);

  // TODO: select points




  writeRawConsole(tiling);
}

/* multimachine routine that is executed on the slave machines */
int main_mmachine_slave(int argc, char* argv[]) {
  using namespace CommonRadial;

  stringstream parser;
  string token;
  unsigned long long temp;
  unsigned idx = 0;

  vec4i origin;
  vec4ilist tiling;
  CommonRadial::dlist output;
  double radius, mdist;

  if (argc != 3) return 1;

  parser.str(argv[1]);
  parser.clear();
  parser >> temp;
  radius = *(reinterpret_cast<const double*>(&temp));

  parser.str(argv[2]);
  parser.clear();
  while (idx < 4) {
    std::getline(parser, token, ',');

    stringstream ss(token);
    ss >> origin[idx];
    ++idx;
  }

  readRawConsole(tiling);
  Octogonal::radialProj(tiling, origin, radius, output, mdist);

  writeRawConsole(output);
}

void reinterpret_double(double d) {
  unsigned long long i;

  i = *(reinterpret_cast<unsigned long long*>(&d));

  cout << i << endl;
}

int main(int argc, char* argv[]) {
  const vec4i init(0, 0, 0, 0);
  uint steps = 40;
  bool sector = false;
  uint mode = 0;
  CommonRadial::vec4ilist tiling, visible;

  vec4i origin;
  bool use_default_origin = true;

  return main_mmachine_slave(argc, argv);
  //return main_mmachine_master(argc, argv);

  /* create dodecagonal tiling with edges
  {
    vec4i::shift.set(1.0e-4, 1.0e-4);
    CommonRadial::edgelist edges;
    DodecagonalRadial::projTilingPatch(init, 15, tiling, edges);
    cout << '{' << tiling << ',' << edges << '}' << endl;
    return 0;
  }
  */

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

  // Parse (non-default) origin. Currently only the octogonal code uses this!
  if (argc >=8) {
    for (uint k = 0; k < 4; ++k) {
      stringstream ss(argv[4 + k]);
      ss >> origin[k];
    }

    if (origin != vec4i(0, 0, 0, 0)) use_default_origin = false;
  }

  // The pentagon windows in the rhombic Penrose case
  // always need to be shifted into a generic position:
  if (mode == 6 || mode == 7) {
    vec4i::shift.set(1.0e-4, 1.0e-4);
  } else {
    // In the pure cyclotomic case the window can be replaced by a circle:
    if (CommonRadial::circularWindow) {
      cerr << "Using (simplified) circular window with same area.\n";
    } else {
      // When using the decagon/dodecagon window from the book, shift it slightly.
      // A singular patch can be constructed by omitting the shift.
      if (CommonRadial::windowBookOrientation && (mode >= 2 && mode <= 5)) {
        cerr << "Using decagon/dodecagon window orientation from the book, applying slight shift to it.\n";
        vec4i::shift.set(1.0e-4, 1.0e-4);
      }
    }
  }

  /* mode 0 and 1 : octogonal tiling (L8 lattice)           *
   * mode 2 and 3 : decagonal tiling (L5 lattice)           *
   * mode 4 and 5 : dodecagonal tiling (L12 lattice)        *
   * mode 6 and 7 : rhombic Penrose tiling (L5 + 4 windows) *
   * Even modes create vertex data in the Z4 lattice space. *
   * Uneven modes are radial projection modes and           *
   * output double-precision floating point data.           */
  switch (mode) {
    case 0: if (use_default_origin)
              Octogonal::projTilingVisLocal(init, steps, sector, tiling, visible);
            else
              Octogonal::projTilingVis(init, origin, steps, false, tiling, visible); // onlySector is ignored
            break;

    case 1: if (use_default_origin)
              Octogonal::projTilingVisLocal(init, steps, sector, tiling, visible);
            else
              Octogonal::projTilingVis(init, origin, steps, true, tiling, visible); // onlySector is ignored
            {
              CommonRadial::dlist output;
              double mean;

              Octogonal::radialProj(visible, output, mean, sector);

              cerr << "mean distance " << mean
                   << " during radial projection of " << (output.size() + 1)
                   << " vertices.\n";

              CommonRadial::writeRawConsole(output);
            }
            break;

    case 2: Decagonal::projTilingVisLocal(init, steps, tiling, visible);
            if (sector) {
              CommonRadial::vec4ilist visSector;
              Decagonal::extractSector(visible, visSector);
              visible = visSector;
              cerr << "Reduced visible tiling to a sector containing "
                   << visible.size() << " vertices.\n";
            }
            break;

    case 3: Decagonal::projTilingVisLocal(init, steps, tiling, visible);
            {
              CommonRadial::vec4ilist visSector;
              CommonRadial::dlist output;
              double mean;

              Decagonal::extractSector(visible, visSector);
              Decagonal::radialProj(visSector, output, mean);

              cerr << "mean distance " << mean
                   << " during radial projection of " << visSector.size()
                   << " vertices.\n";

              CommonRadial::writeRawConsole(output);
            }
            break;

    case 4: {
              bool robustVisibilityTest = true;
              if (argc >= 5) {
                stringstream ss(argv[4]);
                string vismode;
                ss >> vismode;

                if (vismode == "local") robustVisibilityTest = false;
              }
              
              if (robustVisibilityTest)
                Dodecagonal::projTilingVis(init, steps, false, tiling, visible);
              else
                Dodecagonal::projTilingVisLocal(init, steps, tiling, visible);
            }
            if (sector) {
              CommonRadial::vec4ilist visSector;
              Dodecagonal::extractSector(visible, visSector);
              visible = visSector;
              cerr << "Reduced visible tiling to a sector containing "
                   << visible.size() << " vertices.\n";
            }
            break;

    case 5: Dodecagonal::projTilingVis(init, steps, true, tiling, visible);
            {
              CommonRadial::dlist output;
              double mean;

              Dodecagonal::radialProj(visible, output, mean);

              cerr << "mean distance " << mean
                   << " during radial projection of " << visible.size()
                   << " vertices.\n";

              CommonRadial::writeRawConsole(output);
            }
            break;

    case 6: RhombicPenrose::projTilingAll(init, steps, tiling);
            if (sector) {
              CommonRadial::vec4ilist tilingSector;
              RhombicPenrose::extractSector(tiling, tilingSector);
              tiling.swap(tilingSector);
              cerr << "Reduced tiling to a sector containing "
                   << tiling.size() << " vertices.\n";
            }
            RhombicPenrose::selectVisible(tiling, visible, false);
            break;

    case 7: RhombicPenrose::projTilingAll(init, steps, tiling);
            {
              CommonRadial::vec4ilist tilingSector;
              CommonRadial::dlist output;
              double mean;

              RhombicPenrose::extractSector(tiling, tilingSector);
              RhombicPenrose::selectVisible(tilingSector, visible, true);

              RhombicPenrose::radialProj(visible, output, mean);

              cerr << "mean distance " << mean
                   << " during radial projection of " << (output.size() + 1)
                   << " vertices.\n";

              CommonRadial::writeRawConsole(output);
            }
            break;

    default: cerr << "error: unsupported mode selected!\n";
             return 0;
  }

  // Don't output the tiling in radial projection mode:
  // We're only interested in the angular data, which is written
  // to the console in raw mode in this case.
  if (mode % 2 == 0) {
    cerr << "Size of visible point data is around "
         << uint(double(visible.size() * sizeof(vec4i)) / 1024.0)
         << " kilobytes.\n";

    cout << visible;
  }

  return 0;
}

