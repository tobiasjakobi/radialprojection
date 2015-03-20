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

#ifndef _CYCLOTOMIC_DECAGONAL_H_
#define _CYCLOTOMIC_DECAGONAL_H_

#include "common.h"
#include "visibility.h"

namespace Decagonal {

  /* inflation factor of the decagonal tiling */
  const double tau = 0.5 * (1.0 + sqrt(5.0));

  /* decagon radii */
  const double innerRadius = 0.5 * (tau + 1.0);
  const double outerRadius = tau*tau / sqrt(tau*tau + 1.0); /* = sqrt(1 + 2/sqrt(5)) */

  const double innerRadSquared = 0.25 * (3.0 * tau + 2.0);
  const double outerRadSquared = 1.0 + 2.0/sqrt(5.0);

  const double refCircleRadiusSquared = sqrt((5.0 / 4.0) * (11.0 * tau + 7.0)) / Common::pi;

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

  /* Important: When using verticesAlt as window, it needs to be shifted    *
   * by a small epsilon. Otherwise FP precision issues appear since certain *
   * lattice vector are projected onto the window boundary.                 */

  bool checkProjInSector(const vec2d& orthpoint, bool useAlt);
  bool checkProjInWindow(const vec4i& point, bool useCircle);
  bool checkScaledProjInWindow(const vec4i& point, bool useCircle);

  /* projTiling only constructs the tiling, starting from the initpoint */
  void projTiling(const vec4i& initpoint, uint maxstep,
                  Common::vec4ilist& tilingpoints);

  void projTilingVis(const vec4i& initpoint,
                     const vec4i& origin,
                     uint maxstep, bool radialproj,
                     Common::vec4ilist& tilingpoints,
                     Common::vec4ilist& visiblepoints);

  void projTilingVisLocal(const vec4i& initpoint, uint maxstep,
                          Common::vec4ilist& tilingpoints,
                          Common::vec4ilist& visiblepoints);

  void extractSector(const Common::vec4ilist& input,
                     Common::vec4ilist& output);

  void extractVisible(const vec4i& origin, bool radialproj,
                      const Common::vec4ilist& input,
                      Common::vec4ilist& output);

  // Compute estimations based on given 'maxstep' or 'size' value.
  uint estimateGrowth(uint input, bool steps);

  void radialProj(const Common::vec4ilist& input,
                  Common::dlist& output,
                  double& meandist);

  /* radial projection of a given input tiling, but computed from a non-default origin/radius */
  void radialProj(const Common::vec4ilist& input,
                  const vec4i& origin, double radius,
                  Common::dlist& output, double& meandist);

  void testWindow(Common::vec2ilist& output, uint resolution);

  /* Compute inner and outer radius of tiling. */
  void innerOuterRadius(const Common::vec4ilist& tilingpoints,
                        double& inner, double& outer);

  struct LengthSelector {
    static double length(const vec4i& x) {
      return x.paraProjL5().length();
    }
  };

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

  typedef VisTest::VisibleList<VisOp> VisList;

};

#endif /* _CYCLOTOMIC_DECAGONAL_H_ */

