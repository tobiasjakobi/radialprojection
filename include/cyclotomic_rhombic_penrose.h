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

#ifndef _CYCLOTOMIC_RHOMBIC_PENROSE_H_
#define _CYCLOTOMIC__RHOMBIC_PENROSE_H_

#include "common.h"
#include "visibility.h"

namespace RhombicPenrose {

  /* inflation factor of the rhombic penrose tiling */
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

  /* projTiling only constructs the tiling, starting from the initpoint */
  void projTiling(const vec4i& initpoint, uint maxstep,
                  Common::vec4ilist& tilingpoints);

  void projTilingVis(const vec4i& initpoint,
                     const vec4i& origin,
                     uint maxstep, bool radialproj,
                     Common::vec4ilist& tilingpoints,
                     Common::vec4ilist& visiblepoints);

  /* The rhombic Penrose tiling doesn't allow for a local visibility test. */

  void extractSector(const Common::vec4ilist& input,
                     Common::vec4ilist& output);

  void extractVisible(const vec4i& origin, bool radialproj, bool onlySector,
                      const Common::vec4ilist& input,
                      Common::vec4ilist& output);

  // Compute estimations based on given 'maxstep' or 'size' value.
  uint estimateGrowth(uint input, bool steps);

  void radialProj(const Common::vec4ilist& input,
                  Common::dlist& output,
                  double& meandist);

  /* This alternate version only projects the vertices which come from *
   * a specific window (and therefore have a specific 'parity').       */
  void radialProj(const Common::vec4ilist& input,
                  Common::dlist& output,
                  double& meandist, uint window);

  void testWindow(Common::vec2ilist& output, uint resolution, uint window);

  struct VisOp {
    typedef vec4i invectype;
    static const double epsilon;

    static inline double angle(const invectype& a) {
      return a.paraProjL5().angle();
    }

    static inline vec2d toR2(const invectype& a) {
      return a.paraProjL5();
    }

    static bool rayTest(const invectype& a, const invectype& b);
  };

  typedef VisTest::VisibleList<VisOp> VisList;

};

#endif /* _CYCLOTOMIC_RHOMBIC_PENROSE_H_ */

