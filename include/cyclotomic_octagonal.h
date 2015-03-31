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

#ifndef _CYCLOTOMIC_OCTAGONAL_H_
#define _CYCLOTOMIC_OCTAGONAL_H_

#include "common.h"
#include "visibility.h"

namespace Octagonal {

  /* Inflation factor of the octagonal tiling is the silver
   * mean, which is also the unit of Z[Sqrt[2]]. */

  // Octagon radii
  const double innerRadiusSquared = (2.0 * Constants::unitZ2 + 1.0) / 8.0;
  const double outerRadiusSquared = (Constants::unitZ2 + 1.0) / 4.0;

  const double refCircleRadiusSquared = Constants::unitZ2 / Constants::pi;

  bool checkProjInSector(const vec2d& orthpoint);
  bool checkProjInWindow(const vec4i& point, bool useCircle);
  bool checkScaledProjInWindow(const vec4i& point, bool useCircle);

  // Only construct the tiling, starting from the initpoint.
  void projTiling(const vec4i& initpoint, uint maxstep,
                  Common::vec4ilist& tilingpoints);

  /* Constructs the tiling pointset and computes visibility. *
   * An arbitrary origin can be used.                        */
  void projTilingVis(const vec4i& initpoint,
                     const vec4i& origin,
                     uint maxstep, bool radialproj,
                     Common::vec4ilist& tilingpoints,
                     Common::vec4ilist& visiblepoints);

  // Compute visibility with the local test (much faster).
  void projTilingVisLocal(const vec4i& initpoint,
                          uint maxstep, bool onlySector,
                          Common::vec4ilist& tilingpoints,
                          Common::vec4ilist& visiblepoints);

  // Just compute visibility for vertices from 'input' (seen from 'origin').
  void extractVisible(const vec4i& origin, bool radialproj,
                      const Common::vec4ilist& input,
                      Common::vec4ilist& output);

  // Fast computation of (incorrect) visibility, for radial projection.
  void extractVisibleFast(const vec4i& origin, const Common::vec4ilist& input,
                          Common::vec4ilist& output);

  // Compute estimations based on given 'maxstep' or 'size' value.
  uint estimateGrowth(uint input, bool steps);

  void radialProj(const Common::vec4ilist& input,
                  Common::dlist& output,
                  double& meandist, bool onlySector);

  /* Radial projection of a given input tiling, but *
   * computed from a non-default origin/radius.     */
  void radialProj(const Common::vec4ilist& input,
                  const vec4i& origin, double radius,
                  Common::dlist& output, double& meandist);

  // Aliases for function pointer compatibility reasons
  void projTilingVisLocal2(const vec4i& initpoint, uint maxstep,
                           Common::vec4ilist& tilingpoints,
                           Common::vec4ilist& visiblepoints);
  void radialProj2(const Common::vec4ilist& input, Common::dlist& output,
                   double& meandist);

  void testWindow(Common::vec2ilist& output, uint resolution);

  // Compute inner and outer radius of tiling.
  void innerOuterRadius(const Common::vec4ilist& tilingpoints,
                        double& inner, double& outer);

  struct LengthSelector {
    static double length(const vec4i& x) {
      return x.paraProjL8().length();
    }
  };

  struct VisOp {
    typedef Common::vec4ilist list_type;
    static const double epsilon;

    static inline double angle(const vec4i& a) {
      return a.paraProjL8().angle();
    }

    static inline vec2d toR2(const vec4i& a) {
      return a.paraProjL8();
    }

    static bool rayTest(const vec4i& a, const vec4i& b);
  };

  typedef VisTest::VisibleList<VisOp> VisList;

};

#endif /* _CYCLOTOMIC_OCTAGONAL_H_ */

