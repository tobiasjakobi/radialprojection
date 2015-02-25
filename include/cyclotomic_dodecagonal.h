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

#ifndef _CYCLOTOMIC_DODECAGONAL_H_
#define _CYCLOTOMIC_DODECAGONAL_H_

#include "common.h"
#include "visibility.h"

namespace Dodecagonal {

  /* inflation factor of the dodecagonal tiling */
  const double tau = 2.0 + sqrt(3.0);

  /* dodecagon radii */
  const double innerRadius = tau * 0.5;
  const double outerRadius = sqrt(tau);

  const double innerRadSquared = tau - 0.25;
  const double outerRadSquared = tau;

  const double refCircleRadiusSquared = 3.0 * tau / Common::pi;

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

  /* projTiling only constructs the tiling, starting from the initpoint */
  void projTiling(const vec4i& initpoint, uint maxstep,
                  Common::vec4ilist& tilingpoints);

  /* projTilingVisLocal uses a modified local test to determine         *
   * visibility. The test was derived partly from numerical             *
   * simulations and also from direct computation with the              *
   * geometrical visibility condition plus the model set                *
   * description.                                                       */
  void projTilingVisLocal(const vec4i& initpoint, uint maxstep,
                          Common::vec4ilist& tilingpoints,
                          Common::vec4ilist& visiblepoints);

  /* projTilingVis using a regular ray test to determine visibility:     *
   * This method is the most exact one and should be used for            *
   * reference. It is of course also vastly slower than the              *
   * local test.                                                         */
  void projTilingVis(const vec4i& initpoint,
                     const vec4i& origin,
                     uint maxstep, bool radialproj,
                     Common::vec4ilist& tilingpoints,
                     Common::vec4ilist& visiblepoints);

  /* Compute the squared distance of v and w.                                     *
   * Both v and w have to be in standard coordinatization (_not_ direct-sum one). */
  vec4i sqDist(const vec4i& v, const vec4i& w);

  /* Creates a patch of the tiling via the projection method.            *
   * Edges are computed by connecting vertices with a specific distance. */
  void projTilingPatch(const vec4i& initpoint, uint maxstep,
                       Common::vec4ilist& tilingpoints,
                       Common::edgelist& edges);

  void extractSector(const Common::vec4ilist& input,
                     Common::vec4ilist& output);

  void radialProj(const Common::vec4ilist& input,
                  Common::dlist& output,
                  double& meandist);

  void testWindow(Common::vec2ilist& output, uint resolution);

  /* Select the vertices where the gcd of the direct sum representation *
   * has norm two. These vertices are problematic when computing the    *
   * visibility of a given vertex (based on local information).         */
  struct NormTwoSelector {
    static bool eval(const vec4i& x) {
       const vec2i gcd(Coprime::gcdZ3(x.transL12ToDirect()));
       return (gcd.normZ3() == 2);
    }
  };

  struct VisOp {
    typedef vec4i invectype;
    static const double epsilon;

    static inline double angle(const invectype& a) {
      return a.paraProjL12().angle();
    }

    static inline vec2d toR2(const invectype& a) {
      return a.paraProjL12();
    }

    static bool rayTest(const invectype& a, const invectype& b);
  };

  typedef VisTest::VisibleList<VisOp> VisList;

};

#endif /* _CYCLOTOMIC_DODECAGONAL_H_ */

