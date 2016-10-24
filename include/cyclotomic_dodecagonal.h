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

#ifndef _CYCLOTOMIC_DODECAGONAL_H_
#define _CYCLOTOMIC_DODECAGONAL_H_

#include "common.h"
#include "visibility.h"

namespace Dodecagonal {

  // Inflation factor of the dodecagonal tiling is the unit of Z[Sqrt[3]].

  /* projTiling only constructs the tiling, starting from the initpoint */
  void projTiling(const vec4i& initpoint, uint maxstep,
                  Common::vec4ilist& tilingpoints);

  /*
   * projTilingVisLocal() uses a modified local test to determine visibility.
   * The test was derived partly from numerical simulations and also from
   * direct computation with the geometrical visibility condition plus the
   * model set description.
   */
  void projTilingVisLocal(const vec4i& initpoint, uint maxstep,
                          bool sector,
                          Common::vec4ilist& tilingpoints,
                          Common::vec4ilist& visiblepoints);

  /*
   * projTilingVis() using a regular ray test to determine visibility:
   * This method is the most exact one and should be used for reference. It
   * is of course also vastly slower than the local test.
   */
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

  void extractVisible(const vec4i& origin, bool radialproj,
                      const Common::vec4ilist& input,
                      Common::vec4ilist& output);

  // Fast computation of (incorrect) visibility, useful for radial projection.
  void extractVisibleFast(const vec4i& origin, const Common::vec4ilist& input,
                          Common::vec4ilist& output);

  // Compute estimations based on given 'maxstep' or 'size' value.
  uint estimateGrowth(uint input, bool steps);

  void radialProj(const Common::vec4ilist& input,
                  Common::dlist& output,
                  double& meandist);

  // Aliases for function pointer compatibility reasons
  void projTilingVisLocal2(const vec4i& initpoint, uint maxstep,
                           Common::vec4ilist& tilingpoints,
                           Common::vec4ilist& visiblepoints);

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

};

#endif /* _CYCLOTOMIC_DODECAGONAL_H_ */

