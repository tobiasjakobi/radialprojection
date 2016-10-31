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

#ifndef _CYCLOTOMIC_RHOMBIC_PENROSE_H_
#define _CYCLOTOMIC__RHOMBIC_PENROSE_H_

#include "common.h"
#include "visibility.h"

namespace RhombicPenrose {

  /* Inflation factor of the decagonal tiling is the golden
   * mean, which is also the unit of Z[tau]. */

  /* projTiling only constructs the tiling, starting from the initpoint */
  void projTiling(const vec4i& initpoint, uint maxstep,
                  Common::vec4ilist& tilingpoints);

  void projTilingVis(const vec4i& initpoint, const vec4i& origin,
                     uint maxstep, bool sector, bool radproj,
                     Common::vec4ilist& tilingpoints,
                     Common::vec4ilist& visiblepoints);

  /*
   * Construction of the tiling vertices followed by fast computation
   * of visible vertices. This method should be used for radial
   * projection only.
   */
  void projTilingVisFast(const vec4i& initpoint, uint maxstep, bool sector,
                         Common::vec4ilist& tilingpoints,
                         Common::vec4ilist& visiblepoints);

  /* The rhombic Penrose tiling doesn't allow for a local visibility test. */

  // Fast computation of (incorrect) visibility, useful for radial projection.
  void extractVisibleFast(const vec4i& origin, const Common::vec4ilist& input,
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

};

#endif /* _CYCLOTOMIC_RHOMBIC_PENROSE_H_ */

