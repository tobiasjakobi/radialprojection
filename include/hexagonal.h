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

#ifndef _HEXAGONAL_H_
#define _HEXAGONAL_H_

#include "common.h"

#include <algorithm>

namespace Triangular {

  /* The vertices of the triangular tiling are exactly the *
   * elements of the triangular lattice.                   */

  /* The lattice vectors are (1, 0) and (1/2, sqrt(3)/2). */

  const double radiusFactor = sqrt(3.0) * 0.5;

  void tiling(const vec2i& initpoint, uint maxstep,
              Common::vec2ilist& tilingpoints);

  void tilingVisLocal(const vec2i& initpoint, uint maxstep,
                      Common::vec2ilist& tilingpoints,
                      Common::vec2ilist& visiblepoints);

  void extractSector(const Common::vec2ilist& input,
                     Common::vec2ilist& output);

  void radialProj(const Common::vec2ilist& input,
                  Common::dlist& output, double& meandist);

};

/* Hexagonal tiling, also known as honeycomb structure. */
namespace Hexagonal {

  void tiling(const vec2i& initpoint, uint maxstep,
              Common::vec2ilist& tilingpoints);

  void tilingVis(const vec2i& initpoint, uint maxstep,
                 Common::vec2ilist& tilingpoints,
                 Common::vec2ilist& visiblepoints);

  void extractSector(const Common::vec2ilist& input,
                     Common::vec2ilist& output);

  void radialProj(const Common::vec2ilist& input,
                  Common::dlist& output, double& meandist);

};

namespace GenericLattice {

  /* This does not really construct a tiling, but like with the *
   * triangular tiling, it constructs a lattice.                */

  /* The lattice vectors are (1, 0) and the second one *
   * is given by the parameter 'lattice'.              */

  void tiling(const vec2i& initpoint, const vec2d& lattice,
              uint maxstep, Common::vec2ilist& tilingpoints);

  void tilingVisLocal(const vec2i& initpoint, const vec2d& lattice,
                      uint maxstep, Common::vec2ilist& tilingpoints,
                      Common::vec2ilist& visiblepoints);

  void radialProj(const Common::vec2ilist& input,
                  const vec2d& lattice, Common::dlist& output,
                  double& meandist);

};

#endif
