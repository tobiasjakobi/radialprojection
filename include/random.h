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

#ifndef _RANDOM_VIS_H_
#define _RANDOM_VIS_H_

#include "common.h"

namespace RandomVis {

  /* First compute visible points, then randomize the set. */
  void vVisibleRandom(uint steps, double prob, Common::vec2ilist& out);

  /* First randomize the set, then compute visible points. Here of *
   * course the default Z2 visibility test doesn't hold anymore.   */
  void vRandomVisible(uint steps, double prob, Common::vec2ilist& out);

  void radialProjVisRnd(uint steps, double prob, Common::dlist& out);
  void radialProjRndVis(uint steps, double prob, Common::dlist& out);

};

#endif
