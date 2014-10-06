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

#ifndef _ARITH_VISIBILITY_H_
#define _ARITH_VISIBILITY_H_

#include "common.h"

namespace Coprime {
	
  /* Assuming that p = +1 or -1 (mod 8), this finds the tuple *
   * (m,n) that solves the equation algnorm(m,n) = p.         */
  void findTupleZ2(const int p, vec2i& out);

  void squareZ2(const vec2i& a, const vec2i& b, vec2i& out);
  void cubeZ2(const vec2i& a, const vec2i& b, vec2i& out);

  /* The two conditions that apply to our visibility checks: *
   * p = +1 or -1 (mod 8) (first condition)                  *
   * p = +3 or -3 (mod 8) (second condition)                 */
  bool pCond1Z2(const int p);
  bool pCond2Z2(const int p);
};

#endif
