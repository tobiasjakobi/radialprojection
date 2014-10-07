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

  // Very primitive integer factorization (should work for small numbers)
  void factorInteger(uint i, vector<uint>& factorization);

  /* Assuming that p = +1 or -1 (mod 8), this finds the tuple *
   * (m,n) that solves the equation algnorm(m,n) = p.         */
  void findTupleZ2(const int p, vec2i& out);

  /* The two conditions that apply to our visibility checks: *
   * p = +1 or -1 (mod 8) (first condition)                  *
   * p = +3 or -3 (mod 8) (second condition)                 */
  bool pCond1Z2(const int p);
  bool pCond2Z2(const int p);

  /* Do prime factorization in Z[Sqrt[2]]. This takes 'in' as an *
   * element of Z[Sqrt[2]] and returns a list of primes that     *
   * divide 'in' (without mulitplicity).                         */
  void factorZ2(const vec2i& in, vector<vec2i>& factorization);
};

namespace ArithVisibility {
  /* Division test (square-free case) for the primes which *
   * satisfy the first condition.                          */
  bool divTest2Free1(const vec2i& in, const int p);

  /* Division test (square-free case) for the primes which *
   * satisfy the second condition.                         */
  bool divTest2Free2(const vec2i& in, const int p);

  // Check for square-free 'visibility' of an element of Z[Sqrt[2]]
  bool visibility2Free(const vec2i& in);

  /* Let x = in / c, an element of Q(Sqrt[2]), then this computes    *
   * the denominator in the Fourier module Z[Sqrt[2]] * (Sqrt[2]/4). */
  vec2i denomZ2Fourier(const vec2i& in, const int in_c);

  /* Compute the intensity of 'denom', an element of Z[Sqrt[2]].       *
   * The input are going to be denominators of elements of Q(Sqrt[2]). */
  double intensityZ2(const vec2i& denom);
};

void vTableZ2(const uint r, Common::vec2ilist& table);

#endif
