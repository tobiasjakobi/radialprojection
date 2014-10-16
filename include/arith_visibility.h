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

  /* Assuming that p = +1 (mod 4), this finds the tuple (m,n) *
   * that solves the equation algnorm(m,n) = p.               */
  void findTupleGI(const int p, vec2i& out);

  /* The two conditions that apply to our visibility checks (for Z[Sqrt[2]]): *
   * p = +1 or -1 (mod 8) (first condition)                                   *
   * p = +3 or -3 (mod 8) (second condition)                                  */
  bool pCond1Z2(const int p);
  bool pCond2Z2(const int p);

  /* The two conditions (for the Gaussian Integers):  *
   * p = +1 (mod 4) (first condition)                 *
   * p = +3 (mod 4) (second condition)                */
  bool pCond1GI(const int p);
  bool pCond2GI(const int p);

  /* Do prime factorization in Z[Sqrt[2]]. This takes 'in' as an *
   * element of Z[Sqrt[2]] and returns a list of primes that     *
   * divide 'in' (without multiplicity).                         */
  void factorZ2(const vec2i& in, vector<vec2i>& factorization);

  // Do prime factorization in the Gaussian integers.
  void factorGI(const vec2i& in, vector<vec2i>& factorization);
};

namespace ArithVisibility {
  /* Division test (square-free case) for the primes which *
   * satisfy the first condition.                          */
  bool divTest2Free1(const vec2i& in, const int p);

  /* Division test (square-free case) for the primes which *
   * satisfy the second condition.                         */
  bool divTest2Free2(const vec2i& in, const int p);

  // Check for square-free 'visibility' of an element of Z[Sqrt[2]].
  bool visibility2Free(const vec2i& in);

  /* Let x = in / c, an element of Q(Sqrt[2]), then this computes    *
   * the denominator in the Fourier module Z[Sqrt[2]] * (Sqrt[2]/4). */
  vec2i denomZ2Fourier(const vec2i& in, const int in_c);

  /* Compute the intensity of 'denom', an element of Z[Sqrt[2]].       *
   * The input are going to be denominators of elements of Q(Sqrt[2]). */
  double intensityZ2(const vec2i& denom);

  // Division tests for the primes in the cube-free case.
  bool divTest3Free1(const vec2i& in, const int p);
  bool divTest3Free2(const vec2i& in, const int p);

  // Check for cube-free 'visibility' of an element of Z[Sqrt[2]].
  bool visibility3Free(const vec2i& in);

  class vec2iq {
  private:
    vec2i numerator;
    int denominator;

  public:
    vec2iq() {}
    vec2iq(int a, int b, int c) : numerator(a, b), denominator(c) {}

    bool operator==(const vec2iq& v) const {
      return ((numerator == v.numerator) && (denominator == v.denominator));
    }

    /* Lexicographic ordering, this is needed to use sorting *
     * algorithms of STL containers.                         */
    bool operator<(const vec2iq& v) const {
      if (numerator < v.numerator) {
        return true;
      }

      if (numerator == v.numerator) {
        return (denominator < v.denominator);
      }

      return false;
    }

    vec2d minkowskiQ2() const {
      const double inv = 1.0 / double(denominator);

      return (numerator.minkowskiZ2() * inv);
    }

    const vec2i& getNumerator() const {
      return numerator;
    }

    int getDenominator() const {
      return denominator;
    }
  };

  typedef double (&scalefunc)(double);
  typedef bool (&clipfunc)(const vec2d&);

  class bragg {
  private:
    vec2d position;
    double intensity;

  public:
    bragg() {}
    bragg(double x, double y, double i) : position(x, y), intensity(i) {}
    bragg(const vec2d& p, double i) : position(p), intensity(i) {}

    void apply(scalefunc f) {
      intensity = f(intensity);
    }

    const vec2d& getPosition() const {
      return position;
    }

    double getIntensity() const {
      return intensity;
    }
  };

  void diffractionZ2(const vector<vec2iq>& in, vector<bragg>& out,
                     clipfunc f);

   /* Box clipping test that checks if the point 'x' is contained in an *
    * arrangement of fundamental domains of the reciprocal lattice.     *
    * Currently two times two domains (around the origin) are cut.      */
   bool clipFundamental(const vec2d& x);
};

ostream& operator<<(ostream &os, const ArithVisibility::vec2iq& v);

// Format output so that we can use it in SAGE
ostream& operator<<(ostream &os, const ArithVisibility::bragg& b);

void vTableZ2(const uint r, Common::vec2ilist& table);
void vqTableRecipZ2(const uint r, const uint s,
                    vector<ArithVisibility::vec2iq>& table);

#endif
