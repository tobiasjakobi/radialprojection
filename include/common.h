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

#ifndef _COMMON_H_
#define _COMMON_H_

#include <cstdlib>
#include <cmath>
#include <cassert>

#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>
#include <algorithm>

#ifdef COMMON_USE_SSE
#include "alignment_allocator.h"
#include <xmmintrin.h>
#include <smmintrin.h>
#endif

typedef unsigned int uint;
typedef unsigned short ushort;
typedef unsigned char ubyte;
typedef unsigned long long ullong;

using namespace std;

class vec2d;

namespace Constants {
  /* The unit in Z[tau], with tau the golden mean.
   * In particular we have tau = unitGM here. */
  static const double unitGM = (1.0 + sqrt(5.0)) * 0.5;

  // The unit in Z[Sqrt[2]], which is the silver mean.
  static const double unitZ2 = 1.0 + sqrt(2.0);

  // The unit in Z[Sqrt[3]].
  static const double unitZ3 = 2.0 + sqrt(3.0);

  static const double pi = atan(1.0) * 4.0;
  static const double eps = numeric_limits<double>::epsilon();

  // Zeta(2) in different number fields
  static const double zetaZ2 = 48.0 * sqrt(2.0) / (pi * pi * pi * pi);
  static const double zetaES = 1.28519095548414940291751179870;
  static const double zetaGM = 1.1616711956186385497585826363320589131;

  static const double catalanC = 0.915965594177219015054603514932;
};

class vec2i {
public:
  int x, y;

  vec2i() {}
  vec2i(int a, int b) : x(a), y(b) {}

  vec2i(const vec2i& v) : x(v.x), y(v.y) {}

  bool operator==(const vec2i& v) const {
    return ((x == v.x) && (y == v.y));
  }

  /* Lexicographic ordering, this is needed to use sorting *
   * algorithms of STL containers.                         */
  bool operator<(const vec2i& v) const {
    if (x < v.x) {
      return true;
    }

    if (x == v.x) {
      return (y < v.y);
    }

    return false;
  }

  vec2i operator+(const vec2i& v) const {
    return vec2i(x + v.x, y + v.y);
  }

  void set(int a, int b) {
    x = a; y = b;
  }

  bool isZero() const {
    return (x == 0 && y == 0);
  }

  int preNormGI() const {
    return x*x + y*y;
  }

  int normGI() const {
    return abs(x*x + y*y);
  }

  int preNormES() const {
    return x*x - x*y + y*y;
  }

  int normES() const {
    return abs(x*x - x*y + y*y);
  }

  int preNormZ2() const {
    return x*x - 2*y*y;
  }

  int normZ2() const {
    return abs(x*x - 2*y*y);
  }

  int preNormZTau() const {
    return x*x + x*y - y*y;
  }

  int normZTau() const {
    return abs(x*x + x*y - y*y);
  }

  int preNormZ3() const {
    return x*x - 3*y*y;
  }

  int normZ3() const {
    return abs(x*x - 3*y*y);
  }

  double angle() const {
    return atan2(double(y), double(x));
  }

  vec2i conjZ2() const {
    return vec2i(this->x, -this->y);
  }

  vec2i conjGI() const {
    return vec2i(this->x, -this->y);
  }

  vec2i conjES() const {
    return vec2i(this->x - this->y, -this->y);
  }

  vec2i conjGM() const {
    return vec2i(this->x + this->y, -this->y);
  }

  // Squaring in Z[Sqrt[2]]
  vec2i squareZ2() const {
    return vec2i(x*x + 2*y*y, 2*x*y);
  }

  // Cubing in Z[Sqrt[2]]
  vec2i cubeZ2() const {
    return vec2i(x*x*x + 6*x*y*y, 3*x*x*y + 2*y*y*y);
  }

  // Squaring in the Gaussian Integers
  vec2i squareGI() const {
    return vec2i(x*x - y*y, 2*x*y);
  }

  // Cubing in the Gaussian Integers
  vec2i cubeGI() const {
    return vec2i(x*x*x - 3*x*y*y, 3*x*x*y - y*y*y);
  }

  // Squaring in the Eisenstein Integers
  vec2i squareES() const {
    return vec2i(x*x - y*y, 2*x*y - y*y);
  }

  // Cubing in the Eisenstein Integers
  vec2i cubeES() const {
    return vec2i(x*x*x - 3*x*y*y + y*y*y, 3*x*(x-y)*y);
  }

  // Squaring in Z[tau] (tau being the golden mean)
  vec2i squareGM() const {
    return vec2i(x*x + y*y, y*(2*x + y));
  }

  // Cubing in Z[tau] (tau being the golden mean)
  vec2i cubeGM() const {
    return vec2i(x*x*x + 3*x*y*y + y*y*y, y*(3*x*x + 3*x*y + 2*y*y));
  }

  vec2i positiveZ2() const {
    if (double(x) + double(y) * sqrt(2.0) < 0)
      return vec2i(-x, -y);
    else
      return vec2i(x, y);
  }

  bool isPositiveZ2() const {
    return (double(x) + double(y) * sqrt(2.0) >= 0);
  }

  vec2i positiveGM() const {
    if (double(x) + double(y) * Constants::unitGM < 0)
      return vec2i(-x, -y);
    else
      return vec2i(x, y);
  }

  bool isPositiveGM() const {
    return (double(x) + double(y) * Constants::unitGM >= 0);
  }

  bool isDiv(const int d) const {
    return ((x % d == 0) && (y % d == 0));
  }

  // Check if 'this' can be divided by 'd' in Z[Sqrt[2]].
  bool isDivZ2(const vec2i& d) const {
    const int norm = d.preNormZ2();
    const int a = this->x * d.x - 2 * this->y * d.y;
    const int b = d.x * this->y - this->x * d.y;

    return ((a % norm == 0) && (b % norm == 0));
  };

  // Divide 'this' by 'd' in Z[Sqrt[2]], assuming that this is defined.
  vec2i divZ2(const vec2i& d) const {
    const int norm = d.preNormZ2();
    const int a = this->x * d.x - 2 * this->y * d.y;
    const int b = d.x * this->y - this->x * d.y;

    return vec2i(a / norm, b / norm);
  };

  // Check if 'this' can be divided by 'd' in the Gaussian Integers.
  bool isDivGI(const vec2i& d) const {
    const int norm = d.preNormGI();
    const int a = this->x * d.x + this->y * d.y;
    const int b = this->y * d.x - this->x * d.y;

    return ((a % norm == 0) && (b % norm == 0));
  }

  // Divide 'this' by 'd' in the Gaussian Integers, assuming that this is defined.
  vec2i divGI(const vec2i& d) const {
    const int norm = d.preNormGI();
    const int a = this->x * d.x + this->y * d.y;
    const int b = this->y * d.x - this->x * d.y;

    return vec2i(a / norm, b / norm);
  }

  // Check if 'this' can be divided by 'd' in the Eisenstein Integers.
  bool isDivES(const vec2i& d) const {
    const int norm = d.preNormES();
    const int a = this->x * d.x + this->y * d.y - this->x * d.y;
    const int b = this->y * d.x - this->x * d.y;

    return ((a % norm == 0) && (b % norm == 0));
  }

  // Divide 'this' by 'd' in the Eisenstein Integers, assuming that this is defined.
  vec2i divES(const vec2i& d) const {
    const int norm = d.preNormES();
    const int a = this->x * d.x + this->y * d.y - this->x * d.y;
    const int b = this->y * d.x - this->x * d.y;

    return vec2i(a / norm, b / norm);
  }

  // Check if 'this' can be divided by 'd' in Z[tau] (tau being the golden mean).
  bool isDivGM(const vec2i& d) const {
    const int norm = d.preNormZTau();
    const int a = this->x * (d.x + d.y) - this->y * d.y;
    const int b = d.x * this->y - this->x * d.y;

    return ((a % norm == 0) && (b % norm == 0));
  }

  // Divide 'this' by 'd' in in Z[tau], assuming that this is defined.
  vec2i divGM(const vec2i& d) const {
    const int norm = d.preNormZTau();
    const int a = this->x * (d.x + d.y) - this->y * d.y;
    const int b = d.x * this->y - this->x * d.y;

    return vec2i(a / norm, b / norm);
  }

  /* Multiply an element of Z[tau] with the unit tau^k and return the result. */
  vec2i multUnitGM(int k) const;

  /* Reduce an element of Z[tau] (tau being the golden mean) into a the     *
   * fundamental domain [1, tau) (or (-tau, -1] if the element is negative) *
   * by multiplication with the units tau and tau^{-1}.                     */
  vec2i reduceGM(int& k) const;

  // Same as above, but only in Z[Sqrt[2]].
  vec2i multUnitZ2(int k) const;
  vec2i reduceZ2(int& k) const;

  vec2d transTriToR2() const;

  vec2d transGenericToR2(const vec2d& v) const;

  vec2d minkowskiZ2() const;
  vec2d minkowskiES() const;
  vec2d minkowskiGM() const;

  bool coprime() const;

  /* Return primitive version of the element, by dividing *
   * the coordinates by the (positive) GCD.               */
  vec2i primitive() const;

};

class vec2d {
private:
#ifdef COMMON_USE_SSE
  union {
    double a[2];
    __m128d vsse;
  };
#else
  double a[2];
#endif

public:
  vec2d() {}
  vec2d(double x, double y) {
    a[0] = x; a[1] = y;
  }

#ifdef COMMON_USE_SSE
  vec2d(const vec2d& v) : vsse(v.vsse) {}
#else
  vec2d(const vec2d& v) {
    a[0] = v.a[0]; a[1] = v.a[1];
  }
#endif

#ifdef COMMON_USE_SSE
  vec2d(const __m128d& sse) : vsse(sse) {}
#endif

  vec2d operator*(double scale) const {
    return vec2d(scale * a[0], scale * a[1]);
  }

  vec2d operator-(const vec2d& v) const {
    return vec2d(a[0] - v.a[0], a[1] - v.a[1]);
  }

  vec2d operator+(const vec2d& v) const {
    return vec2d(a[0] + v.a[0], a[1] + v.a[1]);
  }

  vec2d& operator+=(const vec2d& v) {
    a[0] += v.a[0]; a[1] += v.a[1];
    return *this;
  }

  vec2d operator-() const {
    return vec2d(-a[0], -a[1]);
  }

  double operator[](uint i) const {
    assert(i < 2);
    return a[i];
  }

  double& operator[](uint i) {
    assert(i < 2);
    return a[i];
  }

  void set(double x, double y) {
    a[0] = x; a[1] = y;
  }

  vec2d abs() const {
#ifdef COMMON_USE_SSE
    return vec2d(_mm_max_pd(_mm_sub_pd(_mm_setzero_pd(), vsse), vsse));

    /*
     * alternate version using sign bitmasking:
     * static const __m128d sign_mask = _mm_set1_pd(-0.); // -0. = 1 << 63
     * return _mm_andnot_pd(sign_mask, vsse); // !sign_mask & x
    */
#else
    return vec2d(std::abs(a[0]), std::abs(a[1]));
#endif
  }

  double lengthSquared() const {
#ifdef COMMON_USE_SSE
    /* Compute dot-product with itself, write result to first *
     * component and extract the component as a double float. */
    return _mm_cvtsd_f64(_mm_dp_pd(vsse, vsse, 0x31));
#else
    return a[0] * a[0] + a[1] * a[1];
#endif
  }

  double length() const {
    return sqrt(lengthSquared());
  }

  bool inSectorL8() const {
    return (a[0] >= 0.0 && a[1] >= a[0]);
  }

  bool inUpperHalfplane() const {
    return (a[1] >= 0.0);
  }

  bool inFirstQuadrant() const {
    return (a[0] >= 0.0 && a[1] >= 0.0);
  }

  bool inFirstQuadOpen() const {
    return (a[0] > 0.0 && a[1] > 0.0);
  }

  // Checks for phi(x,y) <= 2*pi/6 (60 degree):
  bool inSectorL3() const {
    return (a[1]/a[0] <= sectorL3);
  }

  // Checks for phi(x,y) <= 2*pi/5 (72 degree):
  bool inSectorL5() const {
    return (a[1]/a[0] <= sectorL5);
  }

  // Checks for phi(x,y) <= 2*pi/6 (60 degree):
  bool inSectorL12() const {
    return (a[1]/a[0] <= sectorL12);
  }

  // Checks for phi(x,y) <= 2*pi/7 (51 degree):
  bool inSectorL7() const {
    return (a[1]/a[0] <= sectorL7);
  }

  vec2d reduceIntoSectorL12() const {
    const double absvals[2] = {std::abs(a[0]), std::abs(a[1])};

    return (absvals[0] >= absvals[1]) ? vec2d(absvals[0], absvals[1]) :
                                        vec2d(absvals[1], absvals[0]);
  }

  double angle() const {
    return atan2(a[1], a[0]);
  }

  double dot(const vec2d& v) const {
    return a[0] * v.a[0] + a[1] * v.a[1];
  }

  vec2d normalize() const {
    const double invlen = 1.0 / this->length();
    return vec2d(a[0] * invlen, a[1] * invlen);
  }

  // x = cos(alpha), y = sin(alpha)
  vec2d applyRotation(const double x, const double y) const {
    return vec2d(a[0]*x - a[1]*y, a[0]*y + a[1]*x);
  }

private:
  static const double sectorL3;
  static const double sectorL5;
  static const double sectorL12;

  static const double sectorL7;

};

class vec4i {
private:
#ifdef COMMON_USE_SSE
  union {
    int a[4];
    __m128i vsse;
  };
#else
  int a[4];
#endif

public:
  vec4i() {}
  vec4i(int x0, int x1, int x2, int x3) {
#ifdef COMMON_USE_SSE
    vsse = _mm_set_epi32(x0, x1, x2, x3);
#else
    a[0] = x0; a[1] = x1;
    a[2] = x2; a[3] = x3;
#endif
  }

  vec4i(const vec2i& x01, const vec2i& x23) {
    a[0] = x01.x; a[1] = x01.y;
    a[2] = x23.x; a[3] = x23.y;
  }

#ifdef COMMON_USE_SSE
  vec4i(const vec4i& v) : vsse(v.vsse) {}
#else
  vec4i(const vec4i& v) {
    a[0] = v.a[0]; a[1] = v.a[1];
    a[2] = v.a[2]; a[3] = v.a[3];
  }
#endif

#ifdef COMMON_USE_SSE
  vec4i(const __m128i& sse) : vsse(sse) {}
#endif

  vec4i operator+(const vec4i& v) const {
#ifdef COMMON_USE_SSE
    return vec4i(_mm_add_epi32(vsse, v.vsse));
#else
    return vec4i(a[0] + v.a[0], a[1] + v.a[1], a[2] + v.a[2], a[3] + v.a[3]);
#endif
  }

  vec4i operator-(const vec4i& v) const {
#ifdef COMMON_USE_SSE
    return vec4i(_mm_sub_epi32(vsse, v.vsse));
#else
    return vec4i(a[0] - v.a[0], a[1] - v.a[1], a[2] - v.a[2], a[3] - v.a[3]);
#endif
  }

  /* Lexicographic ordering, this is needed to use sorting *
   * algorithms of STL containers.                         */
  bool operator<(const vec4i& v) const {
#ifdef COMMON_USE_SSE
    uint mask = 0x000f;
    const uint ltmask = _mm_movemask_epi8(_mm_cmplt_epi32(vsse, v.vsse));
    const uint eqmask = _mm_movemask_epi8(_mm_cmpeq_epi32(vsse, v.vsse));

    if (ltmask & mask) return true;
    if (!(eqmask & mask)) return false;
    mask <<= 4;

    if (ltmask & mask) return true;
    if (!(eqmask & mask)) return false;
    mask <<= 4;

    if (ltmask & mask) return true;
    if (!(eqmask & mask)) return false;
    mask <<= 4;

    return (ltmask & mask);
#else
    if (a[0] < v[0]) return true;

    if (a[0] == v[0]) {
      if (a[1] < v[1]) return true;

       if (a[1] == v[1]) {
         if (a[2] < v[2]) return true;

         if (a[2] == v[2]) {
           return (a[3] < v[3]);
         }
       }
    }

    return false;
#endif
  }

  bool operator==(const vec4i& v) const {
#ifdef COMMON_USE_SSE
    const __m128i x = _mm_xor_si128(vsse, v.vsse);
    return _mm_testz_si128(x, x);
#else
    return (a[0] == v.a[0] && a[1] == v.a[1] && a[2] == v.a[2] && a[3] == v.a[3]);
#endif
  }

  bool operator!=(const vec4i& v) const {
    return !(*this == v);
  }

  int operator[](uint i) const {
    assert(i < 4);
    return a[i];
  }

  int& operator[](uint i) {
    assert(i < 4);
    return a[i];
  }

  bool isZero() const {
#ifdef COMMON_USE_SSE
    return _mm_testz_si128(vsse, vsse);
#else
    return (a[0] == 0 && a[1] == 0 && a[2] == 0 && a[3] == 0);
#endif
  }

  // Check if first component of direct-sum is zero
  bool isFirstZero() const {
    return (a[0] == 0 && a[1] == 0);
  }

  // Check if second component of direct-sum is zero
  bool isSecondZero() const {
    return (a[2] == 0 && a[3] == 0);
  }

  vec2i getFirst() const {
    return vec2i(a[0], a[1]);
  }

  vec2i getSecond() const {
    return vec2i(a[2], a[3]);
  }

  const vec2i& getFirstDirect() const {
    return *reinterpret_cast<const vec2i*>(&a[0]);
  }

  const vec2i& getSecondDirect() const {
    return *reinterpret_cast<const vec2i*>(&a[2]);
  }

  // Transform from L8 to direct-sum representation
  vec4i transL8ToDirect() const {
    return vec4i(a[0] - a[2], -a[3], a[1] + a[3], a[2]);
  }

  // Transform from L5 to direct-sum representation
  vec4i transL5ToDirect() const {
    return vec4i(a[0] - a[2] + a[3], -a[3], a[1] - a[2] + a[3], a[2] - a[3]);
  }

  // Transform from L12 to direct-sum representation
  vec4i transL12ToDirect() const {
    return vec4i(a[0] - a[2], -a[3], a[1] + 2*a[3], a[2]);
  }

  // Inverse of transL8ToDirect()
  vec4i transDirectToL8() const {
    return vec4i(a[0] + a[3], a[1] + a[2], a[3], -a[1]);
  }

  // Inverse of transL5ToDirect()
  vec4i transDirectToL5() const {
    return vec4i(a[0] + a[3], a[2] + a[3], a[3] - a[1], -a[1]);
  }

  // Inverse of transL12ToDirect()
  vec4i transDirectToL12() const {
    return vec4i(a[0] + a[3], 2*a[1] + a[2], a[3], -a[1]);
  }

  // Perform a multiplication in L12
  vec4i multL12(const vec4i& vec) const {
    const int k0 = a[1] * vec.a[3] + a[3] * vec.a[1] + a[2] * vec.a[2]; // xi^4
    const int k1 = a[2] * vec.a[3] + a[3] * vec.a[2]; // xi^5

    return vec4i(
      a[0] * vec.a[0] - a[3] * vec.a[3] - k0,
      a[0] * vec.a[1] + a[1] * vec.a[0] - k1,
      a[0] * vec.a[2] + a[2] * vec.a[0] + a[1] * vec.a[1] + k0,
      a[0] * vec.a[3] + a[3] * vec.a[0] + a[1] * vec.a[2] + a[2] * vec.a[1] + k1
    );
  }

  // Conjugation in L12
  vec4i conjL12() const {
    return vec4i(a[0] + a[2], a[1], -a[2], -a[1] - a[3]);
  }

  uint kappaL5() const {
    const int sum = (a[0] + a[1] + a[2] + a[3]) % 5;
    return (sum < 0) ? (5 + sum) : sum;
  }

  vec2d orthProjL8() const{
#ifdef COMMON_USE_SSE
    static const __m128 v1 = _mm_set_ps(1.0f/sqrtf(2.0f), -0.5f, 0.0f, 0.5f);
    static const __m128 v2 = _mm_set_ps(0.0f, -0.5f, 1.0f/sqrtf(2.0f), -0.5f);

    // convert int[4] to float[4]
    const __m128 temp = _mm_cvtepi32_ps(vsse);

    /*
     * _mm_dp_ps: dot-product on float[4] (0xAB is the read/write-mask)
     *            A = read-mask (f = all components)
     *            B = write-mask (1 = write result to first component,
     *                            2 = write result to 2nd component)
     *            (issue: this instruction is slow!)
     * _mm_add_ps: add two float[4]
     * _mm_cvtps_pd: convert first 2 components of a float[4] to double[2]
     */
    return vec2d(_mm_cvtps_pd(
      _mm_add_ps(_mm_dp_ps(temp, v1, 0xf1), _mm_dp_ps(temp, v2, 0xf2))));
#else
    static const double v1[4] = {1.0/sqrt(2.0), -0.5, 0.0, 0.5};
    static const double v2[4] = {0.0, -0.5, 1.0/sqrt(2.0), -0.5};

    return vec2d(a[0] * v1[0] + a[1] * v1[1] + a[2] * v1[2] + a[3] * v1[3],
                 a[0] * v2[0] + a[1] * v2[1] + a[2] * v2[2] + a[3] * v2[3]);
#endif
  }

  vec2d paraProjL8() const {
    static const double v1[4] = {1.0/sqrt(2.0), 0.5, 0.0, -0.5};
    static const double v2[4] = {0.0, 0.5, 1.0/sqrt(2.0), 0.5};

    return vec2d(a[0] * v1[0] + a[1] * v1[1] + a[2] * v1[2] + a[3] * v1[3],
                 a[0] * v2[0] + a[1] * v2[1] + a[2] * v2[2] + a[3] * v2[3]);
  }

  vec2d orthProjShiftL8() const {
    return orthProjL8() - shift;
  }

  vec2d orthProjShiftL8(const double scale) const {
    return (orthProjL8() * scale) - shift;
  }

  vec2d orthProjL5() const {
    static const double v1[4] = {1.0, 0.25 * (-1.0 - sqrt(5.0)),
                          0.25 * (-1.0 + sqrt(5.0)), 0.25 * (-1.0 + sqrt(5.0))};
    static const double v2[4] = {0.0, sqrt((5.0 - sqrt(5.0)) / 8.0),
                          -sqrt((5.0 + sqrt(5.0)) / 8.0), sqrt((5.0 + sqrt(5.0)) / 8.0)};

    return vec2d(a[0] * v1[0] + a[1] * v1[1] + a[2] * v1[2] + a[3] * v1[3],
                 a[0] * v2[0] + a[1] * v2[1] + a[2] * v2[2] + a[3] * v2[3]);
  }

  vec2d paraProjL5() const {
    static const double v1[4] = {1.0, 0.25 * (-1.0 + sqrt(5.0)),
                          0.25 * (-1.0 - sqrt(5.0)), 0.25 * (-1.0 - sqrt(5.0))};
    static const double v2[4] = {0.0, sqrt((5.0 + sqrt(5.0)) / 8.0),
                          sqrt((5.0 - sqrt(5.0)) / 8.0), -sqrt((5.0 - sqrt(5.0)) / 8.0)};

    return vec2d(a[0] * v1[0] + a[1] * v1[1] + a[2] * v1[2] + a[3] * v1[3],
                 a[0] * v2[0] + a[1] * v2[1] + a[2] * v2[2] + a[3] * v2[3]);
  }

  vec2d orthProjShiftL5() const {
    return orthProjL5() - shift;
  }

  vec2d orthProjShiftL5(const double scale) const {
    return (orthProjL5() * scale) - shift;
  }

  /* The direction of the shift has to be inverted, when testing for visibility. */
  vec2d orthProjShiftL5(const double scale, bool invert) const {
    return (orthProjL5() * scale) - shift * (invert ? -1.0 : +1.0);
  }

  vec2d orthProjL12() const {
    static const double v1[4] = {1.0, -0.5 * sqrt(3.0), 0.5, 0.0};
    static const double v2[4] = {0.0, 0.5, -0.5 * sqrt(3.0), 1.0};

    return vec2d(a[0] * v1[0] + a[1] * v1[1] + a[2] * v1[2] + a[3] * v1[3],
                 a[0] * v2[0] + a[1] * v2[1] + a[2] * v2[2] + a[3] * v2[3]);
  }

  vec2d paraProjL12() const {
    static const double v1[4] = {1.0, 0.5 * sqrt(3.0), 0.5, 0.0};
    static const double v2[4] = {0.0, 0.5, 0.5 * sqrt(3.0), 1.0};

    return vec2d(a[0] * v1[0] + a[1] * v1[1] + a[2] * v1[2] + a[3] * v1[3],
                 a[0] * v2[0] + a[1] * v2[1] + a[2] * v2[2] + a[3] * v2[3]);
  }

  vec2d orthProjShiftL12() const {
    return orthProjL12() - shift;
  }

  vec2d orthProjShiftL12(const double scale, bool invert) const {
    return (orthProjL12() * scale) - shift * (invert ? -1.0 : +1.0);
  }

  // For comments see vec4s::directL10ToUnique in this header.
  vec4i directL8ToUnique() const;
  vec4i directL5ToUnique() const;

  static vec2d shift;

};

class vec4s {
private:
  short a[4];

public:
  vec4s() {}
  vec4s(short x0, short x1, short x2, short x3) {
    a[0] = x0; a[1] = x1;
    a[2] = x2; a[3] = x3;
  }

  vec4s(const vec2i& x0, const vec2i& x1) {
    a[0] = x0.x; a[1] = x0.y;
    a[2] = x1.x; a[3] = x1.y;
  }

  vec4s operator+(const vec4s& v) const {
    return vec4s(a[0] + v.a[0], a[1] + v.a[1],
                 a[2] + v.a[2], a[3] + v.a[3]);
  }

  vec4s operator-(const vec4s& v) const {
    return vec4s(a[0] - v.a[0], a[1] - v.a[1],
                 a[2] - v.a[2], a[3] - v.a[3]);
  }

  vec4s operator*(short val) const {
    return vec4s(a[0] * val, a[1] * val,
                 a[2] * val, a[3] * val);
  }

  /* Lexicographic ordering, this is needed to use sorting *
   * algorithms of STL containers.                         */
  bool operator<(const vec4s& v) const {
    if (a[0] < v[0]) return true;

    if (a[0] == v[0]) {
      if (a[1] < v[1]) return true;

       if (a[1] == v[1]) {
         if (a[2] < v[2]) return true;

         if (a[2] == v[2]) {
           return (a[3] < v[3]);
         }
       }
    }

    return false;
  }
  
  bool operator==(const vec4s& v) const {
    return (a[0] == v.a[0] && a[1] == v.a[1] &&
            a[2] == v.a[2] && a[3] == v.a[3]);
  }

  short operator[](uint i) const {
    assert(i < 4);
    return a[i];
  }

  short& operator[](uint i) {
    assert(i < 4);
    return a[i];
  }

  vec4s& operator=(const vec4s& v) {
    a[0] = v.a[0]; a[1] = v.a[1];
    a[2] = v.a[2]; a[3] = v.a[3];
    return *this;
  }

  void set(short x0, short x1, short x2, short x3) {
    a[0] = x0; a[1] = x1;
    a[2] = x2; a[3] = x3;
  }

  bool isZero() const {
    return (a[0] == 0 && a[1] == 0 && a[2] == 0 && a[3] == 0);
  }

  // Check if first component of direct-sum is zero
  bool isFirstZero() const {
    return (a[0] == 0 && a[1] == 0);
  }

  // Check if second component of direct-sum is zerp
  bool isSecondZero() const {
    return (a[2] == 0 && a[3] == 0);
  }

  // Interpret as L5 element and return the element scaled by tau (golden mean)
  vec4s scaleTauL5() const {
    return vec4s(a[1] - a[3], a[1] + a[2] - a[3], -a[0] + a[1] + a[2], -a[0] + a[2]);
  }

  // Interpret as L5 element and "shift" by xi^r with
  // xi being the 10-th root of unity (xi = exp(2*Pi*I/10))
  vec4s shiftL5ByL10(ushort r) const {
    assert(r < 10);

    short mult = (r >= 5) ? -1 : 1;
    ushort s = r % 5;
    mult *= ((s % 2) ? -1 : 1);
    s = (s*3) % 5;

    vec4s temp(*this);
    for (uint i = 0; i < s; ++i) {
      temp.set(-temp[3], temp[0] - temp[3],
               temp[1] - temp[3], temp[2] - temp[3]);
    }

    return (temp * mult);
  }

  // Transform from L10 to direct-sum representation
  vec4s transL10ToDirect() const {
    const short b[4] = {a[0] + a[3], a[2] + a[3], a[3], a[3] - a[1]}; // L10 -> L5

    return vec4s(b[0] - b[2] + b[3], -b[3], b[1] - b[2] + b[3], b[2] - b[3]);
  }

  void transL10ToDirect(vec2i& v0, vec2i& v1) const {
    const short b[4] = {a[0] + a[3], a[2] + a[3], a[3], a[3] - a[1]}; // L10 -> L5

    v0.x = b[0] - b[2] + b[3];
    v0.y = -b[3];
    v1.x = b[1] - b[2] + b[3];
    v1.y = b[2] - b[3];
  }

  // Transform from L5 to direct-sum representation
  vec4s transL5ToDirect() const {
    return vec4s(a[0] - a[2] + a[3], -a[3], a[1] - a[2] + a[3], a[2] - a[3]);
  }

  // Counterpart to 'transL10ToDirect'
  vec4s transDirectToL10() const {
    return vec4s(a[0] + a[1], a[3], a[1] + a[2], -a[1] + a[3]);
  }

  // Counterpart to 'transL5ToDirect'
  vec4s transDirectToL5() const {
    return vec4s(a[0] + a[3], a[2] + a[3], -a[1] + a[3], -a[1]);
  }

  // Transfrom from L10 to "physical" R2 space
  vec2d transL10ToR2() const {
    static const double v1[4] = {1.0, (1.0 + sqrt(5.0)) / 4.0,
                                 (-1.0 + sqrt(5.0)) / 4.0, (1.0 - sqrt(5.0)) / 4.0};
    static const double v2[4] = {0.0, sqrt((5.0 - sqrt(5.0)) / 8.0),
                                 sqrt((5.0 + sqrt(5.0)) / 8.0), sqrt((5.0 + sqrt(5.0)) / 8.0)};

    return vec2d(a[0] * v1[0] + a[1] * v1[1] + a[2] * v1[2] + a[3] * v1[3],
                 a[0] * v2[0] + a[1] * v2[1] + a[2] * v2[2] + a[3] * v2[3]);
  }

  vec2d transL5ToR2() const {
    static const double v1[4] = {1.0, (sqrt(5.0) - 1.0) / 4.0,
                                 -(sqrt(5.0) + 1.0) / 4.0, -(sqrt(5.0) + 1.0) / 4.0};
    static const double v2[4] = {0.0, sqrt((5.0 + sqrt(5.0)) / 8.0),
                                 sqrt((5.0 - sqrt(5.0)) / 8.0), -sqrt((5.0 - sqrt(5.0)) / 8.0)};

    return vec2d(a[0] * v1[0] + a[1] * v1[1] + a[2] * v1[2] + a[3] * v1[3],
                 a[0] * v2[0] + a[1] * v2[1] + a[2] * v2[2] + a[3] * v2[3]);
  }

  vec2d directL10ToR2() const {
    static const double zeta5[2] = {(-1.0 + sqrt(5.0)) / 4.0,
                                    sqrt((5.0 + sqrt(5.0)) / 8.0)};

    const double z0 = double(a[0]) + double(a[1]) * Constants::unitGM;
    const double z1 = double(a[2]) + double(a[3]) * Constants::unitGM;

    return vec2d(z0 + z1 * zeta5[0],
                 z1 * zeta5[1]);
  }

  /* Take an element 'x' in direct-sum (L10) representation and return the *
   * primitive reduced version of 'x'. This makes 'x' unique, in the sense *
   * that two elements lie on the same ray through zero iff after this     *
   * transformation they are equal in terms of coordinates.                */
  vec4s directL10ToUnique() const;

};

class vec8s {
private:
  short a[8];

public:
  vec8s() {
    for (uint i = 0; i < 8; ++i) {
      a[i] = 0;
    }
  }

  vec8s(short x0, short x1, short x2, short x3,
        short x4, short x5, short x6, short x7) {
    a[0] = x0; a[1] = x1;
    a[2] = x2; a[3] = x3;
    a[4] = x4; a[5] = x5;
    a[6] = x6; a[7] = x7;
  }

  vec8s(uint i) {
    static const short roots[10][8] = {
      {1, 0, 0, 0, 0, 0, 0, 0},
      {0, 1, 0, 0, 0, 0, 0, 0},
      {0, 0, 1, 0, 0, 0, 0, 0},
      {0, 0, 0, 1, 0, 0, 0, 0},
      {0, 0, 0, 0, 1, 0, 0, 0},
      {0, 0, 0, 0, 0, 1, 0, 0},
      {0, 0, 0, 0, 0, 0, 1, 0},
      {0, 0, 0, 0, 0, 0, 0, 1},
      {-1, 0, 1, 0, -1, 0, 1, 0},
      {0, -1, 0, 1, 0, -1, 0, 1}
    };

    assert(i < 20);
    short mult = 1;
    if (i >= 10) {
      i -= 10;
      mult = -1;
    }

    for (uint j = 0; j < 8; ++j) {
      a[j] = roots[i][j] * mult;
    }
  }

  vec8s(const vec8s& v) {
    for (uint i = 0; i < 8; ++i) {
      a[i] = v.a[i];
    }
  }

  vec8s operator+(const vec8s& v) const {
    return vec8s(a[0] + v.a[0], a[1] + v.a[1],
                 a[2] + v.a[2], a[3] + v.a[3],
                 a[4] + v.a[4], a[5] + v.a[5],
                 a[6] + v.a[6], a[7] + v.a[7]);
  }

  vec8s operator*(short val) const {
    return vec8s(a[0] * val, a[1] * val,
                 a[2] * val, a[3] * val, 
                 a[4] * val, a[5] * val, 
                 a[6] * val, a[7] * val);
  }

  bool operator==(const vec8s& v) const {
    return (v.a[0] == a[0] && v.a[1] == a[1] &&
            v.a[2] == a[2] && v.a[3] == a[3] &&
            v.a[4] == a[4] && v.a[5] == a[5] &&
            v.a[6] == a[6] && v.a[7] == a[7]);
  }

  void set(short x0, short x1, short x2, short x3,
           short x4, short x5, short x6, short x7) {
    a[0] = x0; a[1] = x1;
    a[2] = x2; a[3] = x3;
    a[4] = x4; a[5] = x5;
    a[6] = x6; a[7] = x7;
  }

  bool isInL10(uint& even, uint& odd) const {
    bool state = false;

    if (a[1] == 0 && a[3] == 0 &&
        a[5] == 0 && a[7] == 0) {
      ++even;
      state = true;
    }

    if (a[0] == 0 && a[2] == 0 &&
        a[4] == 0 && a[6] == 0) {
      ++odd;
      state = true;
    }

    return state;
  }

  vec4s reduceToL10(bool even) const {
    return (even ? vec4s(a[0], a[2], a[4], a[6]) :
                   vec4s(a[1], a[3], a[5], a[7]));
  }

  vec8s lambdaScale() const {
    return vec8s(a[1] - a[7], 2*a[0] + a[2], a[1] + a[3] + a[7],
                 -a[0] + a[2] + a[4], a[3] + a[5] - a[7],
                 a[0] + a[4] + a[6], a[5] + 2*a[7], -a[0] + a[6]);
  }

  vec8s shift(ushort r) const {
    bool negate = false;
    r %= 20;
    if (r >= 10) {
      negate = true;
      r -= 10;
    }

    vec8s temp(*this);
    for (uint i = 0; i < r; ++i) {
      temp.set(-temp[7], temp[0], temp[1] + temp[7], temp[2],
               temp[3] - temp[7], temp[4], temp[5] + temp[7], temp[6]);
    }    

    return (negate ? temp * (-1) : temp);
  }

  vec8s& operator=(const vec8s& v) {
    for (uint i = 0; i < 8; ++i) {
      a[i] = v.a[i];
    }
    return *this;
  }

  short operator[](uint i) const {
    assert(i < 8);
    return a[i];
  }

  short& operator[](uint i) {
    assert(i < 8);
    return a[i];
  }

  // Transfrom from L20 to "physical" R2 space
  vec2d transL20ToR2() const {
    static const double v1[8] = {1.0, sqrt((5.0 + sqrt(5.0)) / 8.0),
                                 (1.0 + sqrt(5.0)) / 4.0, sqrt((5.0 - sqrt(5.0)) / 8.0),
                                 (-1.0 + sqrt(5.0)) / 4.0, 0.0,
                                 (1.0 - sqrt(5.0)) / 4.0, -sqrt((5.0 - sqrt(5.0)) / 8.0)};
    static const double v2[8] = {0.0, (-1.0 + sqrt(5.0)) / 4.0,
                                 sqrt((5.0 - sqrt(5.0)) / 8.0), (1.0 + sqrt(5.0)) / 4.0,
                                 sqrt((5.0 + sqrt(5.0)) / 8.0), 1.0,
                                 sqrt((5.0 + sqrt(5.0)) / 8.0), (1.0 + sqrt(5.0)) / 4.0};

    return vec2d(a[0]*v1[0] + a[1]*v1[1] + a[2]*v1[2] + a[3]*v1[3] +
                 a[4]*v1[4] + a[5]*v1[5] + a[6]*v1[6] + a[7]*v1[7],
                 a[0]*v2[0] + a[1]*v2[1] + a[2]*v2[2] + a[3]*v2[3] +
                 a[4]*v2[4] + a[5]*v2[5] + a[6]*v2[6] + a[7]*v2[7]);
  }

};

class vec2s {
private:
  short a[2];

public:
  vec2s() {
    a[0] = 0; a[1] = 0;
  }

  vec2s(short x0, short x1) {
    a[0] = x0; a[1] = x1;
  }

  vec2s(const vec2s& v) {
    a[0] = v.a[0];
    a[1] = v.a[1];
  }

  vec2s operator+(const vec2s& v) const {
    return vec2s(a[0] + v.a[0], a[1] + v.a[1]);
  }

  vec2s operator*(short val) const {
    return vec2s(a[0] * val, a[1] * val);
  }

  bool operator==(const vec2s& v) const {
    return (v.a[0] == a[0] && v.a[1] == a[1]);
  }

  /* Lexicographic ordering, this is needed to use sorting *
   * algorithms of STL containers.                         */
  bool operator<(const vec2s& v) const {
    if (a[0] < v.a[0]) {
      return true;
    }

    if (a[0] == v.a[0]) {
      return (a[1] < v.a[1]);
    }

    return false;
  }

  void set(short x0, short x1) {
    a[0] = x0; a[1] = x1;
  }

  vec2s shift(ushort r) const {
    r %= 4;

    short res[2] = {a[0], a[1]};
    for (uint i = 0; i < r; ++i) {
      const short temp = res[0];
      res[0] = -res[1];
      res[1] = temp;
    }

    return vec2s(res[0], res[1]);
  }

  vec2s& operator=(const vec2s& v) {
    a[0] = v.a[0];
    a[1] = v.a[1];

    return *this;
  }

  short operator[](uint i) const {
    assert(i < 2);
    return a[i];
  }

  short& operator[](uint i) {
    assert(i < 2);
    return a[i];
  }

  vec2d transZ2ToR2() const {
    return vec2d(double(a[0]), double(a[1]));
  }

};

/* The 'vec2iExt' class provides an efficient way to compute (correct)
 * visibility of a set of elements from Z^2.
 *
 * When constructing a vec2iExt object from a vec2i object the GCD g of
 * the two coordinates is computed and stored. The common factor g is then
 * removed from the coordinates.
 *
 * The compare operators just operate on the coordinate (and not the common
 * factor) data, so one can use the standard STL algorithms to sort the set
 * and remove duplicate (in the sense that the coordinate data is the same)
 * elements.
 *
 * To produce properly visible points, a normalization pass has to be
 * applied to the list after the sort. This ensures that the element with the
 * largest distance from zero survives.
 */
class vec2iExt {
public:
  vec2iExt(const vec2i& in);
  vec2iExt(const vec2iExt& c);

  bool operator==(const vec2iExt& rhs) const;
  bool operator<(const vec2iExt& rhs) const;

  operator vec2i() const;
  operator int() const;

  void set(int x) {e = x;}
  const vec2i& get() const {return v;}

private:
  vec2i v;
  int e;
};

class tilingEdge {
private:
  int a[2];

public:
  tilingEdge() {
    a[0] = 0; a[1] = 0;
  }

  tilingEdge(int x0, int x1) {
    a[0] = x0; a[1] = x1;
  }

  tilingEdge(const tilingEdge& e) {
    a[0] = e.a[0];
    a[1] = e.a[1];
  }

  int operator[](uint i) const {
    assert(i < 2);
    return a[i];
  }

  int& operator[](uint i) {
    assert(i < 2);
    return a[i];
  }

};

namespace Coprime {

  static inline double round(double number) {
    return number < 0.0 ? ceil(number - 0.5) : floor(number + 0.5);
  }

  uint gcdZFast(uint u, uint v); /* binary gcd implementation (taken from Wikipedia article) */
  int gcdZ(int a, int b);

  static inline bool coprimeZ(int a, int b) {
    return (gcdZFast(abs(a), abs(b)) == 1);
  }

  void gcdZTest(uint count);

  // Used in the octagonal case:
  vec2i moduloZ2(const vec2i& a, const vec2i& b);
  vec2i gcdZ2(const vec2i& a, const vec2i& b);

  static inline bool coprimeZ2(const vec2i& a, const vec2i& b) {
    return (gcdZ2(a, b).normZ2() == 1);
  }

  static inline bool coprimeZ2(const vec4i& v) {
    return (gcdZ2(v.getFirst(), v.getSecond()).normZ2() == 1);
  }

  void multZ2(const vec2i& a, const vec2i& b, vec2i& out);

  // Used in the decagonal case:
  vec2i moduloZTau(const vec2i& a, const vec2i& b);
  vec2i gcdZTau(const vec2i& a, const vec2i& b);

  static inline vec2i gcdZTau(const vec4s& v) {
    return gcdZTau(vec2i(v[0], v[1]), vec2i(v[2], v[3]));
  }

  static inline bool coprimeZTau(const vec2i& a, const vec2i& b) {
    return (gcdZTau(a, b).normZTau() == 1);
  }

  static inline bool coprimeZTau(const vec4i& v){
    return (gcdZTau(v.getFirst(), v.getSecond()).normZTau() == 1);
  }

  static inline bool coprimeZTau(const vec4s& v) {
    return (gcdZTau(vec2i(v[0], v[1]), vec2i(v[2], v[3])).normZTau() == 1);
  }

  void multZTau(const vec2i& a, const vec2i& b, vec2i& out);

  // Used in the dodecagonal case:
  vec2i moduloZ3(const vec2i& a, const vec2i& b);
  vec2i gcdZ3(const vec2i& a, const vec2i& b);

  static inline vec2i gcdZ3(const vec4i& v) {
    return gcdZ3(v.getFirst(), v.getSecond());
  }

  static inline bool coprimeZ3(const vec2i& a, const vec2i& b){
    return (gcdZ3(a, b).normZ3() == 1);
  }

  static inline bool coprimeZ3(const vec4i& v) {
    return (gcdZ3(v.getFirst(), v.getSecond()).normZ3() == 1);
  }

  void multZ3(const vec2i& a, const vec2i& b, vec2i& out);

  // These routines are primarily used for the arithmetic visibility.
  vec2i moduloGI(const vec2i& a, const vec2i& b);
  vec2i gcdGI(const vec2i& a, const vec2i& b);
  void multGI(const vec2i& a, const vec2i& b, vec2i& out);
  vec2i moduloES(const vec2i& a, const vec2i& b);
  vec2i gcdES(const vec2i& a, const vec2i& b);
  void multES(const vec2i& a, const vec2i& b, vec2i& out);

};

namespace Meta {

  // Some template meta programming (TMP) to determine primality
  // of (constant) integers:
  template <uint p, uint i>
  class isPrime {
  public:
    enum {prim = (p == 2) || ((p % i) &&
                 isPrime<(i > 2 ? p : 0), i - 1>::prim)}; 
  }; 

  template<>
  class isPrime<0, 0> {
  public:
    enum {prim = 1};
  }; 

  template<>
  class isPrime<0, 1> {
  public:
    enum {prim = 1};
  }; 

  template <uint x>
  static inline bool isprime() {
    return (isPrime<x, x-1>::prim == 1);
  }

  template <>
  bool isprime<0>() {
    return false;
  }

  template <>
  bool isprime<1>() {
    return false;
  }

  template <uint u, uint v>
  struct gcd {
    enum { value = gcd<v, u % v>::value };
  };

  template <uint u>
  struct gcd<u, 0> {
    enum { value = u };
  };

  template <>
  struct gcd<0, 0> {
    enum { value = -1 };
  };

  template <uint u, uint v>
  struct coprime {
    enum { value = gcd<u, v>::value == 1 ? 1 : 0 };
  };

  template <uint start, uint end>
  struct eulerPhiLoop {
    enum { value = coprime<start, end>::value +
                   eulerPhiLoop<start + 1, end>::value };
  };

  template <uint end>
  struct eulerPhiLoop<end, end> {
    enum { value = 0 };
  };

  template <uint n>
  struct eulerPhi {
    enum { value = eulerPhiLoop<1, n>::value };
  };

  template <>
  struct eulerPhi<1> {
    enum { value = 1 };
  };

  template <>
  struct eulerPhi<0> {
    enum { value = 1 };
  };

};

namespace Common {

  /* Use the decagon/dodecagon orientation which results from connecting the *
   * ten/twelve roots of unity (false), or the orientation found in the      *
   * book (true), where the right-most edge is aligned with the y-axis.      */
  static const bool windowBookOrientation = true;

  // Replace the window by a circle of the same area.
  static const bool circularWindow = false;

#ifdef COMMON_USE_SSE
  typedef vector<vec2d, AlignmentAllocator<vec2d, 16> > vec2dlist;
  typedef vector<vec4i, AlignmentAllocator<vec4i, 16> > vec4ilist;
#else
  typedef vector<vec2d> vec2dlist;
  typedef vector<vec4i> vec4ilist;
#endif

  typedef vector<vec4s> vec4slist;
  typedef vector<vec2i> vec2ilist;
  typedef vector<vec2iExt> vec2ielist;
  typedef vector<tilingEdge> edgelist;

  typedef vector<double> dlist;
  typedef vector<long double> eflist; // list of extended precision (80-bit IEEE) floats

  void normalize(vec2ielist& in);

  static inline ullong reinterpret_double_to_ullong(const double d) {
    return *(reinterpret_cast<const ullong*>(&d));
  }

  static inline double reinterpret_ullong_to_double(const ullong l) {
    return *(reinterpret_cast<const double*>(&l));
  }

  double checkPosition(const vec2d& a, const vec2d& b, const vec2d& c);
  bool circularCheck(double radSquared, double xSquared);

  static inline uint modulo(int input, uint mod) {
    const int temp = input % int(mod);
    return (temp < 0) ? (int(mod) + temp) : temp;
  }

  /* Only use this if you really have to (very slow!) */
  template <typename T>
  static inline void append_unique(vector<T>& out, const T& elem) {
    if (find(out.begin(), out.end(), elem) == out.end()) {
      out.push_back(elem);
    }
  }

  /* Select vertices based on condition specified in 'S'. */
  template <typename T, typename S>
  static inline void selectVertices(const T& in, T& out) {
    for (typename T::const_iterator i = in.begin(); i != in.end(); ++i) {
      if (S::eval(*i)) out.push_back(*i);
    }
  }

  /* Select the vertices which are contained in a ball of radius R. *
   * The radius is set in radiusSq (=R^2)                           */
  struct RadiusSelector {
    static double radiusSq; /* squared radius */

    static bool eval(const vec4i& x) {
       return (x.paraProjL8().lengthSquared() <= radiusSq);
    }
  };

  /* Select a number of origins from the tiling, based on the following conditions: *
   * Only origins O are selected, where the ball of radius 'sampleRadius' around O  *
   * is still contained in the tiling (which has radius 'tilingRadius').            */
  template <typename T, typename S>
  void selectOrigins(const T& tiling, T& origins, uint samples,
                     double sampleRadius, double tilingRadius) {
    if (sampleRadius >= tilingRadius) {
      cerr << "error: sample radius has to be strictly smaller then tiling radius.\n";
      return;
    }

    const double dist = tilingRadius - sampleRadius;
    T verts;

    for (typename T::const_iterator i = tiling.begin(); i != tiling.end(); ++i) {
      if (S::length(*i) <= dist) verts.push_back(*i);
    }

    const uint numverts = verts.size();
    uint index = 0;

    if (samples > uint(0.1f * float(numverts))) {
      cerr << "error: requested " << samples << " samples, but too few potential origins available.\n";
      return;
    }

    while (samples != 0) {

      /* If n == RAND_MAX, RAND_MAX is already divisible by n   *
       * Else keep searching for an x in a range divisible by n */
      do {
        index = rand();
      } while (numverts < RAND_MAX && index >= RAND_MAX - (RAND_MAX % numverts));

      index %= numverts;

      const typename T::value_type origin(tiling[index]);
      if (find(origins.begin(), origins.end(), origin) != origins.end()) continue;

      origins.push_back(origin);
      --samples;
    }
  }

  // Initialize RNG with external (system) random number
  void srandExt();

  // Create n random double floats in the range [0.0, 1.0]
  // Warning: This call is slow and not suited for large amount of random numbers.
  //          Better use srand and rand in this context!
  void random(uint n, double out[]);

  // Same as random, but creates data in the range [0.0, range]
  void random(uint n, double range, double out[]);

  // VERY naive implementation of a prime number test
  bool isprime(uint x);

  uint eulerPhi(uint i);
  double power(double x, uint y);
  int ipower(int x, uint y);

  static inline void writeRawConsole(const dlist& input) {
    const char* data = reinterpret_cast<const char*>(&(*input.begin()));
    cout.write(data, sizeof(double) * input.size());
  }

  static inline void writeRawConsole(const eflist& input) {
    const char* data = reinterpret_cast<const char*>(&(*input.begin()));
    cout.write(data, sizeof(__float128) * input.size());
  }

  void readRawConsole(dlist& output);
  void writeRawConsole(const vec4ilist& input);
  void readRawConsole(vec4ilist& output);

  /* binning data object:
   * range: binning range, second parameter is ignored in 'tail'-mode
   * step: step size / width of the bins
   * tail: do binning of the 'tail' of the data
   *
   * data: output from binning process
   * catched: number of data points which 'hit' a bin */
  struct BinningData {
    double range[2];
    double step;
    bool tail;

    vector<uint> data;
    uint catched;
  };

  // Compute the minimum/maximum value of the input list
  void minmax(const dlist& input, double& min, double& max);

  /* Compute the index of the 'bin' with the maximum number of entries. *
   * Only works properly if this 'bin' is unique.                       */
  void binmax(const BinningData& input, uint& index);

  // Compute statistics for 'data' and 'bin' and print to console
  void printstats(const dlist& data, const BinningData& bin);

  // Compute binning of 'input' with parameters given in 'output'
  void histogramBinning(const dlist& input, BinningData& output);

  template <typename T>
  void histogramScale(const BinningData& input, vector<T>& output, T scale);

  /* Creates "envelope" data for given histogram input:                 *
   * Can be used for ListPlot to visualize the distributions coming     *
   * from numerical simulations of the radial projection.               */
  void histogramEnvelope(double a, double b, double step, bool stats);

  // Same as histogramEnvelope but processes the "tail" of the data.
  void histoTailEnvelope(double a, double step, bool stats);

  // Histogram envelope routines for large data inputs
  void histogramEnvelopeLD(double a, double b, double step);
  void histoTailEnvelopeLD(double a, double step);

  struct BinningStats {
    double range[2];
    double step;
    bool tail;

    double min, max;
    uint maxbin_index;
    double maxbin_position;
  };

  void binningCopySettings(const BinningStats& input, BinningData& output);

  void histogramStatistics(const dlist& input, BinningStats& output);

  void neighbourDiff(const dlist& input, dlist& output, double& mean);
  void normalizeAngDists(dlist& input, double mean);
  void radialProj(const vec2dlist& input,
                  dlist& output, double& meandist);

  template <typename T>
  static inline uint vectorStats(const vector<T>& input) {
    return ((input.size() * 100) / input.capacity());
  }

  // Common status message we print on screen
  void meanDistanceMessage(uint num, double mean);
};

ostream& operator<<(ostream &os, const vec2i& v);
ostream& operator<<(ostream &os, const vec2d& v);
ostream& operator<<(ostream &os, const vec4i& v);
ostream& operator<<(ostream &os, const vec8s& v);
ostream& operator<<(ostream &os, const vec4s& v);
ostream& operator<<(ostream &os, const vec2s& v);
ostream& operator<<(ostream &os, const vec2iExt& rhs);
ostream& operator<<(ostream &os, const tilingEdge& e);

template <class T, class Allocator>
ostream& operator<<(ostream &os, const vector<T, Allocator>& list) {
  if (!list.empty()) {
    typename vector<T, Allocator>::const_iterator i = list.begin();

    os << '{' << *i;
    ++i;

    while (i != list.end()) {
      os << ',' << *i;
      ++i;
    }

    os << '}';
  } else {
    os << "{}";
  }

  return os;
}

#endif

