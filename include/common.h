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

typedef unsigned int uint;
typedef unsigned short ushort;
typedef unsigned char ubyte;
typedef unsigned long long ullong;

using namespace std;

class vec2d;

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

  // Squaring in Z[Sqrt[2]]
  vec2i squareZ2() const {
    return vec2i(x*x + 2*y*y, 2*x*y);
  }

  // Cubing in Z[Sqrt[2]]
  vec2i cubeZ2() const {
    return vec2i(x*x*x + 6*x*y*y, 3*x*x*y + 2*y*y*y);
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
  
  vec2d transTriToR2() const;

  vec2d transGenericToR2(const vec2d& v) const;

  vec2d minkowskiZ2() const;

  bool coprime() const;

};

class vec2d {
public:
  double x, y;

  vec2d() {}
  vec2d(double a, double b) : x(a), y(b) {}

  vec2d operator*(double scale) const {
    return vec2d(scale * x, scale * y);
  }

  vec2d operator-(const vec2d& v) const {
    return vec2d(x - v.x, y - v.y);
  }

  vec2d operator+(const vec2d& v) const {
    return vec2d(x + v.x, y + v.y);
  }

  vec2d& operator+=(const vec2d& v) {
    x += v.x; y += v.y;
    return *this;
  }

  vec2d operator-() const {
    return vec2d(-x, -y);
  }

  void set(double a, double b) {
    x = a; y = b;
  }

  double lengthSquared() const {
    return x * x + y * y;
  }

  double length() const {
    return sqrt(lengthSquared());
  }

  bool inSectorL8() const {
    return (x >= 0.0 && y >= x);
  }

  bool inUpperHalfplane() const {
    return (y >= 0.0);
  }

  bool inFirstQuadrant() const {
    return (x >= 0.0 && y >= 0.0);
  }

  bool inFirstQuadOpen() const {
    return (x > 0.0 && y > 0.0);
  }

  // Checks for phi(x,y) <= 2*pi/6 (60 degree):
  bool inSectorL3() const {
    return (y/x <= sectorL3);
  }

  // Checks for phi(x,y) <= 2*pi/5 (72 degree):
  bool inSectorL5() const {
    return (y/x <= sectorL5);
  }

  // Checks for phi(x,y) <= 2*pi/6 (60 degree):
  bool inSectorL12() const {
    return (y/x <= sectorL12);
  }

  // Checks for phi(x,y) <= 2*pi/7 (51 degree):
  bool inSectorL7() const {
    return (y/x <= sectorL7);
  }

  vec2d reduceIntoSectorL12() const {
    const double absvals[2] = {abs(x), abs(y)};

    return (absvals[0] >= absvals[1]) ? vec2d(absvals[0], absvals[1]) :
                                        vec2d(absvals[1], absvals[0]);
  }

  double angle() const {
    return atan2(y, x);
  }

  double dot(const vec2d& v) const {
    return x * v.x + y * v.y;
  }

  vec2d normalize() const {
    const double invlen = 1.0 / this->length();
    return vec2d(x * invlen, y * invlen);
  }

  // a = cos(alpha), b = sin(alpha)
  vec2d applyRotation(const double a, const double b) const {
    return vec2d(x*a - y*b, x*b + y*a);
  }

private:
  static const double sectorL3;
  static const double sectorL5;
  static const double sectorL12;

  static const double sectorL7;

};

class vec4i {
private:
  int a[4];

public:
  vec4i() {}
  vec4i(int x0, int x1, int x2, int x3) {
    a[0] = x0; a[1] = x1;
    a[2] = x2; a[3] = x3;
  }

  vec4i operator+(const vec4i& v) const {
    return vec4i(a[0] + v.a[0], a[1] + v.a[1], a[2] + v.a[2], a[3] + v.a[3]);
  }

  vec4i operator-(const vec4i& v) const {
    return vec4i(a[0] - v.a[0], a[1] - v.a[1], a[2] - v.a[2], a[3] - v.a[3]);
  }

  bool operator==(const vec4i& v) const {
    return (a[0] == v.a[0] && a[1] == v.a[1] && a[2] == v.a[2] && a[3] == v.a[3]);
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
    return (a[0] == 0 && a[1] == 0 && a[2] == 0 && a[3] == 0);
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
    static const double v1[4] = {1.0/sqrt(2.0), -0.5, 0.0, 0.5};
    static const double v2[4] = {0.0, -0.5, 1.0/sqrt(2.0), -0.5};

    return vec2d(a[0] * v1[0] + a[1] * v1[1] + a[2] * v1[2] + a[3] * v1[3],
                 a[0] * v2[0] + a[1] * v2[1] + a[2] * v2[2] + a[3] * v2[3]);
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

  // Transform from L5 to direct-sum representation
  vec4s transL5ToDirect() const {
    return vec4s(a[0] - a[2] + a[3], -a[3], a[1] - a[2] + a[3], a[2] - a[3]);
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

  uint gcdZ(uint u, uint v); /* binary gcd implementation (taken from Wikipedia article) */
  int gcdZ(int a, int b);

  static inline bool coprimeZ(int a, int b) {
    return (abs(gcdZ(a, b)) == 1);
  }

  // Used in the octogonal case:
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

};

namespace Meta {

  // Some template meta programming (TMP) to determine primality
  // of (constant) integers:
  template <uint p, uint i>
  class isPrime {
  public:
    enum {prim = (p == 2) || (p % i) &&
                 isPrime<(i > 2 ? p : 0), i - 1>::prim}; 
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

  /* Some constants: */
  static const double pi = atan(1.0) * 4.0;
  static const double eps = numeric_limits<double>::epsilon();

  typedef vector<vec2d> vec2dlist;
  typedef vector<vec4i> vec4ilist;
  typedef vector<vec4s> vec4slist;
  typedef vector<vec2i> vec2ilist;
  typedef vector<tilingEdge> edgelist;

  typedef vector<double> dlist;
  typedef vector<long double> eflist; // list of extended precision (80-bit IEEE) floats

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
  static inline void selectVertices(const vector<T>& in, vector<T>& out) {
    for (typename vector<T>::const_iterator i = in.begin(); i != in.end(); ++i) {
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
  void selectOrigins(const vector<T>& tiling, vector<T>& origins,
                     uint samples, double sampleRadius, double tilingRadius) {
    if (sampleRadius >= tilingRadius) {
      cerr << "error: sample radius has to be strictly smaller then tiling radius.\n";
      return;
    }

    const double dist = tilingRadius - sampleRadius;
    vector<T> verts;

    for (typename vector<T>::const_iterator i = tiling.begin(); i != tiling.end(); ++i) {
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

      const T origin(tiling[index]);
      if (find(origins.begin(), origins.end(), origin) != origins.end()) continue;

      origins.push_back(origin);
      --samples;
    }
  }

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

  uint* histogramBinning(const dlist& input, uint& num_bin, uint& in_bin,
                         const double a, const double b, const double step);
  uint* histoTailBinning(const dlist& input, uint& num_bin, uint& in_bin,
                         const double a, const double step);

  /* Creates "envelope" data for given histogram input:                 *
   * Can be used for ListPlot to visualize the distributions coming     *
   * from numerical simulations of the radial projection.               */
  void histogramEnvelope(const double a, const double b, const double step);

  // Same as histogramEnvelope but processes the "tail" of the data.
  void histoTailEnvelope(const double a, const double step);

  // Histogram envelope routines for large data inputs
  void histogramEnvelopeLD(const double a, const double b, const double step);
  void histoTailEnvelopeLD(const double a, const double step);

  void neighbourDiff(const dlist& input, dlist& output, double& mean);
  void normalizeAngDists(dlist& input, double mean);
  void radialProj(const vec2dlist& input,
                  dlist& output, double& meandist);

  template <typename T>
  static inline uint vectorStats(const vector<T>& input) {
    return ((input.size() * 100) / input.capacity());
  }
};

ostream& operator<<(ostream &os, const vec2i& v);
ostream& operator<<(ostream &os, const vec2d& v);
ostream& operator<<(ostream &os, const vec4i& v);
ostream& operator<<(ostream &os, const vec8s& v);
ostream& operator<<(ostream &os, const vec4s& v);
ostream& operator<<(ostream &os, const vec2s& v);
ostream& operator<<(ostream &os, const tilingEdge& e);

template <class T>
ostream& operator<<(ostream &os, const vector<T>& list) {
  if (!list.empty()) {
    typename vector<T>::const_iterator i = list.begin();

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

