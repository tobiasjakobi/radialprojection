#ifndef _RADIAL_MATH_H_
#define _RADIAL_MATH_H_

#include <cmath>
#include <cstdlib>
#include <cassert>

#include <iostream>
#include <fstream>

#include <iomanip>
#include <vector>
#include <algorithm>
#include <limits>


// Needed for the visibility test routines
#include "pooled_alloc.h"

using namespace std;

class vec2d;

class vec2i {
public:
  int x, y;

  vec2i() {}
  vec2i(int a, int b) : x(a), y(b) {}

  vec2i(const vec2i& v) : x(v.x), y(v.y) {}

  bool operator==(const vec2i& v) const {
    return (x == v.x && y == v.y);
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

  vec2d transHexTo2D() const;

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

  void set(double a, double b) {
    x = a; y = b;
  }

  double lengthSquared() const {
    return x * x + y * y;
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

vec2d vec2i::transHexTo2D() const {
  static const double x1 = sqrt(3.0) * 0.5;

  return vec2d(
    double(x) + double(y) * 0.5,
    double(y) * x1);
}

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

namespace RadialCoprime {

  double round(double number) {
    return number < 0.0 ? ceil(number - 0.5) : floor(number + 0.5);
  }

  // binary gcd implementation (taken from Wikipedia article)
  uint gcdZ(uint u, uint v) {
    uint shift;

    // GCD(0,v) == v; GCD(u,0) == u, GCD(0,0) == 0
    if (u == 0) return v;
    if (v == 0) return u;

    /* Let shift := lg K, where K is the greatest power of 2
       dividing both u and v. */
    for (shift = 0; ((u | v) & 1) == 0; ++shift) {
      u >>= 1;
      v >>= 1;
    }

    while ((u & 1) == 0)
      u >>= 1;

    // From here on, u is always odd.
    do {
      /* remove all factors of 2 in v -- they are not common
         note: v is not zero, so while will terminate */
      while ((v & 1) == 0)  /* Loop X */
        v >>= 1;
 
      /* Now u and v are both odd. Swap if necessary so u <= v,
         then set v = v - u (which is even). For bignums, the
         swapping is just pointer movement, and the subtraction
         can be done in-place. */
      if (u > v) {
        // Swap u and v.
        uint t = v;
        v = u;
        u = t;
      }

      v = v - u; // Here v >= u.
    } while (v != 0);
 
    // restore common factors of 2
    return u << shift;
  }

  int gcdZ(int a, int b) {
    int x = a, y = b, z;

    while (y != 0) {
      z = y;
      y = x % y;
      x = z;
    }

    return x;
  }

  bool coprimeZ(int a, int b) {
    return (abs(gcdZ(a, b)) == 1);
  }

  // Used in the octogonal case:

  vec2i moduloZ2(const vec2i& a, const vec2i& b) {
    const double N = 1.0 / double(b.preNormZ2());

    const int alpha = round(double(a.x * b.x - 2 * a.y * b.y) * N);
    const int beta = round(double(a.y * b.x - a.x * b.y) * N);

    return vec2i(a.x - (alpha * b.x + 2 * beta * b.y),
                 a.y - (alpha * b.y + beta * b.x));
  }

  vec2i gcdZ2(const vec2i& a, const vec2i& b) {
    vec2i x(a);
    vec2i y(b);
    vec2i z;

    while (!y.isZero()) {
      z = y;
      y = moduloZ2(x, y);
      x = z;
    }

    return x;
  }

  bool coprimeZ2(const vec2i& a, const vec2i& b) {
    return (gcdZ2(a, b).normZ2() == 1);
  }

  bool coprimeZ2(const vec4i& v) {
    return (gcdZ2(v.getFirst(), v.getSecond()).normZ2() == 1);
  }

  // a = a.x + a.y * sqrt[2]
  // b = b.x + b.y * sqrt[2]
  // this computes c=a*b in terms of c.x + c.y * sqrt[2]
  // and stores the result in o
  void multZ2(const vec2i& a, const vec2i& b, vec2i& o) {
    o.set(a.x * b.x + 2 * a.y * b.y,
          a.x * b.y + a.y * b.x);
  }

  // Used in the decagonal case:

  vec2i moduloZTau(const vec2i& a, const vec2i& b) {
    const double N = 1.0 / double(b.preNormZTau());

    const int alpha = round(double(a.x*b.x - a.y*b.y + a.x*b.y) * N);
    const int beta = round(double(a.y*b.x - a.x*b.y) * N);

    return vec2i(a.x - (alpha*b.x + beta*b.y),
                 a.y - (alpha*b.y + beta * (b.x + b.y)));
  }

  vec2i gcdZTau(const vec2i& a, const vec2i& b) {
    vec2i x(a);
    vec2i y(b);
    vec2i z;

    while (!y.isZero()) {
      z = y;
      y = moduloZTau(x, y);
      x = z;
    }

    return x;
  }

  vec2i gcdZTau(const vec4s& v) {
    return gcdZTau(vec2i(v[0], v[1]), vec2i(v[2], v[3]));
  }

  bool coprimeZTau(const vec2i& a, const vec2i& b) {
    return (gcdZTau(a, b).normZTau() == 1);
  }

  bool coprimeZTau(const vec4i& v) {
    return (gcdZTau(v.getFirst(), v.getSecond()).normZTau() == 1);
  }

  bool coprimeZTau(const vec4s& v) {
    return (gcdZTau(vec2i(v[0], v[1]), vec2i(v[2], v[3])).normZTau() == 1);
  }

  // a = a.x + a.y * tau
  // b = b.x + b.y * tau (with tau the golden mean)
  // this computes c=a*b in terms of c.x + c.y * tau
  // and stores the result in o
  void multZTau(const vec2i& a, const vec2i& b, vec2i& o) {
    o.set(a.x * b.x + a.y * b.y,
          a.y * (b.x + b.y) + a.x * b.y);
  }

  // Used in the dodecagonal case:
  vec2i moduloZ3(const vec2i& a, const vec2i& b);
  vec2i gcdZ3(const vec2i& a, const vec2i& b);
  bool coprimeZ3(const vec2i& a, const vec2i& b);
  bool coprimeZ3(const vec4i& v);

  vec2i moduloZ3(const vec2i& a, const vec2i& b) {
    const double N = 1.0 / double(b.preNormZ3());

    const int alpha = round(double(a.x * b.x - 3 * a.y * b.y) * N);
    const int beta = round(double(a.y * b.x - a.x * b.y) * N);

    return vec2i(a.x - (alpha * b.x + 3 * beta * b.y),
                 a.y - (alpha * b.y + beta * b.x));
  }

  vec2i gcdZ3(const vec2i& a, const vec2i& b) {
    vec2i x(a);
    vec2i y(b);
    vec2i z;

    while (!y.isZero()) {
      z = y;
      y = moduloZ3(x, y);
      x = z;
    }

    return x;
  }

  vec2i gcdZ3(const vec4i& v) {
    return gcdZ3(v.getFirst(), v.getSecond());
  }

  bool coprimeZ3(const vec2i& a, const vec2i& b) {
    return (gcdZ3(a, b).normZ3() == 1);
  }

  bool coprimeZ3(const vec4i& v) {
    return (gcdZ3(v.getFirst(), v.getSecond()).normZ3() == 1);
  }

  // a = a.x + a.y * sqrt[3]
  // b = b.x + b.y * sqrt[3]
  // this computes c=a*b in terms of c.x + c.y * sqrt[3]
  // and stores the result in o
  void multZ3(const vec2i& a, const vec2i& b, vec2i& o) {
    o.set(a.x * b.x + 3 * a.y * b.y,
          a.x * b.y + a.y * b.x);
  }

};

namespace MetaRadial {

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
  bool isprime() {
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

namespace CommonRadial {

  const double pi = atan(1.0) * 4.0;

  // Create n random double floats in the range [0.0, 1.0]
  // Warning: This call is slow and not suited for large amount of random numbers.
  //          Better use srand and rand in this context!
  void random(uint n, double out[]) {
    if (n == 0) return;

    uint* mem = new uint[n];

    // Use /dev/urandom as RNG (better than the default C/C++ RNG)
    ifstream urandom("/dev/urandom", ios::in|ios::binary);
    urandom.read(reinterpret_cast<char*>(mem), n * sizeof(uint));
    urandom.close();  

    for (uint i = 0; i < n; ++i) {
      out[i] = double(mem[i]) / double(numeric_limits<uint>::max());
    }
    delete [] mem;
  }

  // Same as random, but creates data in the range [0.0, range]
  void random(uint n, double range, double out[]) {
    if (n == 0) return;

    uint* mem = new uint[n];
    const double scaler = range / double(numeric_limits<uint>::max());

    ifstream urandom("/dev/urandom", ios::in|ios::binary);
    urandom.read(reinterpret_cast<char*>(mem), n * sizeof(uint));
    urandom.close();  

    for (uint i = 0; i < n; ++i) {
      out[i] = double(mem[i]) * scaler;
    }
    delete [] mem;
  }

  // VERY naive implementation of a prime number test
  bool isprime(uint x) {
    if (x == 0 || x == 1) return false;

    bool prime = true;
    for (uint i = 2; i < x; ++i) {
      if (x % i == 0) {
        prime = false;
        break;
      }
    }

    return prime;
  }

  uint eulerPhi(uint i) {
    if (i == 1) return 1;

    uint res = i;

    // Check for divisibility by every prime number below the square root.
    // Start with 2.
    if (i % 2 == 0) {
      res -= (res / 2);
      do i /= 2;
      while (i % 2 == 0);
    }

    // Since this doesn't use a list of primes, check every odd number. Ideally, skip past composite numbers.
    for (uint j = 3; j*j <= i; j += 2) {
      if (i % j == 0) {
        res -= (res / j);
        do i /= j;
        while (i % j == 0);
      }
    }

    // If i>1, then it's the last factor at this point.
    if (i > 1) res -= (res / i);

    return res; 
  }

  double power(double x, uint y) {
    if (y == 0) return 1.0;

    double temp = 1.0;
    for (uint i = 0; i < y; ++i) {
      temp *= x;
    }

    return temp;
  }

  int ipower(int x, uint y) {
    if (y == 0) return 1;

    int temp = 1;
    for (uint i = 0; i < y; ++i) {
      temp *= x;
    }

    return temp;
  }

  typedef vector<vec2d> vec2dlist;
  typedef vector<double> dlist;
  typedef vector<long double> eflist; // list of extended precision (80-bit IEEE) floats

  typedef vector<vec4i> vec4ilist;
  typedef vector<vec2i> vec2ilist;

  typedef vector<tilingEdge> edgelist;

  void writeRawConsole(const dlist& input) {
    const char* data = reinterpret_cast<const char*>(&(*input.begin()));
    cout.write(data, sizeof(double) * input.size());
  }

  void writeRawConsole(const eflist& input) {
    const char* data = reinterpret_cast<const char*>(&(*input.begin()));
    cout.write(data, sizeof(__float128) * input.size());
  }

  void readRawConsole(dlist& output) {
    double temp;
    while (true) {
      cin.read(reinterpret_cast<char*>(&temp), sizeof(double));
      if (cin.eof()) break;

      output.push_back(temp);
    }
  }

  void writeRawConsole(const vec4ilist& input) {
    unsigned out;
    const char* data;

    /* write signature first */

    out = sizeof(int); /* element size */
    cout.write(reinterpret_cast<const char*>(&out), sizeof(unsigned) * 1);

    out = 4; /* number of elements per entry */
    cout.write(reinterpret_cast<const char*>(&out), sizeof(unsigned) * 1);

    out = input.size();
    cout.write(reinterpret_cast<const char*>(&out), sizeof(unsigned) * 1);

    data = reinterpret_cast<const char*>(&(*input.begin()));
    cout.write(data, sizeof(vec4i) * input.size());
  }

  void readRawConsole(vec4ilist& output) {
    unsigned in;
    vec4i data;

    cin.read(reinterpret_cast<char*>(&in), sizeof(unsigned));
    if (cin.eof() && cin.fail()) goto readfail;
    if (in != sizeof(int)) goto signfail;

    cin.read(reinterpret_cast<char*>(&in), sizeof(unsigned));
    if (cin.eof() && cin.fail()) goto readfail;
    if (in != 4) goto signfail;

    cin.read(reinterpret_cast<char*>(&in), sizeof(unsigned));
    if (cin.eof() && cin.fail()) goto readfail;
    
    output.reserve(output.size() + in);

    while (true) {
      cin.read(reinterpret_cast<char*>(&data), sizeof(vec4i));
      if (cin.eof()) break;

      output.push_back(data);
    }

    return;
  signfail:
    cerr << "error: verifying signature failed\n";
    return ;

  readfail:
    cerr << "error: reading signature failed\n";
    return ;
  }

  uint* histogramBinning(const dlist& input, uint& num_bin, uint& in_bin,
                         const double a, const double b, const double step) {
    num_bin = uint((b - a) / step);
    in_bin = 0;

    uint* bins = new uint[num_bin];
    for (uint i = 0; i < num_bin; ++i)
      bins[i] = 0;

    for (dlist::const_iterator j = input.begin(); j != input.end(); ++j) {
      const double cur = *j;

      if (cur < a || cur >= b) continue;

      ++bins[uint((cur - a) / step)];
      ++in_bin;
    }

    return bins;
  }

  uint* histoTailBinning(const dlist& input, uint& num_bin, uint& in_bin,
                         const double a, const double step) {
    if (input.empty())
      return NULL;

    // find maximum value of input list
    dlist::const_iterator j = input.begin();
    double input_max = *j; ++j;
    for (; j != input.end(); ++j) {
      if (*j > input_max)
        input_max = *j;
    }

    // check if we don't overflow the (maximum) bin count
    if ((input_max + 1.0 - a) / step > double(numeric_limits<uint>::max())) {
      return NULL;
    }

    num_bin = uint((input_max + 1.0 - a) / step);
    in_bin = 0;

    uint* bins = new uint[num_bin];
    for (uint i = 0; i < num_bin; ++i)
      bins[i] = 0;

    for (j = input.begin(); j != input.end(); ++j) {
      const double cur = *j;

      if (cur < a) continue;

      ++bins[uint((cur - a) / step)];
      ++in_bin;
    }

    return bins;
  }

  /* Creates "envelope" data for given histogram input:                 *
   * Can be used for ListPlot to visualize the distributions coming     *
   * from numerical simulations of the radial projection.               */
  void histogramEnvelope(const double a, const double b, const double step) {
    // fetch input data from console
    dlist inputData;
    readRawConsole(inputData);

    // compute histogram binning
    uint num_bin, in_bin;
    uint* binData = histogramBinning(inputData, num_bin, in_bin, a, b, step);
    dlist envelopeData;

    cerr << "Computing histogram with " << num_bin << " bins (interval = ["
         << a << ',' << b << "); step width = " << step << ")\n";
    cerr << "statistics: " << in_bin << " data points (from " << inputData.size()
         << ") fall into the binning area\n";

    const double scaler = 1.0 / (double(inputData.size()) * step);
    for (uint i = 0; i < num_bin; ++i) {
      envelopeData.push_back(double(binData[i]) * scaler);
    }
    delete [] binData;

    // output to console
    writeRawConsole(envelopeData);
  }

  // Same as histogramEnvelope but processes the "tail" of the data.
  void histoTailEnvelope(const double a, const double step) {
    dlist inputData;
    readRawConsole(inputData);

    uint num_bin, in_bin;
    uint* binData = histoTailBinning(inputData, num_bin, in_bin, a, step);

    // Since the tail is usually quite long, we need extended precision here
    eflist envelopeData;

    if (binData == NULL) {
      cerr << "Binning failed due to range/precision issues.\n";
      return;
    }

    cerr << "Computing histogram with " << num_bin << " bins (tail beginning from "
         << a << "; step width = " << step << ")\n";
    cerr << "statistics: " << in_bin << " data points (from " << inputData.size()
         << ") fall into the binning area\n";

    const long double scaler = 1.0L / (long double)(double(inputData.size()) * step);
    for (uint i = 0; i < num_bin; ++i) {
      envelopeData.push_back((long double)(binData[i]) * scaler);
    }
    delete [] binData;

    writeRawConsole(envelopeData);
  }

  void neighbourDiff(const dlist& input, dlist& output, double& mean) {
    assert(input.size() != 0 && input.size() != 1);

    dlist::const_iterator j = input.begin();
    mean = *j;

    dlist::const_iterator k = j;
    ++k;
    while (k != input.end()) {
      output.push_back(*k - *j);
      ++j;
      ++k;
    }

    --j;
    mean = (*j - mean) / double(input.size() - 1);
  }

  void normalizeAngDists(dlist& input, double mean) {
    const double normalizer = 1.0 / mean;

    for (dlist::iterator i = input.begin(); i != input.end(); ++i) {
      (*i) *= normalizer;
    }
  }

  void radialProj(const vec2dlist& input,
                  dlist& output, double& meandist) {
    output.clear();
    output.reserve(input.size());

    dlist angles;
    angles.reserve(input.size());

    for (vec2dlist::const_iterator i = input.begin(); i != input.end(); ++i)
      angles.push_back(i->angle());

    sort(angles.begin(), angles.end());
    neighbourDiff(angles, output, meandist);
    normalizeAngDists(output, meandist);
  }

  template <typename T>
  uint vectorStats(const vector<T>& input) {
    return ((input.size() * 100) / input.capacity());
  }

};

namespace VisTest {

  template<typename T>
  class visibleList {
  private:
    typedef typename T::invectype _invectype;

    typedef typename PooledList<_invectype>::Type _vlist;
    typedef typename _vlist::iterator _iter;
    typedef typename _vlist::const_iterator _citer;

    OneTimePool* mempool;

    _vlist* l; // internal point list (double-linked)
    _iter c; // last insert position

  public:

    visibleList() : mempool(0), l(0) {}

    visibleList(const visibleList& vl) : mempool(0), l(0) {
      if (vl.mempool != 0) {
        mempool = new OneTimePool(vl.mempool->getByteSize());
      }

      if (vl.l != 0) {
        l = new PooledList<vec4s>::Type(mempool);
        c = l->end();
      }
    }

    ~visibleList() {
      cerr << "info: before destructing the pool used " << mempool->getBytesUsed()
           << " out of " << mempool->getByteSize() << " bytes ("
           << (double(mempool->getBytesUsed()) / double(mempool->getByteSize()) * 100.0)
           << "%).\n";

      delete l;
      delete mempool;
    }

    void reserve(uint elements) {
      if (mempool != 0) return;

      const size_t bytes = elements * PooledList<_invectype>::NodeByteSize;
      cerr << "info: reserving memory pool of " << bytes << " bytes.\n";

      mempool = new OneTimePool(bytes);
    }

    void init() {
      if (l != 0) return;

      l = new typename PooledList<_invectype>::Type(mempool);
      c = l->end();
    }

    void insertSorted(const _invectype& v) {
      const double v_a = T::angle(v);
      const double c_a = T::angle(*c);

      if (l->empty()) {
        insert(v, l->end());
        return;
      }

      bool found = false;
      _iter p = c;

      if (v_a == c_a) {
        found = isearch(v);
      } else {
        if (v_a > c_a) {
          found = fsearch(v, p);
        } else {
          // v_a < c_a
          found = bsearch(v, p);
        }
      }

      if (!found) {
        insert(v, p);
      } else {
        c = p;
      }
    }

    // collect nodes that are "near" i
    void collectNearNodes(const _iter& i, vector<_iter>& list) {
      const double i_a = T::angle(*i);

      _iter j;

      // collect in forward direction
      j = i;
      ++j;
      while (j != l->end()) {
        if (T::angle(*j) > i_a + T::epsilon)
          break;

        list.push_back(j);
        ++j;
      }

      // collect in backward direction
      j = i;
      while (j != l->begin()) {
        --j;

        if (T::angle(*j) < i_a - T::epsilon)
          break;

        list.push_back(j);
      }
    }

    /* This DOESN'T remove invisible points but removes all but one vertex  *
     * on a ray (through zero, the origin) containing multiple vertices     *
     * from the patch. So this produces a vertex set where each vertex      *
     * is visible from the origin.                                          *
     * This suffices for the radial projection approach!                    */
    void removeInvisibleFast() {
      cerr << "info: computing (incorrect) visibility for "
           << l->size() << " vertices." << endl;

      vector<_iter> nodes;

      for (_iter i = l->begin(); i != l->end(); ++i) {
        collectNearNodes(i, nodes);

        for (typename vector<_iter>::iterator j = nodes.begin(); j != nodes.end(); ++j) {
          const _iter k = *j;
          if (T::rayTest(*i, *k)) l->erase(k);
        }

        nodes.clear();
      }

      cerr << "info: " << l->size() << " vertices are visible." << endl;
    }

    // This method really does what the name implies!
    // It only removes but the vertex being nearest to the origin.
    void removeInvisibleProper() {
      cerr << "info: computing (proper) visibility for "
           << l->size() << " vertices." << endl;

      vector<_iter> nodes;

      for (_iter i = l->begin(); i != l->end(); ++i) {
        collectNearNodes(i, nodes);

        for (typename vector<_iter>::iterator j = nodes.begin(); j != nodes.end(); ++j) {
          const _iter k = *j;
          if (!T::rayTest(*i, *k)) continue;

          if (T::toR2(*k).lengthSquared() >= T::toR2(*i).lengthSquared())
            l->erase(k);
        }

        nodes.clear();
      }

      cerr << "info: " << l->size() << " vertices are visible." << endl;
    }

    void toR2(CommonRadial::vec2dlist& output) const {
      for (_citer i = l->begin(); i != l->end(); ++i) {
        output.push_back(T::toR2(*i));
      }
    }

    void dump(vector<_invectype>& output) const {
      for (_citer i = l->begin(); i != l->end(); ++i) {
        output.push_back(*i);
      }
    }

    uint size() const {
      return l->size();
    }

    bool isInitialized() const {
      return l != 0;
    }

    bool isOrdered() const {
      if (l->size() == 0 || l->size() == 1) {
        return true;
      }

      _citer i = l->begin();

      bool ordered = true;
      double a = T::angle(*i);

      ++i;
      while (i != l->end()) {
        const double b = T::angle(*i);

        if (a > b) {
          ordered = false;
          break;
        }

        a = b;
        ++i;
      }

      return ordered;
    }

    void print() const {
      if (l->size() == 0) return;

      _citer i = l->begin();
      cerr << T::angle(*i);;
      ++i;

      while (i != l->end()) {
        cerr << ' ' << T::angle(*i);
        ++i;
      }
      cerr << endl;
    }

  private:

    bool isearch(const _invectype& v) {
      const double v_a = T::angle(v);
      _iter i;

      i = c;
      while (i != l->end()) {
        if (T::angle(*i) != v_a)
          break;

        if (v == *i)
          return true;

        ++i;
      }

      i = c;
      while (i != l->begin()) {
        --i;

        if (T::angle(*i) != v_a)
          break;

        if (v == *i)
          return true;
      }

      return false;
    }

    // search in forward direction for the entry v:
    // returns true if entry found
    // returns false otherwise together with an insert position
    bool fsearch(const _invectype& v, _iter& inspos) {
      const double v_a = T::angle(v);

      _iter i = c;

      while (i != l->end()) {
        if (v == *i) {
          return true;
        }

        if (T::angle(*i) > v_a) {
          inspos = i;
          return false;
        }

        ++i;
      }

      // reached the end of the list without finding anything
      // and v_a is greater than all entries in the list

      inspos = l->end();

      return false;
    }

    bool bsearch(const _invectype& v, _iter& inspos) {
      const double v_a = T::angle(v);

      _iter i = c;
      ++i;

      while (i != l->begin()) {
        --i;

        if (v == *i) {
          return true;
        }

        if (T::angle(*i) < v_a) {
          inspos = i;
          ++inspos;
          return false;
        }
      }

      // reached the start of the list without finding anything
      // and v_a is smaller than all entries in the list

      inspos = l->begin();

      return false;
    }

    void insert(const _invectype& v, const _iter& pos) {
      c = l->insert(pos, v);
    }

  };

};

ostream& operator<<(ostream &os, const vec2i& v)
{
  os << '{' << v.x << ',' << v.y << '}';
  return os;
}

ostream& operator<<(ostream &os, const vec2d& v)
{
  os << '{' << fixed << setprecision(3)
     << v.x << ',' << v.y << '}';
  return os;
}

ostream& operator<<(ostream &os, const vec4i& v)
{
  os << '{' << v[0] << ',' << v[1] << ',' << v[2] << ',' << v[3] << '}';
  return os;
}

ostream& operator<<(ostream &os, const vec8s& v)
{
  os << '{' << int(v[0]);

  for (uint i = 1; i < 8; ++i) {
    os << ',' << int(v[i]);
  }

  os << '}';

  return os;
}

ostream& operator<<(ostream &os, const vec4s& v)
{
  os << '{' << int(v[0]) << ',' << int(v[1])
     << ',' << int(v[2]) << ',' << int(v[3])
     << '}';

  return os;
}

ostream& operator<<(ostream &os, const vec2s& v)
{
  os << '{' << int(v[0])
     << ',' << int(v[1]) << '}';

  return os;
}

ostream& operator<<(ostream &os, const tilingEdge& e)
{
  os << '{' << e[0]
     << ',' << e[1] << '}';

  return os;
}

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

// This is tan(2*pi/6).
const double vec2d::sectorL3 = sqrt(3.0);
// This is tan(2*pi/5).
const double vec2d::sectorL5 = sqrt(5.0 + 2.0*sqrt(5.0));
// This is tan(2*pi/6).
const double vec2d::sectorL12 = sqrt(3.0);

// Used for the higher-order cyclotomic cases:
const double vec2d::sectorL7 = tan(2.0 * CommonRadial::pi / 7.0);

#endif

