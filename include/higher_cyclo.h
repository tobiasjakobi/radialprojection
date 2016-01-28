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

#ifndef _HIGHER_CYCLO_H_
#define _HIGHER_CYCLO_H_

#include "common.h"
#include "visibility.h"

// Used for the internal space of the heptagonal tiling:
class vec4d {
private:
  double a[4];

public:
  vec4d() {}
  vec4d(double x, double y, double z, double w) {
    a[0] = x; a[1] = y;
    a[2] = z; a[3] = w;
  }

  vec4d(const vec2d& v1, const vec2d& v2) {
    a[0] = v1[0]; a[1] = v1[1];
    a[2] = v2[0]; a[3] = v2[1];
  }

  vec4d operator*(double scale) const {
    return vec4d(scale * a[0], scale * a[1],
                 scale * a[2], scale * a[3]);
  }

  vec4d operator-(const vec4d& v) const {
    return vec4d(a[0] - v.a[0], a[1] - v.a[1],
                 a[2] - v.a[2], a[3] - v.a[3]);
  }

  vec4d operator+(const vec4d& v) const {
    return vec4d(a[0] + v.a[0], a[1] + v.a[1],
                 a[2] + v.a[2], a[3] + v.a[3]);
  }

  vec4d& operator+=(const vec4d& v) {
    a[0] += v.a[0]; a[1] += v.a[1];
    a[2] += v.a[2]; a[3] += v.a[3];
    return *this;
  }

  double operator[](uint i) const {
    assert(i < 4);
    return a[i];
  }

  double& operator[](uint i) {
    assert(i < 4);
    return a[i];
  }

  void set(double x, double y, double z, double w) {
    a[0] = x; a[1] = y;
    a[2] = z; a[3] = w;
  }

  double lengthSquared() const {
    return a[0] * a[0] + a[1] * a[1] + a[2] * a[2] + a[3] * a[3];
  }

  double length() const {
    return sqrt(this->lengthSquared());
  }

  double dot(const vec4d& v) const {
    return a[0] * v.a[0] + a[1] * v.a[1] +
           a[2] * v.a[2] + a[3] * v.a[3];
  }

};

// Used for the internal space of the elevenfold tiling:
class vec8d {
private:
  double a[8];

public:
  vec8d() {}
  vec8d(double x1, double x2, double y1, double y2,
        double z1, double z2, double w1, double w2) {
    a[0] = x1; a[1] = x2; a[2] = y1; a[3] = y2;
    a[4] = z1; a[5] = z2; a[6] = w1; a[7] = w2;
  }

  vec8d(const vec2d* array) {
    for (uint i = 0; i < 4; ++i) {
      a[i*2 + 0] = array[i][0];
      a[i*2 + 1] = array[i][1];
    }
  }

  vec8d(const double* array) {
    for (uint i = 0; i < 8; ++i)
      a[i] = array[i];
  }

  vec8d operator*(double scale) const {
    return vec8d(scale * a[0], scale * a[1], scale * a[2], scale * a[3],
                 scale * a[4], scale * a[5], scale * a[6], scale * a[7]);
  }

  vec8d operator-(const vec8d& v) const {
    return vec8d(a[0] - v.a[0], a[1] - v.a[1], a[2] - v.a[2], a[3] - v.a[3],
                 a[4] - v.a[4], a[5] - v.a[5], a[6] - v.a[6], a[7] - v.a[7]);
  }

  vec8d operator+(const vec8d& v) const {
    return vec8d(a[0] + v.a[0], a[1] + v.a[1], a[2] + v.a[2], a[3] + v.a[3],
                 a[4] + v.a[4], a[5] + v.a[5], a[6] + v.a[6], a[7] + v.a[7]);
  }

  vec8d& operator+=(const vec8d& v) {
    for (uint i = 0; i < 8; ++i)
      a[i] += v.a[i];
    return *this;
  }

  double operator[](uint i) const {
    assert(i < 8);
    return a[i];
  }

  double& operator[](uint i) {
    assert(i < 8);
    return a[i];
  }

  void set(double x1, double x2, double y1, double y2,
           double z1, double z2, double w1, double w2) {
    a[0] = x1; a[1] = x2;
    a[2] = y1; a[3] = y2;
    a[4] = z1; a[5] = z2;
    a[6] = w1; a[7] = w2;
  }

  double lengthSquared() const {
    double ret = 0.0;
    for (uint i = 0; i < 8; ++i)
      ret += a[i] * a[i];
    return ret;
  }

  double length() const {
    return sqrt(this->lengthSquared());
  }

  double dot(const vec8d& v) const {
    double ret = 0.0;
    for (uint i = 0; i < 8; ++i)
      ret += a[i] * v.a[i];
    return ret;
  }

};

// Used for creating the tiling vertices of the heptagonal tiling:
class vec6s {
private:
  short a[6];

public:
  vec6s() {}
  vec6s(short x0, short x1, short x2,
        short x3, short x4, short x5) {
    a[0] = x0; a[1] = x1; a[2] = x2;
    a[3] = x3; a[4] = x4; a[5] = x5;
  }

  vec6s operator+(const vec6s& v) const {
    return vec6s(a[0] + v.a[0], a[1] + v.a[1], a[2] + v.a[2],
                 a[3] + v.a[3], a[4] + v.a[4], a[5] + v.a[5]);
  }

  vec6s operator-(const vec6s& v) const {
    return vec6s(a[0] - v.a[0], a[1] - v.a[1], a[2] - v.a[2],
                 a[3] - v.a[3], a[4] - v.a[4], a[5] - v.a[5]);
  }

  /* Lexicographic ordering, this is needed to use sorting *
   * algorithms of STL containers.                         */
  bool operator<(const vec6s& v) const {
    for (uint i = 0; i < 5; ++i) {
      if (a[i] < v[i]) return true;
      if (a[i] != v[i]) return false;
    }

    return (a[5] < v[5]);
  }

  bool operator==(const vec6s& v) const {
    return (a[0] == v.a[0] && a[1] == v.a[1] && a[2] == v.a[2] &&
            a[3] == v.a[3] && a[4] == v.a[4] && a[5] == v.a[5]);
  }

  bool operator!=(const vec6s& v) const {
    return !(*this == v);
  }

  short operator[](uint i) const {
    assert(i < 6);
    return a[i];
  }

  short& operator[](uint i) {
    assert(i < 6);
    return a[i];
  }

  void zero() {
    a[0] = a[1] = a[2] = a[3] = a[4] = a[5] = 0;
  }

  bool isZero() const {
    return (a[0] == 0 && a[1] == 0 && a[2] == 0 &&
            a[3] == 0 && a[4] == 0 && a[5] == 0);
  }

  bool isFirstZero() const {
    return (a[0] == 0 && a[1] == 0 && a[2] == 0);
  }

  bool isSecondZero() const {
    return (a[3] == 0 && a[4] == 0 && a[5] == 0);
  }

  void interleave(vec6s& v) {
    short temp[3];

    temp[0] = a[0]; temp[1] = a[1]; temp[2] = a[2];
    a[0] = v.a[0]; a[1] = v.a[1]; a[2] = v.a[2];
    v.a[0] = temp[0]; v.a[1] = temp[1]; v.a[2] = temp[2];
  }

  vec6s step(uint direction, bool invert) const {
    assert(direction < 7);

    vec6s copy(*this);
    const short x = invert ? -1 : 1;

    if (direction < 6)
      copy.a[direction] += x;
    else {
      for (uint i = 0; i < 6; ++i)
        copy.a[i] -= x;
    }

    return copy;
  }

  vec6s toDirectL7() const {
    return vec6s(a[0] - a[2] + a[4] - a[5], -a[3], -a[4] + a[5],
                 a[1] - a[3] + a[4], a[2] - a[5], a[3] - a[4]);
  }

  vec6s conjugateL7(uint i) const {
    assert(i < 2);

    return (i == 0) ? vec6s(a[0] - a[3], a[4] - a[3], a[1] - a[3],
                            a[5] - a[3], a[2] - a[3], -a[3]) :
                      vec6s(a[0] - a[2], a[5] - a[2], a[3] - a[2],
                            a[1] - a[2], -a[2], a[4] - a[2]);
  }

  vec2d toPhysicalL7() const {
    static const double u = Constants::pi / 14.0;

    static const double v1[6] = {1.0, sin(3.0*u), -sin(u),
                                 -cos(2.0*u), -cos(2.0*u), -sin(u)};
    static const double v2[6] = {0.0, cos(3.0*u),  cos(u),
                                  sin(2.0*u), -sin(2.0*u), -cos(u)};

    return vec2d(
      a[0] * v1[0] + a[1] * v1[1] + a[2] * v1[2] +
      a[3] * v1[3] + a[4] * v1[4] + a[5] * v1[5],
      a[0] * v2[0] + a[1] * v2[1] + a[2] * v2[2] +
      a[3] * v2[3] + a[4] * v2[4] + a[5] * v2[5]);
  }

  vec4d toInternalL7() const {
    const vec2d w1(this->conjugateL7(0).toPhysicalL7());
    const vec2d w2(this->conjugateL7(1).toPhysicalL7());

    return vec4d(w1, w2);
  }

};

// Used for creating the tiling vertices of the elevenfold tiling:
class vec10s {
  private:
  short a[10];

public:
  vec10s() {}
  vec10s(short x0, short x1, short x2, short x3, short x4,
         short x5, short x6, short x7, short x8, short x9) {
    a[0] = x0; a[1] = x1; a[2] = x2; a[3] = x3; a[4] = x4;
    a[5] = x5; a[6] = x6; a[7] = x7; a[8] = x8; a[9] = x9;
  }

  vec10s(const short* array) {
    for (uint i = 0; i < 10; ++i)
      a[i] = array[i];
  }

  vec10s operator+(const vec10s& v) const {
    short res[10];
    for (uint i = 0; i < 10; ++i)
      res[i] = a[i] + v.a[i];
    return vec10s(res);
  }

  vec10s operator-(const vec10s& v) const {
    short res[10];
    for (uint i = 0; i < 10; ++i)
      res[i] = a[i] - v.a[i];
    return vec10s(res);
  }

  /* Lexicographic ordering, this is needed to use sorting *
   * algorithms of STL containers.                         */
  bool operator<(const vec10s& v) const {
    for (uint i = 0; i < 9; ++i) {
      if (a[i] < v[i]) return true;
      if (a[i] != v[i]) return false;
    }

    return (a[9] < v[9]);
  }

  bool operator==(const vec10s& v) const {
    for (uint i = 0; i < 10; ++i) {
      if (a[i] != v.a[i]) return false;
    }
    return true;
  }

  bool operator!=(const vec10s& v) const {
    return !(*this == v);
  }

  short operator[](uint i) const {
    assert(i < 10);
    return a[i];
  }

  short& operator[](uint i) {
    assert(i < 10);
    return a[i];
  }

  void zero() {
    for (uint i = 0; i < 10; ++i)
      a[i] = 0;
  }

  bool isZero() const {
    for (uint i = 0; i < 10; ++i) {
      if (a[i] != 0) return false;
    }
    return true;
  }

  bool isFirstZero() const {
    for (uint i = 0; i < 5; ++i) {
      if (a[i] != 0) return false;
    }
    return true;
  }

  bool isSecondZero() const {
    for (uint i = 5; i < 10; ++i) {
      if (a[i] != 0) return false;
    }
    return true;
  }

  void interleave(vec10s& v) {
    short temp;

    for (uint i = 0; i < 5; ++i) {
      temp = a[i];
      a[i] = v.a[i];
      v.a[i] = temp;
    }
  }

  vec10s step(uint direction, bool invert) const {
    assert(direction < 11);

    vec10s copy(*this);
    const short x = invert ? -1 : 1;

    if (direction < 10)
      copy.a[direction] += x;
    else {
      for (uint i = 0; i < 10; ++i)
        copy.a[i] -= x;
    }

    return copy;
  }

  vec10s toDirectL11() const {
    return vec10s(a[0] - a[2] + a[4] - a[6] + a[7] - a[9],
                 -a[3] + 2 * a[5] - 2 * a[8],
                 -a[4] + 3 * a[6] - 3 * a[7] + a[9],
                 -a[5] + a[8],
                 -a[6] + a[7],
                  a[1] - a[3] + a[5] - a[6] + a[8],
                  a[2] - 2 * a[4] + 2 * a[7] - a[9],
                  a[3] - 3 * a[5] + 3 * a[6] - a[8], a[4] - a[7],
                  a[5] - a[6]);
  }

  vec10s conjugateL11(uint i) const {
    assert(i < 4);

    switch (i) {
      case 0: // first conjugation sigma_2 (mapping zeta to zeta^2)
        return vec10s(a[0] - a[5], -a[5] + a[6], a[1] - a[5], -a[5] + a[7],
                      a[2] - a[5], -a[5] + a[8], a[3] - a[5], -a[5] + a[9],
                      a[4] - a[5], -a[5]);

      case 1: // conjugation sigma_3 (mapping zeta to zeta^3)
        return vec10s(a[0] - a[7],  a[4] - a[7], -a[7] + a[8], a[1] - a[7],
                      a[5] - a[7], -a[7] + a[9],  a[2] - a[7], a[6] - a[7],
                     -a[7], a[3] - a[7]);

      case 2: // conjugation sigma_4 (mapping zeta to zeta^4)
        return vec10s(a[0] - a[8], a[3] - a[8], a[6] - a[8], -a[8] + a[9],
                      a[1] - a[8], a[4] - a[8], a[7] - a[8], -a[8],
                      a[2] - a[8], a[5] - a[8]);

      case 3: // conjugation sigma_5 (mapping zeta to zeta^4)
        return vec10s(a[0] - a[2], -a[2] + a[9], -a[2] + a[7], -a[2] + a[5],
                     -a[2] + a[3],  a[1] - a[2], -a[2], -a[2] + a[8],
                     -a[2] + a[6], -a[2] + a[4]);

      default:
        assert(false);
        return vec10s();
    }
  }

  vec2d toPhysicalL11() const {
    static const double u = Constants::pi / 22.0;

    static const double v1[10] = {1.0, cos(4.0*u), sin(3.0*u), -sin(u),
      -sin(5.0*u), -cos(u), -cos(u), -sin(5.0*u), -sin(u),  sin(3.0*u)};
    static const double v2[10] = {0.0, sin(4.0*u), cos(3.0*u),  cos(u),
       cos(5.0*u),  sin(u), -sin(u), -cos(5.0*u), -cos(u), -cos(3.0*u)};

    double ret[2] = {0.0, 0.0};

    for (uint i = 0; i < 10; ++i) {
      ret[0] += a[i] * v1[i];
      ret[1] += a[i] * v2[i];
    }

    return vec2d(ret[0], ret[1]);
  }

  vec8d toInternalL11() const {
    vec2d conj[4];

    for (uint i = 0; i < 4; ++i)
      conj[i] = this->conjugateL11(i).toPhysicalL11();

    return vec8d(conj);
  }
};

namespace Common {

  typedef vector<vec4d> vec4dlist;
  typedef vector<vec6s> vec6slist;

  typedef vector<vec8d> vec8dlist;
  typedef vector<vec10s> vec10slist;

};

namespace Coprime {

  /* Multiply in Z[heptaLambda] (z = x1 * x2):            *
   * First 3 components of 'in' is the first factor x1,   *
   * the last 3 components the other factor x2.           *
   * Result is stored in the first 3 components of 'out'. */
  void multL7(const vec6s& in, vec6s& out) {
    const short u0 = in[1] * in[3+2] + in[2] * in[3+1];
    const short u1 = in[2] * in[3+2];

    out[0] = in[0] * in[3+0] + u0 - u1;
    out[1] = in[0] * in[3+1] + in[1] * in[3+0] + 2*u0 - u1;
    out[2] = in[0] * in[3+2] + in[1] * in[3+1] + in[2] * in[3+0] - u0 + 3*u1;
  }

  void multL11(const vec10s& in, vec10s& out) {
    const short temp[] = {
      in[0] * in[5+0], // i+j=0
      in[0] * in[5+1] + in[1] * in[5+0],
      in[0] * in[5+2] + in[2] * in[5+0] + in[1] * in[5+1],
      in[0] * in[5+3] + in[3] * in[5+0] + in[1] * in[5+2] + in[2] * in[5+1],
      in[0] * in[5+4] + in[4] * in[5+0] + in[1] * in[5+3] + in[3] * in[5+1] +
        in[2] * in[5+2], // i+j=4
      in[1] * in[5+4] + in[4] * in[5+1] + in[2] * in[5+3] + in[3] * in[5+2],
      in[2] * in[5+4] + in[4] * in[5+2] + in[3] * in[5+3], // i+j=6
      in[3] * in[5+4] + in[4] * in[5+3],
      in[4] * in[5+4] // i+j=8
    };

    out[0] = temp[0] - 1*temp[5] + 1*temp[6] -  5*temp[7] +  6*temp[8];
    out[1] = temp[1] - 3*temp[5] + 2*temp[6] - 14*temp[7] + 13*temp[8];
    out[2] = temp[2] + 3*temp[5] - 6*temp[6] + 17*temp[7] - 32*temp[8];
    out[3] = temp[3] + 4*temp[5] - 1*temp[6] + 14*temp[7] -  7*temp[8];
    out[4] = temp[4] - 1*temp[5] + 5*temp[6] -  6*temp[7] + 20*temp[8];
  }

};

namespace Heptagonal {

  struct VisOp {
    typedef Common::vec6slist list_type;
    static const double epsilon;

    static double angle(const vec6s& a) {
      return a.toPhysicalL7().angle();
    }

    static vec2d toR2(const vec6s& a) {
      return a.toPhysicalL7();
    }

    static bool rayTest(const vec6s& a, const vec6s& b) {
      // convert to direct sum representation
      vec6s pa(a.toDirectL7());
      vec6s pb(b.toDirectL7());

      if (pa.isFirstZero()) {
        return pb.isFirstZero();
      }

      if (pb.isFirstZero()) {
        return pa.isFirstZero();
      }

      if (pa.isSecondZero()) {
        return pb.isSecondZero();
      }

      if (pb.isSecondZero()) {
        return pa.isSecondZero();
      }

      vec6s c, d;
      c.zero(); d.zero();

      pa.interleave(pb);
      Coprime::multL7(pa, c);
      Coprime::multL7(pb, d);

      return (c == d);
    }
  };

  typedef VisTest::VisibleList<VisOp> VisList;

};

namespace Elevenfold {

  struct VisOp {
    typedef Common::vec10slist list_type;
    static const double epsilon;

    static double angle(const vec10s& a) {
      return a.toPhysicalL11().angle();
    }

    static vec2d toR2(const vec10s& a) {
      return a.toPhysicalL11();
    }

    static bool rayTest(const vec10s& a, const vec10s& b) {
      // convert to direct sum representation
      vec10s pa(a.toDirectL11());
      vec10s pb(b.toDirectL11());

      if (pa.isFirstZero()) {
        return pb.isFirstZero();
      }

      if (pb.isFirstZero()) {
        return pa.isFirstZero();
      }

      if (pa.isSecondZero()) {
        return pb.isSecondZero();
      }

      if (pb.isSecondZero()) {
        return pa.isSecondZero();
      }

      vec10s c, d;
      c.zero(); d.zero();

      pa.interleave(pb);
      Coprime::multL11(pa, c);
      Coprime::multL11(pb, d);

      return (c == d);
    }
  };

  typedef VisTest::VisibleList<VisOp> VisList;

};

ostream& operator<<(ostream &os, const vec4d& v);
ostream& operator<<(ostream &os, const vec8d& v);

ostream& operator<<(ostream &os, const vec6s& v)
{
  os << '{' << int(v[0]);

  for (uint i = 1; i < 6; ++i) {
    os << ',' << int(v[i]);
  }

  os << '}';

  return os;
}

ostream& operator<<(ostream &os, const vec10s& v)
{
  os << '{' << int(v[0]);

  for (uint i = 1; i < 10; ++i) {
    os << ',' << int(v[i]);
  }

  os << '}';

  return os;
}

#endif

