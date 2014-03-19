#ifndef _HIGHER_CYCLO_H_
#define _HIGHER_CYCLO_H_

#include "radial_math.h"


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
    a[0] = v1.x; a[1] = v1.y;
    a[2] = v2.x; a[3] = v2.y;
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
    static const double u = CommonRadial::pi / 14.0;

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

namespace CommonRadial {

  typedef vector<vec4d> vec4dlist;

};

namespace RadialCoprime {

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

};

namespace HeptagonalRadial {

  struct HeptagonalOp {
    typedef vec6s invectype;
    static const double epsilon;

    static double angle(const invectype& a) {
      return a.toPhysicalL7().angle();
    }

    static vec2d toR2(const invectype& a) {
      return a.toPhysicalL7();
    }

    static bool rayTest(const invectype& a, const invectype& b) {
      // convert to direct sum representation
      invectype pa(a.toDirectL7());
      invectype pb(b.toDirectL7());

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
      RadialCoprime::multL7(pa, c);
      RadialCoprime::multL7(pb, d);

      return (c == d);
    }
  };

  typedef VisTest::visibleList<HeptagonalOp> heptagonalVisList;

};

ostream& operator<<(ostream &os, const vec6s& v)
{
  os << '{' << int(v[0]);

  for (uint i = 1; i < 6; ++i) {
    os << ',' << int(v[i]);
  }

  os << '}';

  return os;
}

#endif

