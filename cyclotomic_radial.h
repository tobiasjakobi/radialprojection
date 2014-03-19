#ifndef _CYCLOTOMIC_RADIAL_H_
#define _CYCLOTOMIC_RADIAL_H_

#include "radial_math.h"

#include <cassert>

typedef unsigned int uint;
typedef unsigned short ushort;

namespace OctogonalRadial {

  struct OctogonalOp {
    typedef vec4i invectype;
    static const double epsilon;

    static double angle(const invectype& a) {
      return a.paraProjL8().angle();
    }

    static vec2d toR2(const invectype& a) {
      return a.paraProjL8();
    }

    static bool rayTest(const invectype& a, const invectype& b) {
      // let Z<2> be Z[sqrt[2]]
      // transform into the Z<2>*1 + Z<2>*xi
      // representation (this is a direct sum)
      const invectype pa(a.transL8ToDirect());
      const invectype pb(b.transL8ToDirect());

      // first filter the trivial cases
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

      // pa = z_a + w_a * xi
      // pb = z_b + w_b * xi
      // with z_a, z_b, w_a, w_b elements in Z<2>
      vec2i c, d;

      // now compute:
      // c = z_a * w_b
      // d = z_b * w_a
      RadialCoprime::multZ2(vec2i(pa[0], pa[1]),
                            vec2i(pb[2], pb[3]), c);
      RadialCoprime::multZ2(vec2i(pb[0], pb[1]),
                            vec2i(pa[2], pa[3]), d);

      return (c == d);
    }
  };

  typedef VisTest::visibleList<OctogonalOp> octogonalVisList;

};

namespace DodecagonalRadial {

  struct DodecagonalOp {
    typedef vec4i invectype;
    static const double epsilon;

    static double angle(const invectype& a) {
      return a.paraProjL12().angle();
    }

    static vec2d toR2(const invectype& a) {
      return a.paraProjL12();
    }

    static bool rayTest(const invectype& a, const invectype& b) {
      // let Z<3> be Z[sqrt[3]]
      // transform into the Z<3>*1 + Z<3>*xi
      // representation (this is a direct sum)
      const invectype pa(a.transL12ToDirect());
      const invectype pb(b.transL12ToDirect());

      // first filter the trivial cases
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

      // pa = z_a + w_a * xi
      // pb = z_b + w_b * xi
      // with z_a, z_b, w_a, w_b elements in Z<3>
      vec2i c, d;

      // now compute:
      // c = z_a * w_b
      // d = z_b * w_a
      RadialCoprime::multZ3(vec2i(pa[0], pa[1]),
                            vec2i(pb[2], pb[3]), c);
      RadialCoprime::multZ3(vec2i(pb[0], pb[1]),
                            vec2i(pa[2], pa[3]), d);

      return (c == d);
    }
  };

  typedef VisTest::visibleList<DodecagonalOp> dodecagonalVisList;

};

namespace RhombicPenrose {

  struct RhombicOp {
    typedef vec4i invectype;
    static const double epsilon;

    static double angle(const invectype& a) {
      return a.paraProjL5().angle();
    }

    static vec2d toR2(const invectype& a) {
      return a.paraProjL5();
    }

    static bool rayTest(const invectype& a, const invectype& b) {
      // transform into the Z[tau]*1 + Z[tau]*xi
      // representation (this is a direct sum)
      const invectype pa(a.transL5ToDirect());
      const invectype pb(b.transL5ToDirect());

      // first filter the trivial cases
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

      // pa = z_a + w_a * xi
      // pb = z_b + w_b * xi
      // with z_a, z_b, w_a, w_b elements in Z[tau]
      vec2i c, d;

      // now compute:
      // c = z_a * w_b
      // d = z_b * w_a
      RadialCoprime::multZTau(vec2i(pa[0], pa[1]),
                              vec2i(pb[2], pb[3]), c);
      RadialCoprime::multZTau(vec2i(pb[0], pb[1]),
                              vec2i(pa[2], pa[3]), d);

      return (c == d);
    }
  };

  typedef VisTest::visibleList<RhombicOp> rhombicVisList;

};

#endif

