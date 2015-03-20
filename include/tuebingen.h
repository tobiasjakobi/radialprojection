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

#ifndef _TUEBINGEN_H_
#define _TUEBINGEN_H_

#include "common.h"
#include "visibility.h"

// Namespace for both the TTT (Tübingen triangle tiling)
// and the related PRT (Penrose-Robinson tiling)
namespace TuebingenTri {

  /* The visibility test behaves the same for both TTT and PRT. */
  struct VisOp {
    typedef Common::vec4slist list_type;
    static const double epsilon;

    static double angle(const vec4s& a) {
      return a.transL5ToR2().angle();
    }

    static vec2d toR2(const vec4s& a) {
      return a.transL5ToR2();
    }

    static bool rayTest(const vec4s& a, const vec4s& b) {
      // transform into the Z[tau]*1 + Z[tau]*xi
      // representation (this is a direct sum)
      const vec4s pa(a.transL5ToDirect());
      const vec4s pb(b.transL5ToDirect());

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
      Coprime::multZTau(vec2i(pa[0], pa[1]),
                        vec2i(pb[2], pb[3]), c);
      Coprime::multZTau(vec2i(pb[0], pb[1]),
                        vec2i(pa[2], pa[3]), d);

      return (c == d);
    }
  };

  typedef VisTest::VisibleList<VisOp> VisList;
  
  class tri;
  typedef vector<tri> trilist;

  class tri {
  private:
    ubyte type;   // type A = 1, type B = 0
    ubyte orient; // default chirality = 1, other chir. = 0
    ubyte rot;    // rotation of the rhomb around ref (angle = (2*Pi/10) * rot)
    ubyte pad;

    vec4s ref;    // reference point of the triangle

   public:

    tri(ubyte t, ubyte o, ubyte rt, 
        const vec4s& rf) : type(t % 2), orient(o % 2), rot(rt % 10), ref(rf) {}

    tri(const tri& t) : type(t.type), orient(t.orient),
                        rot(t.rot), ref(t.ref) {}

    ubyte getType() const {return type;}
    ubyte getOrient() const {return orient;}
    uint getRot() const {return uint(rot);}
    const vec4s& getRef() const {return ref;}

    // Apply TTT-inflation to triangle and
    // add resulting triangles to list
    void inflateTTT(trilist& list) const;

    // Apply PRT-inflation to triangle and
    // add resulting triangles to list
    void inflatePRT(trilist& list) const;

    // Store the three vertices of the triangle into list
    void getVertices(vec4s* list) const;

  private:
    void inflateTTTa(trilist& list) const;
    void inflateTTTb(trilist& list) const;

    void inflatePRTa(trilist& list) const;
    void inflatePRTb(trilist& list) const;

    void getVerticesA(vec4s* list) const;
    void getVerticesB(vec4s* list) const;

  };

  ostream& operator<<(ostream &os, const tri& t);

};

#endif

