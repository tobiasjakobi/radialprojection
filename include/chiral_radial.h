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

#ifndef _CHIRAL_RADIAL_H_
#define _CHIRAL_RADIAL_H_

#include "common.h"
#include "visibility.h"

namespace ChiralLB {

  struct VisOp {
    typedef Common::vec4slist list_type;
    static const double epsilon;

    static double angle(const vec4s& a) {
      return a.transL10ToR2().angle();
    }

    static vec2d toR2(const vec4s& a) {
      return a.transL10ToR2();
    }

    static bool rayTest(const vec4s& a, const vec4s& b) {
      // transform into the Z[tau]*1 + Z[tau]*xi
      // representation (this is a direct sum)
      const vec4s pa(a.transL10ToDirect());
      const vec4s pb(b.transL10ToDirect());

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

  class rhomb;
  typedef vector<rhomb> rhomblist;

  class rhomb {
  public:
    ushort type; // type A = 0, type B = 1
    ushort rot;  // rotation of the rhomb around ref (angle = psi/2 * rot)

    vec8s ref;   // reference point

    rhomb() : type(0), rot(0), ref() {}

    rhomb(ushort t, ushort r1,
          const vec8s& r2) : type(t % 2), rot(r1 % 20), ref(r2) {}

    rhomb(const rhomb& r) : type(r.type),
                            rot(r.rot), ref(r.ref) {}

    // Inflate rhomb and add resulting rhombs to list
    void inflate(rhomblist& list) const;

    // Store the four vertices of the rhomb into list
    void getVertices(vec8s* list) const;

  private:
    // Construct "midpoint" of the scaled rhomb
    vec8s midpoint() const;
  };

  ostream& operator<<(ostream &os, const rhomb& r);

};

namespace Chair2D {

  class chairL;
  typedef vector<chairL> llist;

  class chairL {
  public:
    uint rot;
    vec2s ref;

    chairL() : rot(0), ref() {}

    chairL(uint r, const vec2s& v) : rot(r), ref(v) {}

    chairL(const chairL& l) : rot(l.rot),
                              ref(l.ref) {}

    void inflate(llist& list) const;
    void getVertices(vec2s* list) const;

  };

  class vismap {
  private:
    typedef vector<bool> blist;

    int gridsize;
    vector<blist>* bmap;

    bool access(const vec2s& v) const {
      const int x = v[0];
      const int y = v[1];

      assert(-gridsize <= x && x <= gridsize &&
             -gridsize <= y && y <= gridsize);
      return (*bmap)[x + gridsize][y + gridsize];
    }

  public:

    // constructor
    vismap(const llist& patch, uint steps) {
      gridsize = 2 * Common::ipower(2, steps);
      bmap = new vector<blist>(2 * gridsize + 1, blist(2 * gridsize + 1, false));

      for (llist::const_iterator i = patch.begin(); i != patch.end(); ++i) {
        vec2s temp[6];
        i->getVertices(temp);

        for (uint j = 0; j < 6; ++j) {
          const vec2s current(temp[j]);

          assert(-gridsize <= current[0] && current[0] <= gridsize &&
                 -gridsize <= current[1] && current[1] <= gridsize);
          (*bmap)[current[0] + gridsize][current[1] + gridsize] = true;
        }
      }
    }

    // destructor
    ~vismap() {
      delete bmap;
    }

    // copy-constructor
    vismap(const vismap& m) {
      gridsize = m.gridsize;
      bmap = new vector<blist>(*m.bmap);
    }

    bool isVisible(const vec2s& v) const {
      if (v[0] == 0 && v[1] == 0) return false;

      const int agcd = Coprime::gcdZFast(abs(v[0]), abs(v[1]));
      if (agcd == 1) return true;

      /* A single check for coprime coordinates isn't correct in this case, since the   *
       * corresponding primitive point doesn't have to be part of the tiling vertices   *
       * (this is different for the entire Z2 lattice).                                 *
       * A point can still be visible even though its coordinates aren't coprime!       *
       *                                                                                *
       * Indeed this does make a difference to the radial projection if this correct    *
       * method is used (instead of the "naive" (and wrong) gcd-only test).             */

      const vec2s primitive(v[0] / agcd, v[1] / agcd);

      for (int k = 1; k < agcd; ++k) {
        if (access(primitive * k)) return false;
      }

      return true;
    }

  };

  ostream& operator<<(ostream &os, const chairL& l);

};

#endif

