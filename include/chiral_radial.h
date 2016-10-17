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

  typedef double (&clipfunc)(uint);

  template <typename T>
  class chairL {
  private:
    typedef chairL<T> item_type;

    uint rot;
    T ref;

  public:
    typedef T coord_type;
    typedef vector<item_type> list_type;

    chairL() : rot(0), ref() {}
    chairL(uint r, const T& v) : rot(r), ref(v) {}
    chairL(const item_type& l) : rot(l.rot), ref(l.ref) {}

    uint getRot() const { return rot; }
    const T& getRef() const { return ref; }

    void inflate(list_type& list) const;
    void getVertices(T* list) const;

    /*
     * Check if at least one vertex of the tile is inside the circle
     * of radius 'r' around the origin.
     */
    bool clip(double r) const;

  };

  template <typename T>
  class VisibilityMap {
  public:
    typedef typename T::coord_type coord_type;
    typedef typename T::list_type list_type;
    typedef typename list_type::const_iterator iter_type;

  private:
    typedef vector<bool> boolvec;
    typedef boolvec::size_type bvsz;

    bvsz range, rowlen, size;
    boolvec vmap;

    bool access(const coord_type& v) const {
      assert((v[0] + int(range)) >= 0 && (v[1] + int(range)) >= 0);

      const bvsz pos = (v[1] + range) * rowlen + (v[0] + range);

      assert(pos < size);

      return vmap[pos];
    }

    void set(const coord_type& v) {
      assert((v[0] + int(range)) >= 0 && (v[1] + int(range)) >= 0);

      const bvsz pos = (v[1] + range) * rowlen + (v[0] + range);

      assert(pos < size);

      vmap[pos] = true;
    }

  public:
    // constructor
    VisibilityMap(const list_type& patch, uint steps) {
      range = 2 * Common::ipower(2, steps);

      assert(range <= bvsz(numeric_limits<int>::max()));

      rowlen = 2 * range + 1;
      size = rowlen * rowlen;

      vmap.resize(size, false);

      for (iter_type i = patch.begin(); i != patch.end(); ++i) {
        coord_type temp[6];
        i->getVertices(temp);

        for (uint j = 0; j < 6; ++j)
          set(temp[j]);
      }
    }

    // destructor
    ~VisibilityMap() {}

    // copy-constructor
    VisibilityMap(const VisibilityMap& vm) : range(vm.range), rowlen(vm.rowlen),
      size(vm.size), vmap(vm.vmap) {}

    bool isVisible(const coord_type& v) const {
      if (v[0] == 0 && v[1] == 0)
        return false;

      const int agcd = Coprime::gcdZFast(abs(v[0]), abs(v[1]));
      if (agcd == 1)
        return true;

      /* A single check for coprime coordinates isn't correct in this case, since the   *
       * corresponding primitive point doesn't have to be part of the tiling vertices   *
       * (this is different for the entire Z2 lattice).                                 *
       * A point can still be visible even though its coordinates aren't coprime!       *
       *                                                                                *
       * Indeed this does make a difference to the radial projection if this correct    *
       * method is used (instead of the "naive" (and wrong) gcd-only test).             */

      const coord_type primitive(v[0] / agcd, v[1] / agcd);

      for (int k = 1; k < agcd; ++k) {
        if (access(primitive * k))
          return false;
      }

      return true;
    }

  };

  template <typename T>
  ostream& operator<<(ostream &os, const chairL<T>& l);

};

#endif

