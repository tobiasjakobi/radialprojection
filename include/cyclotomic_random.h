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

#ifndef _CYCLOTOMIC_RANDOM_H_
#define _CYCLOTOMIC_RANDOM_H_

#include "common.h"

namespace CyclotomicRandom {

  enum random_mode {
    cyclotomic_visrnd,
    cyclotomic_rndvis
  };

  enum processing_mode {
    octagonal_visrnd   = 0, /* octagonal / Ammann-Beenker tiling (L8 lattice) */
    octagonal_rndvis   = 1,
    decagonal_visrnd   = 2, /* decagonal / Tübingen triangle tiling (L5 lattice) */
    decagonal_rndvis   = 3,
    dodecagonal_visrnd = 4, /* dodecagonal / Gähler shield tiling (L12 lattice) */
    dodecagonal_rndvis = 5,
    rhmbpenrose_visrnd = 6, /* rhombic Penrose tiling (L5 with four windows) */
    rhmbpenrose_rndvis = 7,
    processing_mode_end
  };

  bool check_mode(uint mode) {
    return (mode >= processing_mode_end);
  }

  random_mode get_random_mode(uint mode) {
    return (mode % 2 == 0 ? cyclotomic_visrnd : cyclotomic_rndvis);
  }

  void apply_shift(uint mode);

  template <typename T>
  void randomize(const vector<T>& input, vector<T>& output, double prob);

  class RadialFunc {
  public:
    typedef void (*tilingfunc)(const vec4i&, uint, Common::vec4ilist&);
    typedef void (*visfunc)(const vec4i&, uint, bool, Common::vec4ilist&,
                            Common::vec4ilist&);
    typedef void (*extractfunc)(const vec4i&, const Common::vec4ilist&,
                                Common::vec4ilist&);
    typedef void (*radprojfunc)(const Common::vec4ilist&,
                                Common::dlist&, double&);
    typedef uint (*estimatefunc)(uint, bool);

  private:
    tilingfunc projTiling;
    visfunc projTilingVis;
    extractfunc extractVisible;
    radprojfunc radialProj;
    estimatefunc estimateGrowth;

  public:
    RadialFunc(tilingfunc projTiling_, visfunc projTilingVis_,
               extractfunc extractVisible_, radprojfunc radialProj_,
               estimatefunc estimateGrowth_) :
      projTiling(projTiling_), projTilingVis(projTilingVis_),
      extractVisible(extractVisible_), radialProj(radialProj_),
      estimateGrowth(estimateGrowth_) {}

    ~RadialFunc() {}

    RadialFunc(const RadialFunc& rf);

    void call(random_mode mode, uint steps, double prob,
                    Common::dlist& spacings) const;
  };
};

#endif /* _CYCLOTOMIC_RANDOM_H_ */
