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

#ifndef _CYCLOTOMIC_RANDOM_H_
#define _CYCLOTOMIC_RANDOM_H_

#include "common.h"

namespace CyclotomicRandom {

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

};

#endif /* _CYCLOTOMIC_RANDOM_H_ */
