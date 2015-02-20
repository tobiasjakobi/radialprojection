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

#include "cyclotomic_random.h"

#include "cyclotomic_octagonal.h"
#include "cyclotomic_decagonal.h"
#include "cyclotomic_dodecagonal.h"
#include "cyclotomic_rhombic_penrose.h"

#include <sstream>

template <typename T>
void CyclotomicRandom::randomize(const vector<T>& input,
                          vector<T>& output, double prob) {
  Common::srandExt();
  const double rndNorm = 1.0 / double(RAND_MAX);

  for (typename vector<T>::const_iterator i = input.begin();
       i != input.end(); ++i) {
    const double p = double(rand()) * rndNorm;
    if (p >= prob) output.push_back(*i);
  }
}

int main(int argc, char* argv[]) {
  using namespace CyclotomicRandom;

  const vec4i init(0, 0, 0, 0);

  stringstream parser;

  uint mode = 0;
  uint steps = 100;

  Common::vec4ilist tiling, visible;

  if (argc >= 2) {
    parser.str(argv[1]);
    parser.clear();
    parser >> mode;

    if (argc >= 3) {
      parser.str(argv[2]);
      parser.clear();
      parser >> steps;
    }
  }

  if (check_mode(mode)) {
    cerr << "error: unknown mode (" << mode <<  ") selected.\n";
    return 1;
  }

  switch (mode) {
    case octagonal_visrnd:
      Octogonal::projTilingVisLocal(init, steps, false, tiling, visible);
      tiling.clear(); /* original tiling vertices are not used */
      randomize(visible, tiling, 0.5);
      cerr << "info: after randomization " << tiling.size()
           << " vertices remain\n";
    break;

    case octagonal_rndvis:
      /* TODO: implement */
    break;

    case decagonal_visrnd:
      Decagonal::projTilingVisLocal(init, steps, tiling, visible);
      tiling.clear();
      randomize(visible, tiling, 0.5);
      cerr << "info: after randomization " << tiling.size()
           << " vertices remain\n";
    break;

    case decagonal_rndvis:
      /* TODO: implement */
    break;

    case dodecagonal_visrnd:
      Dodecagonal::projTilingVisLocal(init, steps, tiling, visible);
      tiling.clear();
      randomize(visible, tiling, 0.5);
      cerr << "info: after randomization " << tiling.size()
           << " vertices remain\n";
    break;

    case dodecagonal_rndvis:
      /* TODO: implement */
    break;

    case rhmbpenrose_visrnd:
      RhombicPenrose::projTilingVis(init, init, steps, false, tiling, visible);
      tiling.clear();
      randomize(visible, tiling, 0.5);
      cerr << "info: after randomization " << tiling.size()
           << " vertices remain\n";
    break;

    case rhmbpenrose_rndvis:
      /* TODO: implement */
    break;

    default:
      assert(false);
    break;
  }

  return 0;
}
